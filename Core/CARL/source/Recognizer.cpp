/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <carl/Recognizer.h>

#include "Descriptor.h"
#include "DynamicTimeWarping.h"
#include "SessionImpl.h"

#include <atomic>

namespace carl::action
{
    namespace
    {
        // Slightly less naive holistic resampling. Still linear, though.
        std::vector<InputSample> resample(gsl::span<const InputSample> samples, NumberT startTimestamp, NumberT endTimestamp, NumberT frameDuration)
        {
            size_t samplesIdx = 0;

            std::vector<InputSample> newSamples{};
            size_t newSamplesCount = std::max<size_t>(static_cast<size_t>(std::ceil((endTimestamp - startTimestamp) / frameDuration)), 1);
            newSamples.reserve(newSamplesCount);

            for (size_t idx = 0; idx < newSamplesCount; ++idx)
            {
                NumberT timestamp = startTimestamp + frameDuration * idx;

                while (samplesIdx < samples.size() && samples[samplesIdx].Timestamp < timestamp)
                {
                    ++samplesIdx;
                }

                if (samplesIdx == 0)
                {
                    newSamples.emplace_back(samples.front());
                }
                else if (samplesIdx == samples.size())
                {
                    newSamples.emplace_back(samples.back());
                }
                else
                {
                    auto& earlierSample = samples[samplesIdx - 1];
                    auto& laterSample = samples[samplesIdx];
                    NumberT t = (timestamp - earlierSample.Timestamp) / (laterSample.Timestamp - earlierSample.Timestamp);
                    newSamples.emplace_back(InputSample::Lerp(earlierSample, laterSample, t));
                }
            }

            return newSamples;
        }

        std::vector<InputSample> resample(const action::Example& original, NumberT frameDuration)
        {
            return resample(original.getRecording().getSamples(), original.getStartTimestamp(), original.getEndTimestamp(), frameDuration);
        }

        std::vector<action::Example> expandExamples(gsl::span<const action::Example> examples)
        {
            std::vector<action::Example> expandedExamples{};

            for (const auto& example : examples)
            {
                const auto& samples = example.getRecording().getSamples();
                const auto startTimestamp = example.getStartTimestamp();
                const auto endTimestamp = example.getEndTimestamp();
                
                size_t idx = 0;
                while (idx < samples.size() && samples[idx].Timestamp < startTimestamp)
                {
                    ++idx;
                }
                do
                {
                    InProgressRecording inProgressRecording{};
                    inProgressRecording.addSample(samples[idx]);
                    expandedExamples.emplace_back(std::move(inProgressRecording), samples[idx].Timestamp, samples[idx].Timestamp);
                    ++idx;
                } while (idx < samples.size() && samples[idx].Timestamp < endTimestamp);
            }

            return expandedExamples;
        }

        template<typename T>
        class has_get_timestamp
        {
            template<typename InternalT>
            static constexpr auto test(uint8_t) -> decltype(std::declval<InternalT>().getTimestamp(), bool())
            {
                return true;
            }

            template<typename...>
            static constexpr bool test(uint32_t)
            {
                return false;
            }

        public:
            static constexpr inline bool value{ test<T>(uint8_t{}) };
        };
    }

    class action::Recognizer::Impl
    {
        friend class action::Recognizer;

    public:
        Impl(Session& session, NumberT sensitivity)
            : m_sessionImpl{ Session::Impl::getFromSession(session) }
            , m_sensitivity{ sensitivity }
        {
        }

        virtual ~Impl() = default;

    protected:
        virtual Example createAutoTrimmedExample(const Recording&) const = 0;
        virtual void analyzeRecording(const Recording& recording, std::ostream& output) const = 0;

    protected:
        Session::Impl& m_sessionImpl;
        arcana::weak_table<Signal<bool>::HandlerT> m_whenRecognitionChangedHandlers{};
        std::atomic<NumberT> m_currentScore{};
        std::unique_ptr<Recording> m_canonicalRecording{};
        NumberT m_sensitivity{};
    };

    namespace
    {
        template <typename DescriptorT>
        class RecognizerImpl final : public action::Recognizer::Impl
        {
        public:
            RecognizerImpl(Session& session, gsl::span<const action::Example> examples, gsl::span<const action::Example> counterexamples, NumberT sensitivity)
                : Recognizer::Impl{ session, sensitivity }
                , m_ticket{ Session::Impl::getFromSession(session).addHandler<DescriptorT>(
                    [this](gsl::span<const DescriptorT> sequence, size_t newDescriptorsCount) { handleSequence(sequence, newDescriptorsCount); }) }
                , m_distanceFunction{ [this](const auto& a, const auto& a0, const auto& b, const auto& b0) { return DescriptorT::Distance(a, a0, b, b0, m_tuning); } }
            {
                m_tuning = DescriptorT::DEFAULT_TUNING;
                initializeTemplates(examples, counterexamples);
                calculateTuning(examples);
                // TODO: Figure out why the tuning resampling negatively impacts recognition, then substitute the following for the above
                // calculateTuning(examples);
                // initializeTemplates(examples, counterexamples);
                createScoringFunction();
                createCanonicalRecording(examples);
            }

        protected:
            Example createAutoTrimmedExample(const Recording& recording) const override
            {
                std::vector<descriptor::TimestampedDescriptor<DescriptorT>> timestampedSequence{};
                auto samples = recording.getSamples();
                auto mostRecentSample = samples.front();
                for (const auto& sample : recording.getSamples())
                {
                    descriptor::extendSequence<descriptor::TimestampedDescriptor<DescriptorT>>(sample, timestampedSequence, mostRecentSample, m_tuning);
                }

                std::vector<DescriptorT> sequence{};
                sequence.reserve(timestampedSequence.size());
                for (const auto& desc : timestampedSequence)
                {
                    sequence.push_back(desc.getUnderlyingDescriptor());
                }

                auto bestMatchResult = DynamicTimeWarping::Match<const DescriptorT>(sequence, m_templates[0], m_distanceFunction);
                for (size_t idx = 1; idx < m_templates.size(); ++idx)
                {
                    auto matchResult = DynamicTimeWarping::Match<const DescriptorT>(sequence, m_templates[idx], m_distanceFunction);
                    if (matchResult.MatchCost < bestMatchResult.MatchCost)
                    {
                        bestMatchResult = matchResult;
                    }
                }

                auto firstDescriptorT = timestampedSequence[bestMatchResult.ImageStartIdx].getTimestamp();
                auto lastDescriptorT = timestampedSequence[bestMatchResult.ImageStartIdx + bestMatchResult.ImageSize - 1].getTimestamp();
                constexpr auto epsilonT = std::numeric_limits<double>::epsilon();

                auto startT = std::numeric_limits<NumberT>::lowest() + 2 * epsilonT;
                auto endT = std::numeric_limits<NumberT>::max() - 2 * epsilonT;
                for (const auto& sample : samples)
                {
                    if (sample.Timestamp < firstDescriptorT && sample.Timestamp > startT)
                    {
                        startT = sample.Timestamp;
                    }

                    if (sample.Timestamp > lastDescriptorT && sample.Timestamp < endT)
                    {
                        endT = sample.Timestamp;
                    }
                }

                return{ recording, startT - epsilonT, endT + epsilonT};
            }

            void analyzeRecording(const Recording& recording, std::ostream& output) const override
            {
                auto inspector = recording.getInspector();
                auto targetSequence = descriptor::createDescriptorSequenceFromRecording<DescriptorT>(recording, inspector.startTimestamp(), inspector.endTimestamp(), DescriptorT::DEFAULT_TUNING);

                const auto outputAnalysis = [&output, &targetSequence](auto analysis, std::string label, const auto& querySequence, size_t idx) {
                    for (const auto& [name, identicality, tuning, rows] : analysis)
                    {
                        output << idx << ",\"" << label << "\",\"" << name << "\"," << identicality << "," << tuning << "," << std::endl;

                        if constexpr (has_get_timestamp<DescriptorT>::value)
                        {
                            output << "\"TIMESTAMPS\",";
                            for (const auto& descriptor : targetSequence)
                            {
                                output << descriptor.getTimestamp() << ",";
                            }
                            output << std::endl;
                        }

                        for (size_t idx = 0; idx < rows.size(); ++idx)
                        {
                            if constexpr (has_get_timestamp<DescriptorT>::value)
                            {
                                output << querySequence[idx].getTimestamp() << ",";
                            }

                            const auto& row = rows[idx];
                            for (const auto& entry : row)
                            {
                                output << entry.MatchCost / entry.Connections << ",";
                            }

                            output << std::endl;
                        }

                        output << std::endl;
                    }
                };

                for (size_t idx = 0; idx < m_templates.size(); ++idx)
                {
                    const auto& t = m_templates[idx];
                    outputAnalysis(DescriptorT::template Analyze<false>(targetSequence, t, m_tuning), "Raw Distance", t, idx);
                    outputAnalysis(DescriptorT::template Analyze<true>(targetSequence, t, m_tuning), "Normalized Distance", t, idx);
                }
            }

        private:
            typename Signal<gsl::span<const DescriptorT>, size_t>::TicketT m_ticket;
            std::vector<std::vector<DescriptorT>> m_templates{};
            std::vector<std::vector<DescriptorT>> m_countertemplates{};
            size_t m_trimmedSequenceLength{};
            std::array<NumberT, DescriptorT::DEFAULT_TUNING.size()> m_tuning{};
            std::function<NumberT(const DescriptorT&, const DescriptorT&, const DescriptorT&, const DescriptorT&)> m_distanceFunction{};
            std::function<NumberT(NumberT)> m_scoringFunction{};
            bool m_recognition{ false };

            void initializeTemplates(gsl::span<const action::Example> examples, gsl::span<const action::Example> counterexamples)
            {
                m_templates = descriptor::createDescriptorSequencesFromExamples<DescriptorT>(examples, m_tuning);
                m_countertemplates = descriptor::createDescriptorSequencesFromExamples<DescriptorT>(counterexamples, m_tuning);

                // TODO: Parameterize this calculation, instead of hard-coding 5/4ths?
                for (const auto& t : m_templates)
                {
                    m_trimmedSequenceLength = std::max(m_trimmedSequenceLength, (5 * t.size()) / 4);
                }
                for (const auto& ct : m_countertemplates)
                {
                    m_trimmedSequenceLength = std::max(m_trimmedSequenceLength, (5 * ct.size()) / 4);
                }

                m_sessionImpl.supportSequenceOfLength<DescriptorT>(2 * m_trimmedSequenceLength);
            }

            void calculateTuning(gsl::span<const Example> examples)
            {
                m_tuning = DescriptorT::CalculateTuning(examples);
            }

            void createUnitScoringFunction()
            {
                m_scoringFunction = [this](NumberT distance)
                {
                    distance /= DescriptorT::DEFAULT_TUNING.size();
                    return std::max(1. - std::pow(distance / (3.16228 * m_sensitivity), 2.), 0.);
                };
            }

            void createScoringFunction()
            {
                createUnitScoringFunction();
            }

            /// <summary>
            /// Creates a new recording based on the "centroid" example. Used for things like tutorials.
            /// </summary>
            /// <param name="examples">The examples from which this recognizer was created. MUST be in the same order as when the templates were initialized.</param>
            void createCanonicalRecording(gsl::span<const Example> examples)
            {
                // Select which recording is the "centroid"
                // TODO: This does not currently take counterexamples into account very well since they won't apply AT the site of the score.
                NumberT maxScore = calculateScore(m_templates[0], m_templates[0].size());
                size_t maxIdx = 0;
                for (size_t idx = 1; idx < m_templates.size(); ++idx)
                {
                    NumberT score = calculateScore(m_templates[idx], m_templates[idx].size());
                    if (score > maxScore)
                    {
                        maxScore = score;
                        maxIdx = idx;
                    }
                }

                // Create a new recording based on the "centroid" example, discarding unneeded information
                const auto& canonicalExample = examples[maxIdx];
                InProgressRecording recording{};
                for (const auto& sample : canonicalExample.getRecording().getSamples())
                {
                    if (sample.Timestamp >= canonicalExample.getStartTimestamp() && sample.Timestamp <= canonicalExample.getEndTimestamp())
                    {
                        auto sampleCopy{ sample };
                        /* TODO: Is this necessary?
                        if constexpr (DescriptorT::HANDEDNESS == descriptor::Handedness::LeftHanded)
                        {
                            sampleCopy.RightWristPose.reset();
                            sampleCopy.RightHandJointPoses.reset();
                        }
                        else if constexpr (DescriptorT::HANDEDNESS == descriptor::Handedness::RightHanded)
                        {
                            sampleCopy.LeftHandJointPoses.reset();
                            sampleCopy.LeftHandJointPoses.reset();
                        }
                        */
                        recording.addSample(std::move(sampleCopy));
                    }
                }
                m_canonicalRecording = std::make_unique<Recording>(std::move(recording));
            }

            NumberT calculateMatchDistance(
                gsl::span<const DescriptorT> target,
                gsl::span<const DescriptorT> query,
                size_t minimumImageEndIdx = 0) const
            {
                auto result = DynamicTimeWarping::Match(target, query, m_distanceFunction, minimumImageEndIdx);
                return result.MatchCost / result.Connections;
            }

            NumberT calculateScore(gsl::span<const DescriptorT> sequence, size_t newDescriptorsCount) const
            {
                // Early-out if the provided sequence is too short.
                if (sequence.size() < m_trimmedSequenceLength)
                {
                    return 0;
                }

                gsl::span<const DescriptorT> trimmedSequence{ &sequence[sequence.size() - m_trimmedSequenceLength], m_trimmedSequenceLength };

                // Early-out if no template's end is scored as appearing in the sequence
                NumberT score = std::numeric_limits<NumberT>::lowest();
                for (const auto& t : m_templates)
                {
                    auto endSpan = gsl::make_span<const DescriptorT>(&t.back(), 1);
                    auto distance = calculateMatchDistance(trimmedSequence, endSpan);
                    score = std::max(m_scoringFunction(distance), score);
                }
                if (score < std::numeric_limits<NumberT>::epsilon())
                {
                    return 0;
                }

                const size_t firstNewDescriptorIdx = newDescriptorsCount <= trimmedSequence.size() ? trimmedSequence.size() - newDescriptorsCount : 0;

                // Calculate the base score based on proximity to templates
                score = std::numeric_limits<NumberT>::lowest();
                NumberT minDistance = std::numeric_limits<NumberT>::max();
                for (const auto& t : m_templates)
                {
                    NumberT distance = calculateMatchDistance(trimmedSequence, t, firstNewDescriptorIdx);
                    score = std::max(m_scoringFunction(distance), score);
                    minDistance = std::min(distance, minDistance);
                }

                // Clamp score to [0, 1] range.
                score = std::clamp<NumberT>(score, 0, 1);

                // Penalize score based on relative proximity to countertemplates
                NumberT minCounterDistance = std::numeric_limits<NumberT>::max();
                for (const auto& t : m_countertemplates)
                {
                    NumberT distance = calculateMatchDistance(trimmedSequence, t, firstNewDescriptorIdx);
                    minCounterDistance = std::min(distance, minCounterDistance);
                }
                if (minDistance >= minCounterDistance)
                {
                    score *= minCounterDistance / minDistance;
                }

                return score;
            }

            void handleSequence(gsl::span<const DescriptorT> sequence, size_t newDescriptorsCount)
            {
                auto score = calculateScore(sequence, newDescriptorsCount);
                if (!m_recognition && score > 0.7)
                {
                    m_recognition = true;
                    m_whenRecognitionChangedHandlers.apply_to_all([this](auto& callable) { callable(m_recognition); });
                }
                else if (m_recognition)
                {
                    m_recognition = false;
                    m_whenRecognitionChangedHandlers.apply_to_all([this](auto& callable) { callable(m_recognition); });
                }
                m_currentScore.store(score);
            }
        };

        std::unique_ptr<action::Recognizer::Impl> createImpl(
            Session& session,
            const action::Definition& definition)
        {
            std::vector<action::Example> examples{};
            std::vector<action::Example> counterexamples{};
            switch (definition.getDescriptorType())
            {
            case action::Definition::ActionType::LeftHandPose:
                examples = expandExamples(definition.getExamples());
                counterexamples = expandExamples(definition.getCounterexamples());
                return std::make_unique<RecognizerImpl<descriptor::HandPose<descriptor::Handedness::LeftHanded>>>(session, examples, counterexamples, definition.DefaultSensitivity);
            case action::Definition::ActionType::LeftHandGesture:
                return std::make_unique<RecognizerImpl<descriptor::HandGesture<descriptor::Handedness::LeftHanded>>>(session, definition.getExamples(), definition.getCounterexamples(), definition.DefaultSensitivity);
            case action::Definition::ActionType::RightHandPose:
                examples = expandExamples(definition.getExamples());
                counterexamples = expandExamples(definition.getCounterexamples());
                return std::make_unique<RecognizerImpl<descriptor::HandPose<descriptor::Handedness::RightHanded>>>(session, examples, counterexamples, definition.DefaultSensitivity);
            case action::Definition::ActionType::RightHandGesture:
                return std::make_unique<RecognizerImpl<descriptor::HandGesture<descriptor::Handedness::RightHanded>>>(session, definition.getExamples(), definition.getCounterexamples(), definition.DefaultSensitivity);
            case action::Definition::ActionType::TwoHandGesture:
                return std::make_unique<RecognizerImpl<descriptor::TwoHandGesture>>(session, definition.getExamples(), definition.getCounterexamples(), definition.DefaultSensitivity);
            case action::Definition::ActionType::LeftControllerGesture:
                return std::make_unique<RecognizerImpl<descriptor::ControllerGesture<descriptor::Handedness::LeftHanded>>>(session, definition.getExamples(), definition.getCounterexamples(), definition.DefaultSensitivity);
            case action::Definition::ActionType::RightControllerGesture:
                return std::make_unique<RecognizerImpl<descriptor::ControllerGesture<descriptor::Handedness::RightHanded>>>(session, definition.getExamples(), definition.getCounterexamples(), definition.DefaultSensitivity);
            case action::Definition::ActionType::TwoControllerGesture:
                return std::make_unique<RecognizerImpl<descriptor::TwoControllerGesture>>(session, definition.getExamples(), definition.getCounterexamples(), definition.DefaultSensitivity);
            case action::Definition::ActionType::LeftWristTrajectory:
                return std::make_unique<RecognizerImpl<descriptor::WristTrajectory<descriptor::Handedness::LeftHanded>>>(session, definition.getExamples(), definition.getCounterexamples(), definition.DefaultSensitivity);
            case action::Definition::ActionType::RightWristTrajectory:
                return std::make_unique<RecognizerImpl<descriptor::WristTrajectory<descriptor::Handedness::RightHanded>>>(session, definition.getExamples(), definition.getCounterexamples(), definition.DefaultSensitivity);
            case action::Definition::ActionType::LeftHandShape:
                examples = expandExamples(definition.getExamples());
                counterexamples = expandExamples(definition.getCounterexamples());
                return std::make_unique<RecognizerImpl<descriptor::HandShape<descriptor::Handedness::LeftHanded>>>(session, examples, counterexamples, definition.DefaultSensitivity);
            case action::Definition::ActionType::RightHandShape:
                examples = expandExamples(definition.getExamples());
                counterexamples = expandExamples(definition.getCounterexamples());
                return std::make_unique<RecognizerImpl<descriptor::HandShape<descriptor::Handedness::RightHanded>>>(session, examples, counterexamples, definition.DefaultSensitivity);
            default:
                throw std::runtime_error{ "Unknown definition type" };
            }
        }
    }

    action::Recognizer::Recognizer(Session& session, const action::Definition& definition)
        : m_impl{ createImpl(session, definition) }
        , whenRecognitionChangedSignal{ m_impl->m_whenRecognitionChangedHandlers }
    {
    }

    action::Recognizer::~Recognizer()
    {
    }

    double action::Recognizer::currentScore()
    {
        return static_cast<NumberT>(m_impl->m_currentScore);
    }

    RecordingInspector Recognizer::getCanonicalRecordingInspector() const
    {
        return m_impl->m_canonicalRecording->getInspector();
    }

    void Recognizer::setSensitivity(double sensitivity)
    {
        m_impl->m_sensitivity = static_cast<NumberT>(sensitivity);
    }

    Example Recognizer::createAutoTrimmedExample(const Recording& recording) const
    {
        return m_impl->createAutoTrimmedExample(recording);
    }

    void Recognizer::analyzeRecording(const Recording& recording, std::ostream& output) const
    {
        return m_impl->analyzeRecording(recording, output);
    }
}

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
        // This is an extremely naive hack for dealing with the problem of inifinitessimal differences. This
        // problem most prominently emerges when comparing single-frame actions, which by their nature will
        // contain no motion and will thus will have no difference in translation space. This will cause a 0
        // average distance, which in turn will put infinity into the tunings, causing all sorts of
        // problems. This hack is _sort of_ on the path to one of the correct approaches to solving this,
        // but not really. Two different approaches should be solved to try this. First, there SHOULD be
        // average distance epsilons (essentially default [or possibly minimum] average distances), but they
        // should be hand-tuned per descriptor tuning dimension; this is because descriptors themselves are
        // the only place with a reasonable idea of expected variance (this is sort of a way to allow
        // descriptors to self-normalize their distance space). Second, an alternate creation path should be
        // introduced for static pose gestures where the system dices up a longer segment into chunks.
        constexpr NumberT kAverageDistanceEpsilon{ 0.000000001 };

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

        std::vector<action::Example> expandExamples(gsl::span<const action::Example> examples, NumberT durationT = 0.2)
        {
            std::vector<action::Example> expandedExamples{};

            for (const auto& example : examples)
            {
                if (example.getEndTimestamp() - example.getStartTimestamp() < durationT)
                {
                    expandedExamples.emplace_back(example);
                }
                else
                {
                    for (auto endT = example.getStartTimestamp() + durationT; endT < example.getEndTimestamp(); endT += durationT)
                    {
                        expandedExamples.emplace_back(example.getRecording(), endT - durationT, endT);
                    }
                }
            }

            return expandedExamples;
        }
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
        class RecognizerImpl : public action::Recognizer::Impl
        {
        public:
            RecognizerImpl(Session& session, gsl::span<const action::Example> examples, gsl::span<const action::Example> counterexamples, NumberT sensitivity)
                : Recognizer::Impl{ session, sensitivity }
                , m_ticket{ Session::Impl::getFromSession(session).addHandler<DescriptorT>(
                    [this](gsl::span<const DescriptorT> sequence) { handleSequence(sequence); }) }
            {
                initializeTemplates(examples, counterexamples);
                calculateTuning();
                createScoringFunction();
                createCanonicalRecording(examples);
            }

        protected:
            // TODO: This has not been reworked or retested sufficiently since paradigms shifted underneath it! The underlying 
            // logic is probably most still sound, but the implementation needs to be fixed and tested now that samples and 
            // descriptors can differ in count.
            Example createAutoTrimmedExample(const Recording& recording) const override
            {
                using DescT = descriptor::TimestampedDescriptor<DescriptorT>;
                using SignalT = Signal<const InputSample&>;
                arcana::weak_table<SignalT::HandlerT> inputSamplesHandlers{};
                SignalT inputSampleSignal{ inputSamplesHandlers };
                typename DescriptorSequence<DescT>::Provider descriptorSequenceProvider{ inputSampleSignal };
                descriptorSequenceProvider.supportSequenceOfLength(2 * m_trimmedSequenceLength);
                Signal<gsl::span<const DescT>>& descriptorSignal{ descriptorSequenceProvider };

                auto samples = recording.getSamples();
                size_t idx = 0;

                auto maxScore = std::numeric_limits<NumberT>::lowest();
                size_t maxScoreIdx = 0;

                std::vector<DescriptorT> descriptorSequence{};
                std::vector<DescT> maxScoreTrimmedSequence{};

                using OptionalTicketT = std::optional<typename Signal<gsl::span<const DescT>>::TicketT>;
                OptionalTicketT maxScoreDescriptorHandlerTicket{ descriptorSignal.addHandler(
                    [this, &idx, &samples, &maxScore, &maxScoreIdx, &maxScoreTrimmedSequence, &descriptorSequence](gsl::span<const DescT> sequence) {
                        if (sequence.size() <= m_trimmedSequenceLength)
                        {
                            return;
                        }

                        descriptorSequence.clear();
                        for (const auto& element : sequence)
                        {
                            descriptorSequence.push_back(element.getUnderlyingDescriptor());
                        }
                        auto score = calculateScore(descriptorSequence);
                        if (score > maxScore)
                        {
                            maxScore = score;
                            maxScoreIdx = idx;

                            gsl::span<const DescT> trimmedSequence{
                                &sequence[sequence.size() - m_trimmedSequenceLength], m_trimmedSequenceLength};

                            // Note that this logic depends on the fact that descriptors are trivially copyable.
                            static_assert(std::is_trivially_copyable<DescT>::value);
                            maxScoreTrimmedSequence.resize(trimmedSequence.size());
                            std::memcpy(maxScoreTrimmedSequence.data(), trimmedSequence.data(), sizeof(DescT) * trimmedSequence.size());
                        }
                    }) };

                for (; idx < samples.size(); ++idx)
                {
                    auto& sample = samples[idx];
                    inputSamplesHandlers.apply_to_all([&sample](auto& callable) mutable { callable(sample); });
                }

                maxScoreDescriptorHandlerTicket.reset();

                auto startT = samples[0].Timestamp;
                auto endT = samples[samples.size() - 1].Timestamp;
                auto distance = std::numeric_limits<NumberT>::max();

                descriptorSequence.clear();
                for (const auto& element : maxScoreTrimmedSequence)
                {
                    descriptorSequence.push_back(element.getUnderlyingDescriptor());
                }

                for (const auto& t : m_templates)
                {
                    auto [normalizedDistance, imageSize] = calculateSequenceDistance(descriptorSequence, t);
                    if (normalizedDistance < distance)
                    {
                        // TODO: Rework this to use descriptor timestamps!
                        startT = samples[maxScoreIdx - imageSize].Timestamp;
                        endT = samples[maxScoreIdx].Timestamp;
                        distance = normalizedDistance;
                    }
                }

                return{ recording, startT, endT };
            }

            void analyzeRecording(const Recording& recording, std::ostream& output) const override
            {
                using SignalT = Signal<const InputSample&>;
                arcana::weak_table<SignalT::HandlerT> inputSamplesHandlers{};
                SignalT inputSampleSignal{ inputSamplesHandlers };
                typename DescriptorSequence<DescriptorT>::Provider descriptorSequenceProvider{ inputSampleSignal };
                descriptorSequenceProvider.supportSequenceOfLength(2 * m_trimmedSequenceLength);
                Signal<gsl::span<const DescriptorT>>& descriptorSignal{ descriptorSequenceProvider };

                auto samples = recording.getSamples();
                size_t idx = 0;

                using OptionalTicketT = std::optional<typename Signal<gsl::span<const DescriptorT>>::TicketT>;
                OptionalTicketT maxScoreDescriptorHandlerTicket{ descriptorSignal.addHandler(
                    [this, &idx, &samples, &output](gsl::span<const DescriptorT> sequence) {
                        if (sequence.size() <= m_trimmedSequenceLength)
                        {
                            return;
                        }

                        gsl::span<const DescriptorT> trimmedSequence{
                            &sequence[sequence.size() - m_trimmedSequenceLength], m_trimmedSequenceLength};

                        for (size_t templateIdx = 0; templateIdx < m_templates.size(); ++templateIdx)
                        {
                            output << samples[idx].Timestamp << ",";
                            output << templateIdx << ",";

                            gsl::span<const DescriptorT> t{ m_templates[templateIdx] };

                            decltype(m_tuning) tuning{};
                            std::fill(tuning.begin(), tuning.end(), descriptor::NULL_TUNING);
                            for (size_t tuningIdx = 0; tuningIdx < tuning.size(); ++tuningIdx)
                            {
                                tuning[tuningIdx] = 1;

                                auto distanceFunction{ [&tuning](const DescriptorT& a, const DescriptorT& b) {
                                    return DescriptorT::Distance(a, b, tuning);
                                } };
                                auto [distance, imageSize] = DynamicTimeWarping::InjectiveDistanceAndImageSize(trimmedSequence, t, distanceFunction, m_minimumImageRatio);
                                output << distance << ",";

                                tuning[tuningIdx] = descriptor::NULL_TUNING;
                            }

                            output << std::endl;
                        }
                    }) };

                for (; idx < samples.size(); ++idx)
                {
                    auto& sample = samples[idx];
                    inputSamplesHandlers.apply_to_all([&sample](auto& callable) mutable { callable(sample); });
                }

                maxScoreDescriptorHandlerTicket.reset();
            }

        private:
            typename Signal<gsl::span<const DescriptorT>>::TicketT m_ticket;
            std::vector<std::vector<DescriptorT>> m_templates{};
            std::vector<std::vector<DescriptorT>> m_countertemplates{};
            size_t m_trimmedSequenceLength{};
            NumberT m_minimumImageRatio{ 0.8 }; // TODO: Parameterize?
            std::array<NumberT, DescriptorT::DEFAULT_TUNING.size()> m_tuning{};
            std::function<NumberT(NumberT)> m_scoringFunction{};
            bool m_recognition{ false };

            void initializeTemplates(gsl::span<const action::Example> examples, gsl::span<const action::Example> counterexamples)
            {
                auto initializeTemplatesFromExamples = [this](gsl::span<const action::Example> examples, std::vector<std::vector<DescriptorT>>& templates) {
                    for (const auto& example : examples)
                    {
                        templates.emplace_back();
                        auto& sequence = templates.back();

                        auto samples = example.getRecording().getSamples();
                        InputSample mostRecentSample{};

                        size_t idx = 0;
                        while (idx < samples.size() - 1 && samples[idx + 1].Timestamp < example.getStartTimestamp())
                        {
                            ++idx;
                        }
                        while (idx < samples.size() - 1 && samples[idx + 1].Timestamp < example.getEndTimestamp())
                        {
                            descriptor::extendSequence(samples[idx], sequence, mostRecentSample, DescriptorT::DEFAULT_TUNING);
                            ++idx;
                        }
                    }
                };

                initializeTemplatesFromExamples(examples, m_templates);
                initializeTemplatesFromExamples(counterexamples, m_countertemplates);

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

            void calculateTuning()
            {
                std::array<NumberT, DescriptorT::DEFAULT_TUNING.size()> averageDistances{};
                const NumberT templatesChooseTwo = (m_templates.size() * (m_templates.size() - 1)) / 2.;
                std::fill(m_tuning.begin(), m_tuning.end(), descriptor::NULL_TUNING);
                if (m_templates.size() > 1)
                {
                    for (size_t idx = 0; idx < m_tuning.size(); ++idx)
                    {
                        m_tuning[idx] = 1.;
                        for (size_t l = 0; l < m_templates.size(); ++l)
                        {
                            for (size_t r = l + 1; r < m_templates.size(); ++r)
                            {
                                auto distance = calculateNormalizedSequenceDistance(m_templates[l], m_templates[r]);
                                averageDistances[idx] += distance;
                            }
                        }
                        averageDistances[idx] /= templatesChooseTwo;
                        m_tuning[idx] = descriptor::NULL_TUNING;
                    }
                }

                constexpr NumberT defaultTuningWeight{ 1. };
                NumberT t = m_templates.size() > 1 ? 1. / (templatesChooseTwo + 1) : defaultTuningWeight;
                for (size_t idx = 0; idx < m_tuning.size(); ++idx)
                {
                    m_tuning[idx] = 1;/*t * DescriptorT::DEFAULT_TUNING[idx] +
                        (1. - t) * (1. / std::max(averageDistances[idx], kAverageDistanceEpsilon));*/
                }
            }

            void createAverageMinDistanceScoringFunction()
            {
                std::vector<NumberT> minDistances{};
                minDistances.resize(m_templates.size(), std::numeric_limits<NumberT>::max());
                for (size_t l = 0; l < m_templates.size(); ++l)
                {
                    for (size_t r = l + 1; r < m_templates.size(); ++r)
                    {
                        NumberT distance = calculateNormalizedSequenceDistance(m_templates[l], m_templates[r]);
                        minDistances[l] = std::min(distance, minDistances[l]);
                        minDistances[r] = std::min(distance, minDistances[r]);
                    }
                }

                NumberT averageMinDistance = 0.;
                for (size_t idx = 0; idx < m_templates.size(); ++idx)
                {
                    averageMinDistance += minDistances[idx];
                }
                averageMinDistance /= m_templates.size();

                m_scoringFunction = [this, averageMinDistance](NumberT distance)
                {
                    return std::max(1. - std::pow(distance / (m_sensitivity * 3. * averageMinDistance), 2.), 0.);
                };
            }

            void createAverageAverageDistanceScoringFunction()
            {
                std::vector<NumberT> averageDistances{};
                averageDistances.resize(m_templates.size());
                for (size_t l = 0; l < m_templates.size(); ++l)
                {
                    for (size_t r = l + 1; r < m_templates.size(); ++r)
                    {
                        NumberT distance = calculateNormalizedSequenceDistance(m_templates[l], m_templates[r]);
                        averageDistances[l] += distance;
                        averageDistances[r] += distance;
                    }
                    averageDistances[l] /= m_templates.size();
                }

                constexpr NumberT DEFAULT_AVERAGE_DISTANCE{ 0.1 };
                constexpr NumberT DEFAULT_AVERAGE_DISTANCE_WEIGHT{ 1. };
                NumberT averageAverageDistance{ DEFAULT_AVERAGE_DISTANCE * DEFAULT_AVERAGE_DISTANCE_WEIGHT };
                for (size_t idx = 0; idx < m_templates.size(); ++idx)
                {
                    averageAverageDistance += averageDistances[idx];
                }
                averageAverageDistance /= (m_templates.size() + DEFAULT_AVERAGE_DISTANCE_WEIGHT);

                m_scoringFunction = [this, averageAverageDistance](NumberT distance)
                {
                    return std::max(1. - std::pow(distance / (m_sensitivity * 3. * averageAverageDistance), 2.), 0.);
                };
            }

            void createUnitScoringFunction()
            {
                m_scoringFunction = [this](NumberT distance)
                {
                    return std::max(1. - std::pow(distance / (3.16228 * m_sensitivity), 2.), 0.);
                };
            }

            void createScoringFunction()
            {
                // createAverageMinDistanceScoringFunction();
                // createAverageAverageDistanceScoringFunction();
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
                NumberT maxScore = calculateScore(m_templates[0]);
                size_t maxIdx = 0;
                for (size_t idx = 1; idx < m_templates.size(); ++idx)
                {
                    NumberT score = calculateScore(m_templates[idx]);
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
                        recording.addSample(std::move(sampleCopy));
                    }
                }
                m_canonicalRecording = std::make_unique<Recording>(std::move(recording));
            }

            std::tuple<NumberT, size_t> calculateSequenceDistance(
                gsl::span<const DescriptorT> a,
                gsl::span<const DescriptorT> b) const
            {
                bool aLongerThanB = a.size() > b.size();
                auto& longer = aLongerThanB ? a : b;
                auto& shorter = aLongerThanB ? b : a;

                auto distanceFunction{ [this](const DescriptorT& a, const DescriptorT& b) {
                  return DescriptorT::Distance(a, b, m_tuning);
                } };
                return DynamicTimeWarping::InjectiveDistanceAndImageSize(longer, shorter, distanceFunction, m_minimumImageRatio);
            }

            NumberT calculateNormalizedSequenceDistance(
                gsl::span<const DescriptorT> a,
                gsl::span<const DescriptorT> b) const
            {
                auto [distance, imageSize] = calculateSequenceDistance(a, b);
                return distance;
            }

            NumberT calculateScore(gsl::span<const DescriptorT> sequence) const
            {
                // Early-out if the provided sequence is too short.
                if (sequence.size() <= m_trimmedSequenceLength)
                {
                    return 0.;
                }

                gsl::span<const DescriptorT> trimmedSequence{ &sequence[sequence.size() - m_trimmedSequenceLength], m_trimmedSequenceLength };

                // Calculate the base score based on proximity to templates
                NumberT score{ 0 };
                NumberT minDistance = std::numeric_limits<NumberT>::max();
                for (const auto& t : m_templates)
                {
                    NumberT distance = calculateNormalizedSequenceDistance(trimmedSequence, t);
                    score = std::max(m_scoringFunction(distance), score);
                    minDistance = std::min(distance, minDistance);
                }

                // Clamp score to [0, 1] range.
                score = std::clamp<NumberT>(score, 0, 1);

                // Penalize score based on relative proximity to countertemplates
                NumberT minCounterDistance = std::numeric_limits<NumberT>::max();
                for (const auto& t : m_countertemplates)
                {
                    NumberT distance = calculateNormalizedSequenceDistance(trimmedSequence, t);
                    minCounterDistance = std::min(distance, minCounterDistance);
                }
                if (minDistance >= minCounterDistance)
                {
                    score *= minCounterDistance / minDistance;
                }

                return score;
            }

            void handleSequence(gsl::span<const DescriptorT> sequence)
            {
                auto score = calculateScore(sequence);
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

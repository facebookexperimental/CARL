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
        constexpr double kAverageDistanceEpsilon{ 0.000000001 };

        // Slightly less naive holistic resampling. Still linear, though.
        action::Recording resample(const action::Example& original, double frameDuration)
        {
            auto samples = original.getRecording().getSamples();
            size_t samplesIdx = 0;

            std::vector<InputSample> newSamples{};
            size_t newSamplesCount = std::max<size_t>(
                static_cast<size_t>(
                    (original.getEndTimestamp() - original.getStartTimestamp()) / frameDuration),
                1);
            newSamples.reserve(newSamplesCount);

            for (size_t idx = 0; idx < newSamplesCount; ++idx)
            {
                double timestamp = original.getStartTimestamp() + frameDuration * idx;

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
                    double t = (timestamp - original.getStartTimestamp()) /
                        (original.getEndTimestamp() - original.getStartTimestamp());
                    newSamples.emplace_back(InputSample::lerp(samples[samplesIdx - 1], samples[samplesIdx], t));
                }
            }

            return { action::InProgressRecording{std::move(newSamples)} };
        }

        std::vector<action::Example> expandExamples(
            gsl::span<const action::Example> examples,
            double durationT = 0.2)
        {
            std::vector<action::Example> expandedExamples{};

            for (const auto& example : examples)
            {
                for (auto endT = example.getStartTimestamp() + durationT; endT < example.getEndTimestamp(); endT += durationT)
                {
                    expandedExamples.emplace_back(example.getRecording(), endT - durationT, endT);
                }
            }

            return expandedExamples;
        }
    }

    class action::Recognizer::Impl
    {
        friend class action::Recognizer;

    public:
        Impl(double sensitivity)
            : m_sensitivity{ sensitivity }
        {
        }

        virtual ~Impl() = default;

    protected:
        virtual Example createAutoTrimmedExample(const Recording&) const = 0;

    protected:
        arcana::weak_table<Signal<bool>::HandlerT> m_whenRecognitionChangedHandlers{};
        std::atomic<double> m_currentScore{};
        std::unique_ptr<Recording> m_canonicalRecording{};
        double m_sensitivity{};
    };

    namespace
    {
        template <typename DescriptorT>
        class RecognizerImpl : public action::Recognizer::Impl
        {
        public:
            RecognizerImpl(Session& session, gsl::span<const action::Example> examples, double sensitivity)
                : Recognizer::Impl{ sensitivity },
                m_ticket{ Session::Impl::getFromSession(session).addHandler<DescriptorT>(
                    [this](gsl::span<const DescriptorT> sequence) { handleSequence(sequence); }) }
            {
                initializeTemplates(session, examples);
                calculateTuning();
                createScoringFunction();
                createCanonicalRecording(examples);
            }

        protected:
            Example createAutoTrimmedExample(const Recording& recording) const override
            {
                using SignalT = Signal<gsl::span<const InputSample>>;
                arcana::weak_table<SignalT::HandlerT> inputSamplesHandlers{};
                SignalT inputSampleSignal{ inputSamplesHandlers };
                typename DescriptorSequence<DescriptorT, 1024>::Provider descriptorSequenceProvider{ inputSampleSignal }; // TODO: Magic number
                Signal<gsl::span<const DescriptorT>>& descriptorSignal{ descriptorSequenceProvider };

                size_t idx = 0;
                auto samples = recording.getSamples();

                auto startT = samples[0].Timestamp;
                auto endT = samples[samples.size() - 1].Timestamp;
                auto distance = std::numeric_limits<double>::max();

                auto descriptorHandlerTicket{ descriptorSignal.addHandler(
                    [this, &idx, &samples, &startT, &endT, &distance](gsl::span<const DescriptorT> sequence) {
                        if (sequence.size() <= m_trimmedSequenceLength)
                        {
                          return;
                        }

                        gsl::span<const DescriptorT> trimmedSequence{
                            &sequence[sequence.size() - m_trimmedSequenceLength], m_trimmedSequenceLength};

                        for (const auto& t : m_templates)
                        {
                            auto [imageDistance, imageSize] = calculateSequenceDistance(trimmedSequence, t);
                            auto normalizedDistance = imageDistance / imageSize;
                            
                            if (normalizedDistance < distance)
                            {
                                startT = samples[idx - imageSize].Timestamp;
                                endT = samples[idx].Timestamp;
                                distance = normalizedDistance;
                            }
                        }
                    }) };

                for (; idx < samples.size() - 1; ++idx)
                {
                    auto span = gsl::make_span(samples.data() + idx, 2);
                    inputSamplesHandlers.apply_to_all([span](auto& callable) mutable { callable(span); });
                }

                return{ recording, startT, endT };
            }

        private:
            typename Signal<gsl::span<const DescriptorT>>::TicketT m_ticket;
            std::vector<std::vector<DescriptorT>> m_templates{};
            size_t m_trimmedSequenceLength{};
            double m_minimumImageRatio{ 0.8 }; // TODO: Parameterize?
            std::array<double, DescriptorT::DEFAULT_TUNING.size()> m_tuning{};
            std::function<double(double)> m_scoringFunction{};
            bool m_recognition{ false };

            void initializeTemplates(Session& session, gsl::span<const action::Example> examples)
            {
                auto& sessionImpl = Session::Impl::getFromSession(session);

                std::vector<action::Recording> resampledRecordings{};
                resampledRecordings.reserve(examples.size());

                for (const auto& example : examples)
                {
                    auto recording = resample(example, sessionImpl.frameDuration);
                    auto samples = recording.getSamples();

                    m_templates.emplace_back();
                    auto& sequence = m_templates.back();
                    sequence.reserve(samples.size() - 1);
                    if (samples.size() == 1)
                    {
                        // For recordings short enough that they only contain one sample, consider them to
                        // represent a static pose.
                        auto descriptor = DescriptorT::TryCreate(samples[0], samples[0]);
                        if (descriptor.has_value())
                        {
                            sequence.emplace_back(std::move(descriptor.value()));
                        }
                    }
                    else
                    {
                        for (size_t idx = 0; idx < samples.size() - 1; ++idx)
                        {
                            auto descriptor = DescriptorT::TryCreate(samples[idx + 1], samples[idx]);
                            if (descriptor.has_value())
                            {
                                sequence.emplace_back(std::move(descriptor.value()));
                            }
                        }
                    }

                    // TODO: Parameterize this calculation, instead of hard-coding 5/4ths?
                    m_trimmedSequenceLength = std::max(m_trimmedSequenceLength, (5 * sequence.size()) / 4);
                }
            }

            void calculateTuning()
            {
                std::array<double, DescriptorT::DEFAULT_TUNING.size()> averageDistances{};
                const double templatesChooseTwo = (m_templates.size() * (m_templates.size() - 1)) / 2.;
                if (m_templates.size() > 1)
                {
                    for (size_t idx = 0; idx < m_tuning.size(); ++idx)
                    {
                        m_tuning[idx] = 1.;
                        for (size_t l = 0; l < m_templates.size(); ++l)
                        {
                            for (size_t r = l + 1; r < m_templates.size(); ++r)
                            {
                                averageDistances[idx] += calculateNormalizedSequenceDistance(m_templates[l], m_templates[r]);
                            }
                        }
                        averageDistances[idx] /= templatesChooseTwo;
                        m_tuning[idx] = 0.;
                    }
                }

                constexpr double defaultTuningWeight{ 1. };
                double t = m_templates.size() > 1 ? 1. / (templatesChooseTwo + 1) : defaultTuningWeight;
                for (size_t idx = 0; idx < m_tuning.size(); ++idx)
                {
                    m_tuning[idx] = t * DescriptorT::DEFAULT_TUNING[idx] +
                        (1. - t) * (1. / std::max(averageDistances[idx], kAverageDistanceEpsilon));
                }
            }

            void createAverageMinDistanceScoringFunction()
            {
                std::vector<double> minDistances{};
                minDistances.resize(m_templates.size(), std::numeric_limits<double>::max());
                for (size_t l = 0; l < m_templates.size(); ++l)
                {
                    for (size_t r = l + 1; r < m_templates.size(); ++r)
                    {
                        double distance = calculateNormalizedSequenceDistance(m_templates[l], m_templates[r]);
                        minDistances[l] = std::min(distance, minDistances[l]);
                        minDistances[r] = std::min(distance, minDistances[r]);
                    }
                }

                double averageMinDistance = 0.;
                for (size_t idx = 0; idx < m_templates.size(); ++idx)
                {
                    averageMinDistance += minDistances[idx];
                }
                averageMinDistance /= m_templates.size();

                m_scoringFunction = [this, averageMinDistance](double distance)
                {
                    return std::max(1. - std::pow(distance / (m_sensitivity * 3. * averageMinDistance), 2.), 0.);
                };
            }

            void createAverageAverageDistanceScoringFunction()
            {
                std::vector<double> averageDistances{};
                averageDistances.resize(m_templates.size());
                for (size_t l = 0; l < m_templates.size(); ++l)
                {
                    for (size_t r = l + 1; r < m_templates.size(); ++r)
                    {
                        double distance = calculateNormalizedSequenceDistance(m_templates[l], m_templates[r]);
                        averageDistances[l] += distance;
                        averageDistances[r] += distance;
                    }
                    averageDistances[l] /= m_templates.size();
                }

                constexpr double DEFAULT_AVERAGE_DISTANCE{ 10. };
                constexpr double DEFAULT_AVERAGE_DISTANCE_WEIGHT{ 1. };
                double averageAverageDistance{ DEFAULT_AVERAGE_DISTANCE * DEFAULT_AVERAGE_DISTANCE_WEIGHT };
                for (size_t idx = 0; idx < m_templates.size(); ++idx)
                {
                    averageAverageDistance += averageDistances[idx];
                }
                averageAverageDistance /= (m_templates.size() + DEFAULT_AVERAGE_DISTANCE_WEIGHT);

                m_scoringFunction = [this, averageAverageDistance](double distance)
                {
                    return std::max(1. - std::pow(distance / (m_sensitivity * 3. * averageAverageDistance), 2.), 0.);
                };
            }

            void createScoringFunction()
            {
                // createAverageMinDistanceScoringFunction();
                createAverageAverageDistanceScoringFunction();
            }

            void createCanonicalRecording(gsl::span<const Example> examples)
            {
                double maxScore = calculateScore(m_templates[0]);
                size_t maxIdx = 0;
                for (size_t idx = 1; idx < m_templates.size(); ++idx)
                {
                    double score = calculateScore(m_templates[idx]);
                    if (score > maxScore)
                    {
                        maxScore = score;
                        maxIdx = idx;
                    }
                }

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

            std::tuple<double, size_t> calculateSequenceDistance(
                gsl::span<const DescriptorT> a,
                gsl::span<const DescriptorT> b) const
            {
                bool aLongerThanB = a.size() > b.size();
                auto& longer = aLongerThanB ? a : b;
                auto& shorter = aLongerThanB ? b : a;

                auto distanceFunction{ [this](const DescriptorT& a, const DescriptorT& b) {
                  return DescriptorT::Distance(a, b, m_tuning);
                } };
                return DynamicTimeWarping::InjectiveDistance(longer, shorter, distanceFunction, m_minimumImageRatio);
            }

            double calculateNormalizedSequenceDistance(
                gsl::span<const DescriptorT> a,
                gsl::span<const DescriptorT> b) const
            {
                auto [distance, imageSize] = calculateSequenceDistance(a, b);
                return distance / imageSize;
            }

            double calculateScore(gsl::span<const DescriptorT> sequence)
            {
                if (sequence.size() <= m_trimmedSequenceLength)
                {
                    return 0.;
                }
                gsl::span<const DescriptorT> trimmedSequence{ &sequence[sequence.size() - m_trimmedSequenceLength], m_trimmedSequenceLength };

                double score = 0.;
                for (const auto& t : m_templates)
                {
                    double distance = calculateNormalizedSequenceDistance(trimmedSequence, t);
                    score += m_scoringFunction(distance);
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
            switch (definition.getDescriptorType())
            {
            case action::Definition::ActionType::LeftHandPose:
                examples = expandExamples(definition.getExamples());
                return std::make_unique<RecognizerImpl<descriptor::HandPose<descriptor::Handedness::LeftHanded>>>(session, examples, definition.DefaultSensitivity);
            case action::Definition::ActionType::LeftHandGesture:
                return std::make_unique<RecognizerImpl<descriptor::HandGesture<descriptor::Handedness::LeftHanded>>>(session, definition.getExamples(), definition.DefaultSensitivity);
            case action::Definition::ActionType::RightHandPose:
                examples = expandExamples(definition.getExamples());
                return std::make_unique<RecognizerImpl<descriptor::HandPose<descriptor::Handedness::RightHanded>>>(session, examples, definition.DefaultSensitivity);
            case action::Definition::ActionType::RightHandGesture:
                return std::make_unique<RecognizerImpl<descriptor::HandGesture<descriptor::Handedness::RightHanded>>>(session, definition.getExamples(), definition.DefaultSensitivity);
            case action::Definition::ActionType::TwoHandGesture:
                return std::make_unique<RecognizerImpl<descriptor::TwoHandGesture>>(session, definition.getExamples(), definition.DefaultSensitivity);
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
        return m_impl->m_currentScore;
    }

    RecordingInspector Recognizer::getCanonicalRecordingInspector() const
    {
        return m_impl->m_canonicalRecording->getInspector();
    }

    void Recognizer::setSensitivity(double sensitivity)
    {
        m_impl->m_sensitivity = sensitivity;
    }

    Example Recognizer::createAutoTrimmedExample(const Recording& recording) const
    {
        return m_impl->createAutoTrimmedExample(recording);
    }
}

/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <carl/DefinitionBuilder.h>

#include "AutoTrimImpl.h"

#include <carl/Session.h>
#include <carl/ActionTypeDefinitions.h>
#include <carl/descriptor/SequenceOperations.h>
#include <carl/descriptor/TimestampedDescriptor.h>
#include <carl/DynamicTimeWarping.h>

#include <algorithm>
#include <numeric>

namespace carl::action
{
    namespace
    {
        constexpr double STATIC_GESTURE_SEED_DURATION = 0.5;
        constexpr double MOTION_ENERGY_PEAK_RATIO = 1.5;

        template<typename DescriptorT>
        struct TimestampedSequenceInfo
        {
            std::vector<descriptor::TimestampedDescriptor<DescriptorT>> sequence{};
            double startTimestamp{};
            double endTimestamp{};
        };

        template<typename DescriptorT>
        TimestampedSequenceInfo<DescriptorT> buildTimestampedSequence(const Recording& recording)
        {
            constexpr auto tryCreate = [](const InputSample& current, const InputSample& prior) {
                return descriptor::TimestampedDescriptor<DescriptorT>::TryCreate(current, prior);
            };

            std::vector<descriptor::TimestampedDescriptor<DescriptorT>> sequence{};
            auto samples = recording.getSamples();
            InputSample mostRecentSample = samples.front();
            for (const auto& sample : samples)
            {
                descriptor::extendSequence<descriptor::TimestampedDescriptor<DescriptorT>>(
                    sample, sequence, mostRecentSample, DescriptorT::DEFAULT_TUNING, tryCreate);
            }

            double start = samples.front().Timestamp;
            double end = samples.back().Timestamp;

            return { std::move(sequence), start, end };
        }

        /// Finds the time window of the given duration that contains the
        /// most descriptors (highest motion energy). Returns the center
        /// timestamp of the densest window and its descriptor count.
        template<typename DescriptorT>
        std::pair<double, size_t> findDensestWindow(
            const std::vector<descriptor::TimestampedDescriptor<DescriptorT>>& sequence,
            double windowDuration)
        {
            if (sequence.empty())
            {
                return { 0., 0 };
            }

            double bestCenter = sequence.front().getTimestamp() + windowDuration / 2.;
            size_t bestCount = 0;

            size_t windowEnd = 0;
            for (size_t windowStart = 0; windowStart < sequence.size(); ++windowStart)
            {
                double startTime = sequence[windowStart].getTimestamp();
                double endTime = startTime + windowDuration;

                while (windowEnd < sequence.size() && sequence[windowEnd].getTimestamp() <= endTime)
                {
                    ++windowEnd;
                }

                size_t count = windowEnd - windowStart;
                if (count > bestCount)
                {
                    bestCount = count;
                    bestCenter = startTime + windowDuration / 2.;
                }
            }

            return { bestCenter, bestCount };
        }

        /// Determines the window duration to use for motion-energy seeding.
        double computeWindowDuration(gsl::span<const Recording> recordings, double expectedActionDuration)
        {
            if (expectedActionDuration > 0.)
            {
                return expectedActionDuration;
            }

            double shortest = std::numeric_limits<double>::max();
            for (const auto& recording : recordings)
            {
                auto samples = recording.getSamples();
                double duration = samples.back().Timestamp - samples.front().Timestamp;
                shortest = std::min(shortest, duration);
            }

            return shortest / 2.;
        }

        /// Picks a half-second chunk from the middle of the recording
        /// when no significant motion is detected (static gesture case).
        std::pair<double, double> staticGestureSeed(const Recording& recording)
        {
            auto samples = recording.getSamples();
            double start = samples.front().Timestamp;
            double end = samples.back().Timestamp;
            double mid = (start + end) / 2.;
            double halfDur = STATIC_GESTURE_SEED_DURATION / 2.;
            return { std::max(mid - halfDur, start), std::min(mid + halfDur, end) };
        }

        template<typename DescriptorT>
        std::vector<Example> createExamplesImpl(
            gsl::span<const Recording> recordings,
            double expectedActionDuration)
        {
            if (recordings.empty())
            {
                return {};
            }

            // Step 1: Motion-energy seed
            double windowDuration = computeWindowDuration(recordings, expectedActionDuration);

            struct SeedCandidate
            {
                size_t recordingIdx;
                double center;
                size_t peakCount;
                double avgCount;
            };

            std::vector<TimestampedSequenceInfo<DescriptorT>> tsSequences{};
            tsSequences.reserve(recordings.size());
            SeedCandidate bestSeed{ 0, 0., 0, 0. };

            for (size_t i = 0; i < recordings.size(); ++i)
            {
                auto info = buildTimestampedSequence<DescriptorT>(recordings[i]);
                auto [center, peakCount] = findDensestWindow<DescriptorT>(info.sequence, windowDuration);

                double totalDuration = info.endTimestamp - info.startTimestamp;
                double avgCount = totalDuration > 0. ? static_cast<double>(info.sequence.size()) / totalDuration * windowDuration : 0.;

                double peakRatio = avgCount > 0. ? static_cast<double>(peakCount) / avgCount : 0.;

                if (i == 0 || peakRatio > (bestSeed.avgCount > 0. ? static_cast<double>(bestSeed.peakCount) / bestSeed.avgCount : 0.))
                {
                    bestSeed = { i, center, peakCount, avgCount };
                }

                tsSequences.push_back(std::move(info));
            }

            // Determine seed boundaries
            double seedStart, seedEnd;
            bool isStaticGesture = bestSeed.peakCount == 0 ||
                (bestSeed.avgCount > 0. && static_cast<double>(bestSeed.peakCount) / bestSeed.avgCount < MOTION_ENERGY_PEAK_RATIO);

            if (isStaticGesture)
            {
                auto [s, e] = staticGestureSeed(recordings[bestSeed.recordingIdx]);
                seedStart = s;
                seedEnd = e;
            }
            else
            {
                seedStart = bestSeed.center - windowDuration / 2.;
                seedEnd = bestSeed.center + windowDuration / 2.;
                auto samples = recordings[bestSeed.recordingIdx].getSamples();
                seedStart = std::max(seedStart, static_cast<double>(samples.front().Timestamp));
                seedEnd = std::min(seedEnd, static_cast<double>(samples.back().Timestamp));
            }

            Example seedExample{ recordings[bestSeed.recordingIdx], seedStart, seedEnd };

            // Single recording: just return the seed
            if (recordings.size() == 1)
            {
                std::vector<Example> result{};
                result.push_back(std::move(seedExample));
                return result;
            }

            // Step 2: First bootstrap pass
            constexpr auto tryCreate = [](const InputSample& current, const InputSample& prior) {
                return DescriptorT::TryCreate(current, prior);
            };
            constexpr auto tryCreateTimestamped = [](const InputSample& current, const InputSample& prior) {
                return descriptor::TimestampedDescriptor<DescriptorT>::TryCreate(current, prior);
            };

            std::vector<Example> seedExamples{};
            seedExamples.push_back(seedExample);

            auto seedTemplates = descriptor::createDescriptorSequencesFromExamples<DescriptorT>(
                gsl::span<const Example>{ seedExamples }, DescriptorT::DEFAULT_TUNING, tryCreate);

            std::vector<Example> firstPassExamples{};
            firstPassExamples.reserve(recordings.size());

            for (size_t i = 0; i < recordings.size(); ++i)
            {
                if (i == bestSeed.recordingIdx)
                {
                    firstPassExamples.push_back(std::move(seedExample));
                    continue;
                }

                gsl::span<const std::vector<DescriptorT>> templateSpan{ seedTemplates };
                auto [startT, endT] = autoTrimRecording<DescriptorT>(
                    recordings[i], templateSpan, DescriptorT::DEFAULT_TUNING, tryCreateTimestamped);
                firstPassExamples.emplace_back(recordings[i], startT, endT);
            }

            // Step 3: Refinement pass with real tuning
            auto tuning = DescriptorT::CalculateTuning(gsl::span<const Example>{ firstPassExamples });

            auto refinedTemplates = descriptor::createDescriptorSequencesFromExamples<DescriptorT>(
                gsl::span<const Example>{ firstPassExamples }, tuning, tryCreate);

            std::vector<Example> result{};
            result.reserve(recordings.size());

            for (size_t i = 0; i < recordings.size(); ++i)
            {
                gsl::span<const std::vector<DescriptorT>> templateSpan{ refinedTemplates };
                auto [startT, endT] = autoTrimRecording<DescriptorT>(
                    recordings[i], templateSpan, tuning, tryCreateTimestamped);
                result.emplace_back(recordings[i], startT, endT);
            }

            return result;
        }
    }

    std::vector<Example> DefinitionBuilder::createExamplesFromRecordings(
        ActionType actionType,
        gsl::span<const Recording> recordings,
        double expectedActionDuration)
    {
        return ActionTypes::invokeByActionType(actionType, [&](auto* tag) -> std::vector<Example> {
            using Association = std::remove_pointer_t<decltype(tag)>;
            using DescriptorT = typename Association::Descriptor;
            return createExamplesImpl<DescriptorT>(recordings, expectedActionDuration);
        });
    }
}

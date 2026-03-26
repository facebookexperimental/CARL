/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include <carl/descriptor/DescriptorUtils.h>
#include <carl/descriptor/SequenceOperations.h>
#include <carl/descriptor/Trivial.h>

namespace carl::descriptor
{
    // TODO: Review theory of delta descriptors and assess how they work with delta sampling. Make sure creation isn't oscillating.
    template<Handedness Handedness>
    class WristRotation
    {
        static inline constexpr const char* ANALYSIS_DIMENSION_NAME = "Wrist Rotation";

    public:
        static constexpr bool ANCHOR_INDEPENDENT = true;
        static constexpr std::array<NumberT, 1> DEFAULT_TUNING{ 1. };

        static NumberT AnchorFreeDistance(
            const WristRotation& a,
            const WristRotation& b,
            gsl::span<const NumberT> tuning)
        {
            return Distance(a, a, b, b, tuning);
        }

        static NumberT AnchorDependentDistance(
            const WristRotation&,
            const WristRotation&,
            const WristRotation&,
            const WristRotation&,
            gsl::span<const NumberT>)
        {
            return 0;
        }

        static std::optional<WristRotation> TryCreate(
            const InputSample& sample,
            const InputSample& priorSample)
        {
            if (!sample.HmdPose.has_value())
            {
                return{};
            }

            auto optionalWristPose = getWristPose<Handedness>(sample);
            auto optionalPriorWristPose = getWristPose<Handedness>(priorSample);
            if (optionalWristPose.has_value() && optionalPriorWristPose.has_value())
            {
                return WristRotation{ optionalWristPose.value(), optionalPriorWristPose.value() };
            }
            return{};
        }

        static NumberT Distance(
            const WristRotation& a,
            const WristRotation&,
            const WristRotation& b,
            const WristRotation&,
            gsl::span<const NumberT> tuning)
        {
            return InternalNormalizedDistance(a, b, tuning);
        }

        static WristRotation Lerp(const WristRotation& a, const WristRotation& b, NumberT t)
        {
            QuaternionT aQuat{ a.m_deltaOrientation };
            QuaternionT bQuat{ b.m_deltaOrientation };
            WristRotation result{};
            result.m_deltaOrientation = aQuat.slerp(t, bQuat);
            return result;
        }

        template<typename ExamplesT>
        static std::array<NumberT, DEFAULT_TUNING.size()> CalculateTuning(const ExamplesT& examples)
        {
            if (examples.size() < 2)
            {
                return DEFAULT_TUNING;
            }

            constexpr auto tryCreate = [](const InputSample& current, const InputSample& prior) { return TryCreate(current, prior); };
            auto sequences = createDescriptorSequencesFromExamples<WristRotation<Handedness>>(examples, DEFAULT_TUNING, tryCreate);
            auto extendedSequences = createDescriptorSequencesFromExamples<WristRotation<Handedness>>(examples, DEFAULT_TUNING, tryCreate, 1.);

            auto tuning = DEFAULT_TUNING;
            NumberT maxAverageConnectionCost = 0;
            for (const auto& extendedSequence : extendedSequences)
            {
                for (const auto& sequence : sequences)
                {
                    auto distanceFunction = [](const auto& a, const auto&, const auto& b, const auto&) {
                        return InternalRawDistance(a, b);
                        };
                    auto result = DynamicTimeWarping::Match<const WristRotation<Handedness>>(extendedSequence, sequence, distanceFunction);
                    maxAverageConnectionCost = std::max<NumberT>(result.MaxConnectionCost / result.Connections, maxAverageConnectionCost);
                }
            }
            tuning[0] = std::max(maxAverageConnectionCost / IDENTICALITY_THRESHOLD, DEFAULT_TUNING[0]);
            return tuning;
        }

        template<bool NormalizeDistance>
        static auto Analyze(
            gsl::span<const WristRotation<Handedness>> target,
            gsl::span<const WristRotation<Handedness>> query,
            gsl::span<const NumberT> tuning)
        {
            std::array<AnalysisT, DEFAULT_TUNING.size()> results{};
            results[0] = { ANALYSIS_DIMENSION_NAME, IDENTICALITY_THRESHOLD, tuning[0], {} };
            auto& rows = std::get<3>(results[0]);
            auto distanceFunction = [tuning](const auto& a, const auto&, const auto& b, const auto&) {
                if constexpr (NormalizeDistance)
                {
                    return InternalNormalizedDistance(a, b, tuning);
                }
                else
                {
                    return InternalRawDistance(a, b);
                }
                };
            auto rowsCallback = [&rows](std::vector<DynamicTimeWarping::MatchResult<NumberT>> row) { rows.push_back(std::move(row)); };
            DynamicTimeWarping::Match<const WristRotation<Handedness>, decltype(distanceFunction), NumberT, true, decltype(rowsCallback)>(target, query, distanceFunction, 0, rowsCallback);
            return results;
        }

        WristRotation() = default;

    private:
        static inline constexpr NumberT IDENTICALITY_THRESHOLD{ 0.06 };
        static inline constexpr auto normalizeDistance{ createDistanceNormalizationFunction(IDENTICALITY_THRESHOLD) };
        trivial::Quaternion m_deltaOrientation{};

        WristRotation(const TransformT& wristPose, const TransformT& priorWristPose)
        {
            AngleAxisT angleAxis{ priorWristPose.rotation().inverse() * wristPose.rotation() };
            m_deltaOrientation = QuaternionT{ angleAxis };
        }

        static NumberT InternalRawDistance(const WristRotation& a, const WristRotation& b)
        {
            return a.m_deltaOrientation.angularDistance(b.m_deltaOrientation);
        }

        static NumberT InternalNormalizedDistance(const WristRotation& a, const WristRotation& b, gsl::span<const NumberT> tuning)
        {
            return normalizeDistance(InternalRawDistance(a, b), tuning.front());
        }
    };
}

/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include "DescriptorUtils.h"
#include "EgocentricTemporalSpace.h"
#include "SequenceOperations.h"
#include "Trivial.h"

namespace carl::descriptor
{
    template<Handedness Handedness>
    class EgocentricWristDisplacement
    {
        static inline constexpr const char* ANALYSIS_DIMENSION_NAME = "Egocentric Wrist Displacement";

    public:
        static constexpr std::array<NumberT, 1> DEFAULT_TUNING{ 1. };

        static std::optional<EgocentricWristDisplacement> TryCreate(
            const InputSample& sample,
            const InputSample& priorSample)
        {
            if (!sample.HmdPose.has_value())
            {
                return{};
            }

            if constexpr (Handedness == Handedness::LeftHanded)
            {
                if (sample.LeftWristPose.has_value())
                {
                    return EgocentricWristDisplacement{ sample, priorSample };
                }
            }
            else if constexpr (Handedness == Handedness::RightHanded)
            {
                if (sample.RightWristPose.has_value())
                {
                    return EgocentricWristDisplacement{ sample, priorSample };
                }
            }
            return{};
        }

        static NumberT Distance(
            const EgocentricWristDisplacement& a,
            const EgocentricWristDisplacement& a0,
            const EgocentricWristDisplacement& b,
            const EgocentricWristDisplacement& b0,
            gsl::span<const NumberT> tuning)
        {
            return InternalNormalizedDistance(a, a0, b, b0, tuning);
        }

        static EgocentricWristDisplacement Lerp(
            const EgocentricWristDisplacement& a,
            const EgocentricWristDisplacement& b,
            NumberT t)
        {
            EgocentricWristDisplacement result{};
            result.m_position = (1 - t) * VectorT{ a.m_position } + t * VectorT{ b.m_position };
            result.m_inverseEts = math::Lerp(a.m_inverseEts, b.m_inverseEts, t);
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
            auto sequences = createDescriptorSequencesFromExamples<EgocentricWristDisplacement<Handedness>>(examples, DEFAULT_TUNING, tryCreate);
            auto extendedSequences = createDescriptorSequencesFromExamples<EgocentricWristDisplacement<Handedness>>(examples, DEFAULT_TUNING, tryCreate, 1.);

            auto tuning = DEFAULT_TUNING;
            NumberT maxAverageConnectionCost = 0;
            for (const auto& extendedSequence : extendedSequences)
            {
                for (const auto& sequence : sequences)
                {
                    auto distanceFunction = [](const auto& a, const auto& a0, const auto& b, const auto& b0) {
                        return InternalRawDistance(a, a0, b, b0);
                        };
                    auto result = DynamicTimeWarping::Match<const EgocentricWristDisplacement<Handedness>>(extendedSequence, sequence, distanceFunction);
                    maxAverageConnectionCost = std::max<NumberT>(result.MaxConnectionCost / result.Connections, maxAverageConnectionCost);
                }
            }
            tuning[0] = std::max(maxAverageConnectionCost / IDENTICALITY_THRESHOLD, DEFAULT_TUNING[0]);
            return tuning;
        }

        template<bool NormalizeDistance>
        static auto Analyze(
            gsl::span<const EgocentricWristDisplacement<Handedness>> target,
            gsl::span<const EgocentricWristDisplacement<Handedness>> query,
            gsl::span<const NumberT> tuning)
        {
            std::array<AnalysisT, DEFAULT_TUNING.size()> results{};
            results[0] = { ANALYSIS_DIMENSION_NAME, IDENTICALITY_THRESHOLD, tuning.front(), {} };
            auto& rows = std::get<3>(results[0]);
            auto distanceFunction = [tuning](const auto& a, const auto& a0, const auto& b, const auto& b0) {
                if constexpr (NormalizeDistance)
                {
                    return InternalNormalizedDistance(a, a0, b, b0, tuning);
                }
                else
                {
                    return InternalRawDistance(a, a0, b, b0);
                }
                };
            auto rowsCallback = [&rows](std::vector<DynamicTimeWarping::MatchResult<NumberT>> row) { rows.push_back(std::move(row)); };
            DynamicTimeWarping::Match<const EgocentricWristDisplacement<Handedness>, decltype(distanceFunction), NumberT, true, decltype(rowsCallback)>(target, query, distanceFunction, 0, rowsCallback);
            return results;
        }

        EgocentricWristDisplacement() = default;

    private:
        static inline constexpr NumberT IDENTICALITY_THRESHOLD{ 0.02 };
        static inline constexpr auto normalizeDistance{ createDistanceNormalizationFunction(IDENTICALITY_THRESHOLD) };
        trivial::Point m_position{};
        trivial::Transform m_inverseEts{};

        EgocentricWristDisplacement(const InputSample& sample, const InputSample&)
        {
            constexpr bool isLeftHanded = Handedness == Handedness::LeftHanded;
            auto& wristPose = isLeftHanded ? sample.LeftWristPose.value() : sample.RightWristPose.value();

            auto ets = EgocentricTemporalSpace::getPose(wristPose.translation(), sample.HmdPose.value());
            m_position = ets.translation();
            m_inverseEts = ets.inverse();
        }

        static NumberT InternalRawDistance(
            const EgocentricWristDisplacement& a,
            const EgocentricWristDisplacement& a0,
            const EgocentricWristDisplacement& b,
            const EgocentricWristDisplacement& b0)
        {
            auto aPos = a0.m_inverseEts * a.m_position;
            auto bPos = b0.m_inverseEts * b.m_position;
            return aPos.distance(bPos);
        }

        static NumberT InternalNormalizedDistance(
            const EgocentricWristDisplacement& a,
            const EgocentricWristDisplacement& a0,
            const EgocentricWristDisplacement& b,
            const EgocentricWristDisplacement& b0,
            gsl::span<const NumberT> tuning)
        {
            return normalizeDistance(InternalRawDistance(a, a0, b, b0), tuning.front());
        }
    };
}

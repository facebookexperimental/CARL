/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include <carl/descriptor/DescriptorUtils.h>
#include <carl/descriptor/EgocentricTemporalSpace.h>
#include <carl/descriptor/SequenceOperations.h>
#include <carl/descriptor/Trivial.h>

namespace carl::descriptor
{
    class EgocentricRelativeWristPosition
    {
        static inline constexpr const char* ANALYSIS_DIMENSION_NAME = "Egocentric Relative Wrist Position";

    public:
        static constexpr std::array<NumberT, 1> DEFAULT_TUNING{ 1. };

        static NumberT AbsoluteDistance(
            const EgocentricRelativeWristPosition& a,
            const EgocentricRelativeWristPosition& b,
            gsl::span<const NumberT> tuning)
        {
            return Distance(a, a, b, b, tuning);
        }

        static std::optional<EgocentricRelativeWristPosition> TryCreate(
            const InputSample& sample,
            const InputSample& priorSample)
        {
            if (sample.HmdPose.has_value() &&
                sample.LeftWristPose.has_value() &&
                sample.RightWristPose.has_value())
            {
                return EgocentricRelativeWristPosition{ sample, priorSample };
            }
            return{};
        }

        static NumberT Distance(
            const EgocentricRelativeWristPosition& a,
            const EgocentricRelativeWristPosition&,
            const EgocentricRelativeWristPosition& b,
            const EgocentricRelativeWristPosition&,
            gsl::span<const NumberT> tuning)
        {
            return InternalNormalizedDistance(a, b, tuning);
        }

        static EgocentricRelativeWristPosition Lerp(
            const EgocentricRelativeWristPosition& a,
            const EgocentricRelativeWristPosition& b,
            NumberT t)
        {
            VectorT aVec{ a.m_egocentricRelativeWristPosition };
            VectorT bVec{ b.m_egocentricRelativeWristPosition };
            EgocentricRelativeWristPosition result{};
            result.m_egocentricRelativeWristPosition = (static_cast<NumberT>(1) - t) * aVec + t * bVec;
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
            auto sequences = createDescriptorSequencesFromExamples<EgocentricRelativeWristPosition>(examples, DEFAULT_TUNING, tryCreate);
            auto extendedSequences = createDescriptorSequencesFromExamples<EgocentricRelativeWristPosition>(examples, DEFAULT_TUNING, tryCreate, 1.);

            auto tuning = DEFAULT_TUNING;
            NumberT maxAverageConnectionCost = 0;
            for (const auto& extendedSequence : extendedSequences)
            {
                for (const auto& sequence : sequences)
                {
                    auto absoluteDistFn = [](const auto& a, const auto& b) {
                        return InternalRawDistance(a, b);
                        };
                    auto deltaDistFn = [](const auto&, const auto&, const auto&, const auto&) -> NumberT {
                        return 0;
                        };
                    auto result = DynamicTimeWarping::Match<const EgocentricRelativeWristPosition>(extendedSequence, sequence, absoluteDistFn, deltaDistFn);
                    maxAverageConnectionCost = std::max<NumberT>(result.MaxConnectionCost / result.Connections, maxAverageConnectionCost);
                }
            }
            tuning[0] = std::max(maxAverageConnectionCost / IDENTICALITY_THRESHOLD, DEFAULT_TUNING[0]);
            return tuning;
        }

        template<bool NormalizeDistance>
        static auto Analyze(
            gsl::span<const EgocentricRelativeWristPosition> target,
            gsl::span<const EgocentricRelativeWristPosition> query,
            gsl::span<const NumberT> tuning)
        {
            std::array<AnalysisT, DEFAULT_TUNING.size()> results{};
            results[0] = { ANALYSIS_DIMENSION_NAME, IDENTICALITY_THRESHOLD, tuning.front(), {} };
            auto& rows = std::get<3>(results[0]);
            auto absoluteDistFn = [tuning](const auto& a, const auto& b) {
                if constexpr (NormalizeDistance)
                {
                    return InternalNormalizedDistance(a, b, tuning);
                }
                else
                {
                    return InternalRawDistance(a, b);
                }
                };
            auto deltaDistFn = [](const auto&, const auto&, const auto&, const auto&) -> NumberT {
                return 0;
                };
            auto rowsCallback = [&rows](std::vector<DynamicTimeWarping::MatchResult<NumberT>> row) { rows.push_back(std::move(row)); };
            DynamicTimeWarping::Match<const EgocentricRelativeWristPosition, decltype(absoluteDistFn), decltype(deltaDistFn), NumberT, true, decltype(rowsCallback)>(target, query, absoluteDistFn, deltaDistFn, 0, rowsCallback);
            return results;
        }

        EgocentricRelativeWristPosition() = default;

    private:
        static inline constexpr NumberT IDENTICALITY_THRESHOLD{ 0.02 };
        static inline constexpr auto normalizeDistance{ createDistanceNormalizationFunction(IDENTICALITY_THRESHOLD) };
        trivial::Point m_egocentricRelativeWristPosition{};

        EgocentricRelativeWristPosition(const InputSample& sample, const InputSample&)
        {
            auto& leftWristPose = sample.LeftWristPose.value();
            auto& rightWristPose = sample.RightWristPose.value();
            auto& hmdPose = sample.HmdPose.value();

            auto ets = EgocentricTemporalSpace::getPose(leftWristPose.translation(), hmdPose);
            m_egocentricRelativeWristPosition = ets.inverse() * rightWristPose.translation();
        }

        static NumberT InternalRawDistance(const EgocentricRelativeWristPosition& a, const EgocentricRelativeWristPosition& b)
        {
            return a.m_egocentricRelativeWristPosition.distance(b.m_egocentricRelativeWristPosition);
        }

        static NumberT InternalNormalizedDistance(const EgocentricRelativeWristPosition& a, const EgocentricRelativeWristPosition& b, gsl::span<const NumberT> tuning)
        {
            return normalizeDistance(InternalRawDistance(a, b), tuning.front());
        }
    };
}

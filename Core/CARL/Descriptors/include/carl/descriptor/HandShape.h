/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include <carl/descriptor/DescriptorUtils.h>
#include <carl/descriptor/Trivial.h>
#include <carl/descriptor/SequenceOperations.h>

namespace carl::descriptor
{
    template<Handedness Handedness>
    class HandShape
    {
        static constexpr std::array<size_t, 5> JOINTS{
            static_cast<size_t>(InputSample::Joint::XR_HAND_JOINT_THUMB_TIP_EXT),
            static_cast<size_t>(InputSample::Joint::XR_HAND_JOINT_INDEX_TIP_EXT),
            static_cast<size_t>(InputSample::Joint::XR_HAND_JOINT_MIDDLE_TIP_EXT),
            static_cast<size_t>(InputSample::Joint::XR_HAND_JOINT_RING_TIP_EXT),
            static_cast<size_t>(InputSample::Joint::XR_HAND_JOINT_LITTLE_TIP_EXT),
        };
        static constexpr std::array<const char*, JOINTS.size()> ANALYSIS_DIMENSION_NAMES{
            "Thumb",
            "Index Finger",
            "Middle Finger",
            "Ring Finger",
            "Little Finger",
        };
        static constexpr size_t NORMALIZATION_JOINT{
            static_cast<size_t>(InputSample::Joint::XR_HAND_JOINT_INDEX_PROXIMAL_EXT) };

    public:
        static constexpr std::array<NumberT, JOINTS.size()> DEFAULT_TUNING{ 1., 1., 1., 1., 1. };

        static std::optional<HandShape> TryCreate(
            const InputSample& sample,
            const InputSample& priorSample)
        {
            auto optionalWristPose = getWristPose<Handedness>(sample);
            if constexpr (Handedness == Handedness::LeftHanded)
            {
                if (optionalWristPose.has_value() && sample.LeftHandJointPoses.has_value())
                {
                    return HandShape{ sample, optionalWristPose.value() };
                }
            }
            else if constexpr (Handedness == Handedness::RightHanded)
            {
                if (optionalWristPose.has_value() && sample.RightHandJointPoses.has_value())
                {
                    return HandShape{ sample, optionalWristPose.value() };
                }
            }
            return {};
        }

        static NumberT Distance(
            const HandShape& a,
            const HandShape&,
            const HandShape& b,
            const HandShape&,
            gsl::span<const NumberT> tuning)
        {
            NumberT distance = 0;
            for (size_t idx = 0; idx < a.m_positions.size(); ++idx) {
                distance += InternalNormalizedDistance(idx, a, b, tuning);
            }
            return distance;
        }

        static HandShape Lerp(const HandShape& a, const HandShape& b, NumberT t)
        {
            HandShape result{};
            for (size_t idx = 0; idx < JOINTS.size(); ++idx)
            {
                VectorT aVec{ a.m_positions[idx] };
                VectorT bVec{ b.m_positions[idx] };
                result.m_positions[idx] = (static_cast<NumberT>(1) - t) * aVec + t * bVec;
            }
            return result;
        }

        static std::array<NumberT, DEFAULT_TUNING.size()> CalculateTuning(gsl::span<const action::Example> examples)
        {
            if (examples.size() < 2)
            {
                return DEFAULT_TUNING;
            }

            constexpr auto tryCreate = [](const InputSample& current, const InputSample& prior) { return TryCreate(current, prior); };
            auto sequences = createDescriptorSequencesFromExamples<HandShape<Handedness>>(examples, DEFAULT_TUNING, tryCreate);
            auto extendedSequences = createDescriptorSequencesFromExamples<HandShape<Handedness>>(examples, DEFAULT_TUNING, tryCreate, 1.);

            auto tuning = DEFAULT_TUNING;
            for (size_t idx = 0; idx < tuning.size(); ++idx)
            {
                NumberT maxAverageConnectionCost = 0;
                for (const auto& extendedSequence : extendedSequences)
                {
                    for (const auto& sequence : sequences)
                    {
                        auto distanceFunction = [idx](const auto& a, const auto&, const auto& b, const auto&) {
                            return InternalRawDistance(idx, a, b);
                            };
                        auto result = DynamicTimeWarping::Match<const HandShape<Handedness>>(extendedSequence, sequence, distanceFunction);
                        maxAverageConnectionCost = std::max<NumberT>(result.MaxConnectionCost / result.Connections, maxAverageConnectionCost);
                    }
                }
                tuning[idx] = std::max(maxAverageConnectionCost / IDENTICALITY_THRESHOLD, DEFAULT_TUNING[idx]);
            }
            return tuning;
        }

        template<bool NormalizeDistance>
        static auto Analyze(
            gsl::span<const HandShape<Handedness>> target,
            gsl::span<const HandShape<Handedness>> query,
            gsl::span<const NumberT> tuning)
        {
            std::array<AnalysisT, DEFAULT_TUNING.size()> results{};
            for (size_t idx = 0; idx < DEFAULT_TUNING.size(); ++idx)
            {
                results[idx] = { ANALYSIS_DIMENSION_NAMES[idx], IDENTICALITY_THRESHOLD, tuning[idx], {} };
                auto& rows = std::get<3>(results[idx]);
                auto distanceFunction = [idx, tuning](const auto& a, const auto&, const auto& b, const auto&) {
                    if constexpr (NormalizeDistance)
                    {
                        return InternalNormalizedDistance(idx, a, b, tuning);
                    }
                    else
                    {
                        return InternalRawDistance(idx, a, b);
                    }
                    };
                auto rowsCallback = [&rows](std::vector<DynamicTimeWarping::MatchResult<NumberT>> row) { rows.push_back(std::move(row)); };
                DynamicTimeWarping::Match<const HandShape<Handedness>, decltype(distanceFunction), NumberT, true, decltype(rowsCallback)>(target, query, distanceFunction, 0, rowsCallback);
            }
            return results;
        }

        HandShape() = default;

    private:
        static inline constexpr NumberT IDENTICALITY_THRESHOLD{ 0.005 };
        static inline constexpr auto normalizeDistance{ createDistanceNormalizationFunction(IDENTICALITY_THRESHOLD) };
        static inline constexpr NumberT CANONICAL_NORMALIZATION_LENGTH{ 0.1 };
        std::array<trivial::Point, JOINTS.size()> m_positions{};

        HandShape(const InputSample& sample, const TransformT& wristPose)
        {
            constexpr bool isLeftHanded = Handedness == Handedness::LeftHanded;
            const auto& jointPoses =
                isLeftHanded ? sample.LeftHandJointPoses.value() : sample.RightHandJointPoses.value();

            auto inverseWristPose = wristPose.inverse();

            NumberT normalizationFactor = CANONICAL_NORMALIZATION_LENGTH /
                (inverseWristPose * jointPoses[NORMALIZATION_JOINT].translation()).norm();
            for (size_t idx = 0; idx < JOINTS.size(); ++idx)
            {
                m_positions[idx] = normalizationFactor * (inverseWristPose * jointPoses[JOINTS[idx]].translation());
            }
        }

        static NumberT InternalRawDistance(size_t dimension, const HandShape& a, const HandShape& b)
        {
            return a.m_positions[dimension].distance(b.m_positions[dimension]);
        }

        static NumberT InternalNormalizedDistance(size_t dimension, const HandShape& a, const HandShape& b, gsl::span<const NumberT> tuning)
        {
            return normalizeDistance(InternalRawDistance(dimension, a, b), tuning[dimension]);
        }
    };
}
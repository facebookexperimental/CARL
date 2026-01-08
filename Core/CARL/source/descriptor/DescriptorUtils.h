/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include "carl/Example.h"
#include "DynamicTimeWarping.h"

#include <gsl/span>

namespace carl::descriptor
{
    enum class Handedness : uint64_t
    {
        LeftHanded,
        RightHanded,
    };

    template<Handedness Handedness>
    inline std::optional<TransformT> getWristPose(const InputSample& sample)
    {
        if constexpr (Handedness == Handedness::LeftHanded)
        {
            if (!sample.LeftWristPose.has_value())
            {
                return{};
            }
        }
        else
        {
            if (!sample.RightWristPose.has_value())
            {
                return{};
            }
        }

        constexpr size_t Y_AXIS_JOINT{
            static_cast<size_t>(InputSample::Joint::XR_HAND_JOINT_LITTLE_PROXIMAL_EXT) };
        constexpr size_t Z_AXIS_JOINT{
            static_cast<size_t>(InputSample::Joint::XR_HAND_JOINT_INDEX_PROXIMAL_EXT) };

        constexpr bool isLeftHanded = Handedness == Handedness::LeftHanded;
        const auto& wristPose =
            isLeftHanded ? sample.LeftWristPose.value() : sample.RightWristPose.value();

        const auto& jointPosesOptional =
            isLeftHanded ? sample.LeftHandJointPoses : sample.RightHandJointPoses;
        if (jointPosesOptional.has_value())
        {
            const auto& jointPoses = jointPosesOptional.value();
            auto zAxis = (jointPoses[Z_AXIS_JOINT].translation() - wristPose.translation()).normalized();
            auto yAxis = zAxis.cross((jointPoses[Y_AXIS_JOINT].translation() - wristPose.translation()).cross(zAxis)).normalized();
            // In some cases where OpenXR reports valid joint tracking but returns bogus data, all joint positions
            // will be tracked to identity, causing these calculations to produce 0-length axes. This allows NaNs to 
            // get into the descriptor sequences, wreaking havoc on recognition. The following check prevents us from
            // using these calculations when it is unsave to do so.
            if (yAxis.norm() < std::numeric_limits<NumberT>::epsilon() || zAxis.norm() < std::numeric_limits<NumberT>::epsilon())
            {
                return{};
            }
            return math::LookTransform(zAxis, yAxis, wristPose.translation());
        }
        else
        {
            // TODO: This implicit behavior change based on the presence or absence of joint data is a bug farm. 
            // Reevaluate the benefit provided by getWristPose, and consider eliminating it altogether if it's
            // not truly essential.
            return wristPose;
        }
    }

    constexpr NumberT NULL_TUNING{ 1000000 };
    constexpr auto createDistanceNormalizationFunction(NumberT identicalityThreshold, NumberT irreconcilabilityThreshold)
    {
        return [identicalityThreshold, irreconcilabilityThreshold](NumberT distance, NumberT tuning) {
            auto lowerBound = tuning * identicalityThreshold;
            auto upperBound = tuning * irreconcilabilityThreshold;
            return std::max<NumberT>(0, std::pow<NumberT>((distance - lowerBound) / (upperBound - lowerBound), 3));
            };
    }

    constexpr auto createDistanceNormalizationFunction(NumberT identicalityThreshold)
    {
        return createDistanceNormalizationFunction(identicalityThreshold, 2 * identicalityThreshold);
    }

    using AnalysisT = std::tuple<std::string, NumberT, NumberT, std::vector<std::vector<DynamicTimeWarping::MatchResult<NumberT>>>>;
}

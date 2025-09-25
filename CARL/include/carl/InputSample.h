/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include "carl/Serialization.h"

#include <Eigen/Geometry>

namespace carl
{
    /// <summary>
    /// A snapshot of all data relevant to CARL at a specific point in time.
    /// All 3D data described in an InputSample must be in the same reference
    /// space.
    /// </summary>
    struct InputSample
    {
        enum class Joint : uint64_t
        {
            XR_HAND_JOINT_PALM_EXT = 0,
            XR_HAND_JOINT_WRIST_EXT = 1,
            XR_HAND_JOINT_THUMB_METACARPAL_EXT = 2,
            XR_HAND_JOINT_THUMB_PROXIMAL_EXT = 3,
            XR_HAND_JOINT_THUMB_DISTAL_EXT = 4,
            XR_HAND_JOINT_THUMB_TIP_EXT = 5,
            XR_HAND_JOINT_INDEX_METACARPAL_EXT = 6,
            XR_HAND_JOINT_INDEX_PROXIMAL_EXT = 7,
            XR_HAND_JOINT_INDEX_INTERMEDIATE_EXT = 8,
            XR_HAND_JOINT_INDEX_DISTAL_EXT = 9,
            XR_HAND_JOINT_INDEX_TIP_EXT = 10,
            XR_HAND_JOINT_MIDDLE_METACARPAL_EXT = 11,
            XR_HAND_JOINT_MIDDLE_PROXIMAL_EXT = 12,
            XR_HAND_JOINT_MIDDLE_INTERMEDIATE_EXT = 13,
            XR_HAND_JOINT_MIDDLE_DISTAL_EXT = 14,
            XR_HAND_JOINT_MIDDLE_TIP_EXT = 15,
            XR_HAND_JOINT_RING_METACARPAL_EXT = 16,
            XR_HAND_JOINT_RING_PROXIMAL_EXT = 17,
            XR_HAND_JOINT_RING_INTERMEDIATE_EXT = 18,
            XR_HAND_JOINT_RING_DISTAL_EXT = 19,
            XR_HAND_JOINT_RING_TIP_EXT = 20,
            XR_HAND_JOINT_LITTLE_METACARPAL_EXT = 21,
            XR_HAND_JOINT_LITTLE_PROXIMAL_EXT = 22,
            XR_HAND_JOINT_LITTLE_INTERMEDIATE_EXT = 23,
            XR_HAND_JOINT_LITTLE_DISTAL_EXT = 24,
            XR_HAND_JOINT_LITTLE_TIP_EXT = 25,
            COUNT = 26,
        };

        InputSample() = default;
        InputSample(Deserialization& deserialization);

        void serialize(Serialization& serialization) const;

        /// <summary>
        /// Timestamp at which the observations within this InputSample were valid.
        /// </summary>
        double Timestamp{};

        std::optional<TransformT> HmdPose{};
        std::optional<TransformT> LeftWristPose{};
        std::optional<TransformT> RightWristPose{};
        std::optional<std::array<TransformT, static_cast<size_t>(Joint::COUNT)>> LeftHandJointPoses{};
        std::optional<std::array<TransformT, static_cast<size_t>(Joint::COUNT)>> RightHandJointPoses{};

        static InputSample Lerp(const InputSample& a, const InputSample& b, double t);
    };
}

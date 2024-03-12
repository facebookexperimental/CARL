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
            UNUSED_HandJointId_HandThumb0,
            ThumbFingerBase,
            UNUSED_HandJointId_HandThumb2,
            UNUSED_HandJointId_HandThumb3,
            ThumbFingerTip,
            IndexFingerBase,
            UNUSED_HandJointId_HandIndex2,
            UNUSED_HandJointId_HandIndex3,
            IndexFingerTip,
            MiddleFingerBase,
            UNUSED_HandJointId_HandMiddle2,
            UNUSED_HandJointId_HandMiddle3,
            MiddleFingerTip,
            RingFingerBase,
            UNUSED_HandJointId_HandRing2,
            UNUSED_HandJointId_HandRing3,
            RingFingerTip,
            UNUSED_HandJointId_HandPinky0,
            LittleFingerBase,
            UNUSED_HandJointId_HandPinky2,
            UNUSED_HandJointId_HandPinky3,
            LittleFingerTip,
            COUNT
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

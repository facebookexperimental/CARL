/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include <cstdint>

namespace carl::action
{
    enum class ActionType : uint64_t
    {
        LeftHandPose = 0,
        LeftHandGesture = 1,
        RightHandPose = 2,
        RightHandGesture = 3,
        TwoHandGesture = 4,
        LeftControllerGesture = 5,
        RightControllerGesture = 6,
        TwoControllerGesture = 7,
        LeftWristTrajectory = 8,
        RightWristTrajectory = 9,
        LeftHandShape = 10,
        RightHandShape = 11,
        Custom = ~uint64_t{0x0},
    };
}

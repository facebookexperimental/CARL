/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include <carl/ActionTypeSet.h>
#include <carl/Descriptor.h>

namespace carl::action
{
    using ActionTypes = ActionTypeSet<
        ActionTypeDefinition<ActionType::LeftHandShape, descriptor::HandShape<descriptor::Handedness::LeftHanded>, true>,
        ActionTypeDefinition<ActionType::LeftHandPose, descriptor::HandPose<descriptor::Handedness::LeftHanded>, true>,
        ActionTypeDefinition<ActionType::LeftHandGesture, descriptor::HandGesture<descriptor::Handedness::LeftHanded>>,
        ActionTypeDefinition<ActionType::LeftControllerGesture, descriptor::ControllerGesture<descriptor::Handedness::LeftHanded>>,
        ActionTypeDefinition<ActionType::LeftWristTrajectory, descriptor::WristTrajectory<descriptor::Handedness::LeftHanded>>,
        ActionTypeDefinition<ActionType::RightHandShape, descriptor::HandShape<descriptor::Handedness::RightHanded>, true>,
        ActionTypeDefinition<ActionType::RightHandPose, descriptor::HandPose<descriptor::Handedness::RightHanded>, true>,
        ActionTypeDefinition<ActionType::RightHandGesture, descriptor::HandGesture<descriptor::Handedness::RightHanded>>,
        ActionTypeDefinition<ActionType::RightControllerGesture, descriptor::ControllerGesture<descriptor::Handedness::RightHanded>>,
        ActionTypeDefinition<ActionType::RightWristTrajectory, descriptor::WristTrajectory<descriptor::Handedness::RightHanded>>,
        ActionTypeDefinition<ActionType::TwoHandGesture, descriptor::TwoHandGesture>,
        ActionTypeDefinition<ActionType::TwoControllerGesture, descriptor::TwoControllerGesture>
    >;
}

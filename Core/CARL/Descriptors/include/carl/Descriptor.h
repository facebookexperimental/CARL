/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include <carl/descriptor/CombinedDescriptor.h>
#include <carl/descriptor/Custom.h>
#include <carl/descriptor/EgocentricWristDisplacement.h>
#include <carl/descriptor/ControllerState.h>
#include <carl/descriptor/EgocentricWristTranslation.h>
#include <carl/descriptor/EgocentricRelativeWristPosition.h>
#include <carl/descriptor/TimePoint.h>
#include <carl/descriptor/EgocentricWristOrientation.h>
#include <carl/descriptor/WristRotation.h>
#include <carl/descriptor/HandShape.h>
#include <carl/descriptor/TimestampedDescriptor.h>

namespace carl::descriptor
{
    template<Handedness Handedness>
    using HandPose = CombinedDescriptorT<HandShape<Handedness>, EgocentricWristOrientation<Handedness>>;

    template<Handedness Handedness>
    using WristTrajectory = CombinedDescriptorT<TimePoint, EgocentricWristDisplacement<Handedness>>;

    template<Handedness Handedness>
    using HandGesture = CombinedDescriptorT<HandPose<Handedness>, WristTrajectory<Handedness>, WristRotation<Handedness>, EgocentricWristTranslation<Handedness>>;

    using TwoHandGesture = CombinedDescriptorT<HandGesture<Handedness::LeftHanded>, HandGesture<Handedness::RightHanded>, EgocentricRelativeWristPosition>;

    template<Handedness Handedness>
    using ControllerGesture = CombinedDescriptorT<ControllerState<Handedness>, WristTrajectory<Handedness>, EgocentricWristOrientation<Handedness>, WristRotation<Handedness>, EgocentricWristTranslation<Handedness>>;

    using TwoControllerGesture = CombinedDescriptorT<ControllerGesture<Handedness::LeftHanded>, ControllerGesture<Handedness::RightHanded>, EgocentricRelativeWristPosition>;
}

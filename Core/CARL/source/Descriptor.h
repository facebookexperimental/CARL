/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include "descriptor/CombinedDescriptor.h"
#include "descriptor/Custom.h"
#include "descriptor/EgocentricWristDisplacement.h"
#include "descriptor/ControllerState.h"
#include "descriptor/EgocentricWristTranslation.h"
#include "descriptor/EgocentricRelativeWristPosition.h"
#include "descriptor/TimePoint.h"
#include "descriptor/EgocentricWristOrientation.h"
#include "descriptor/WristRotation.h"
#include "descriptor/HandShape.h"
#include "descriptor/TimestampedDescriptor.h"

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

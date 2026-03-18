/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "InputSampleConversion.h"

#include <map>

namespace carl::capi
{
    static const std::map<carl_InputSample::HAND_JOINT, carl::InputSample::Joint> openXrToCarlJointMap
    {
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_THUMB_METACARPAL_EXT, carl::InputSample::Joint::XR_HAND_JOINT_THUMB_METACARPAL_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_THUMB_PROXIMAL_EXT, carl::InputSample::Joint::XR_HAND_JOINT_THUMB_PROXIMAL_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_THUMB_DISTAL_EXT, carl::InputSample::Joint::XR_HAND_JOINT_THUMB_DISTAL_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_THUMB_TIP_EXT, carl::InputSample::Joint::XR_HAND_JOINT_THUMB_TIP_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_INDEX_PROXIMAL_EXT, carl::InputSample::Joint::XR_HAND_JOINT_INDEX_PROXIMAL_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_INDEX_INTERMEDIATE_EXT, carl::InputSample::Joint::XR_HAND_JOINT_INDEX_INTERMEDIATE_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_INDEX_DISTAL_EXT, carl::InputSample::Joint::XR_HAND_JOINT_INDEX_DISTAL_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_INDEX_TIP_EXT, carl::InputSample::Joint::XR_HAND_JOINT_INDEX_TIP_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_MIDDLE_PROXIMAL_EXT, carl::InputSample::Joint::XR_HAND_JOINT_MIDDLE_PROXIMAL_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_MIDDLE_INTERMEDIATE_EXT, carl::InputSample::Joint::XR_HAND_JOINT_MIDDLE_INTERMEDIATE_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_MIDDLE_DISTAL_EXT, carl::InputSample::Joint::XR_HAND_JOINT_MIDDLE_DISTAL_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_MIDDLE_TIP_EXT, carl::InputSample::Joint::XR_HAND_JOINT_MIDDLE_TIP_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_RING_PROXIMAL_EXT, carl::InputSample::Joint::XR_HAND_JOINT_RING_PROXIMAL_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_RING_INTERMEDIATE_EXT, carl::InputSample::Joint::XR_HAND_JOINT_RING_INTERMEDIATE_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_RING_DISTAL_EXT, carl::InputSample::Joint::XR_HAND_JOINT_RING_DISTAL_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_RING_TIP_EXT, carl::InputSample::Joint::XR_HAND_JOINT_RING_TIP_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_LITTLE_METACARPAL_EXT, carl::InputSample::Joint::XR_HAND_JOINT_LITTLE_METACARPAL_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_LITTLE_PROXIMAL_EXT, carl::InputSample::Joint::XR_HAND_JOINT_LITTLE_PROXIMAL_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_LITTLE_INTERMEDIATE_EXT, carl::InputSample::Joint::XR_HAND_JOINT_LITTLE_INTERMEDIATE_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_LITTLE_DISTAL_EXT, carl::InputSample::Joint::XR_HAND_JOINT_LITTLE_DISTAL_EXT},
        {carl_InputSample::HAND_JOINT::XR_HAND_JOINT_LITTLE_TIP_EXT, carl::InputSample::Joint::XR_HAND_JOINT_LITTLE_TIP_EXT},
    };

    carl::InputSample convert(const carl_InputSample& input)
    {
        carl::InputSample sample{};

        sample.Timestamp = input.Timestamp;

        constexpr auto cvt = [](const carl_InputSample::OptionalTransform& inT, carl::TransformT& outT) {
            outT.fromPositionOrientationScale(
                carl::VectorT
                {
                    static_cast<carl::NumberT>(inT.Position.X),
                    static_cast<carl::NumberT>(inT.Position.Y),
                    static_cast<carl::NumberT>(inT.Position.Z)
                },
                carl::QuaternionT
                {
                    static_cast<carl::NumberT>(inT.Orientation.W),
                    static_cast<carl::NumberT>(inT.Orientation.X),
                    static_cast<carl::NumberT>(inT.Orientation.Y),
                    static_cast<carl::NumberT>(inT.Orientation.Z)
                },
                carl::UNIT_SCALE);
        };

        constexpr auto cvtOptional = [cvt](const carl_InputSample::OptionalTransform& inT) {
            std::optional<carl::TransformT> outT{};
            if (inT.Valid)
            {
                outT.emplace();
                cvt(inT, *outT);
            }
            return outT;
        };

        sample.HmdPose = cvtOptional(input.HmdPose);
        sample.LeftWristPose = cvtOptional(input.LeftWristPose);
        sample.RightWristPose = cvtOptional(input.RightWristPose);
        // TODO: Fix mismatched joint count mapping so that palms, etc. stop mapping to nothing and causing weirdness.
        if (input.LeftHandJointPoses[2].Valid)
        {
            sample.LeftHandJointPoses.emplace();
            for (size_t idx = 0; idx < carl_InputSample::HAND_JOINT::COUNT; ++idx)
            {
                auto found = openXrToCarlJointMap.find(static_cast<carl_InputSample::HAND_JOINT>(idx));
                if (found != openXrToCarlJointMap.end())
                {
                    cvt(input.LeftHandJointPoses[idx], sample.LeftHandJointPoses->at(static_cast<size_t>(found->second)));
                }
            }
        }
        if (input.RightHandJointPoses[2].Valid)
        {
            sample.RightHandJointPoses.emplace();
            for (size_t idx = 0; idx < carl_InputSample::HAND_JOINT::COUNT; ++idx)
            {
                auto found = openXrToCarlJointMap.find(static_cast<carl_InputSample::HAND_JOINT>(idx));
                if (found != openXrToCarlJointMap.end())
                {
                    cvt(input.RightHandJointPoses[idx], sample.RightHandJointPoses->at(static_cast<size_t>(found->second)));
                }
            }
        }
        if (input.LeftControllerState.Valid)
        {
            sample.LeftControllerInput = 
                std::array<carl::NumberT, static_cast<size_t>(carl::InputSample::ControllerInput::COUNT)>{
                static_cast<carl::NumberT>(input.LeftControllerState.PrimaryClick),
                static_cast<carl::NumberT>(input.LeftControllerState.SecondaryClick),
                static_cast<carl::NumberT>(input.LeftControllerState.ThumbstickX),
                static_cast<carl::NumberT>(input.LeftControllerState.ThumbstickY),
                static_cast<carl::NumberT>(input.LeftControllerState.ThumbstickClick),
                static_cast<carl::NumberT>(input.LeftControllerState.SqueezeValue),
                static_cast<carl::NumberT>(input.LeftControllerState.TriggerValue),
            };
        }
        if (input.RightControllerState.Valid)
        {
            sample.RightControllerInput =
                std::array<carl::NumberT, static_cast<size_t>(carl::InputSample::ControllerInput::COUNT)>{
                static_cast<carl::NumberT>(input.RightControllerState.PrimaryClick),
                static_cast<carl::NumberT>(input.RightControllerState.SecondaryClick),
                static_cast<carl::NumberT>(input.RightControllerState.ThumbstickX),
                static_cast<carl::NumberT>(input.RightControllerState.ThumbstickY),
                static_cast<carl::NumberT>(input.RightControllerState.ThumbstickClick),
                static_cast<carl::NumberT>(input.RightControllerState.SqueezeValue),
                static_cast<carl::NumberT>(input.RightControllerState.TriggerValue),
            };
        }
        
        return sample;
    }

    carl_InputSample convert(const carl::InputSample& input)
    {
        carl_InputSample sample{};

        sample.Timestamp = input.Timestamp;

        constexpr auto cvt = [](const carl::TransformT& inT, carl_InputSample::OptionalTransform& outT) {
            carl::VectorT position{ inT.translation() };
            carl::QuaternionT orientation{ inT.rotation() };
            outT.Valid = true;
            outT.Position.X = position.x();
            outT.Position.Y = position.y();
            outT.Position.Z = position.z();
            outT.Orientation.W = orientation.w();
            outT.Orientation.X = orientation.x();
            outT.Orientation.Y = orientation.y();
            outT.Orientation.Z = orientation.z();
        };

        constexpr auto cvtOptional = [cvt](const std::optional<carl::TransformT>& inT) {
            carl_InputSample::OptionalTransform outT{};
            if (inT.has_value())
            {
                outT.Valid = true;
                cvt(*inT, outT);
            }
            return outT;
        };

        sample.HmdPose = cvtOptional(input.HmdPose);
        sample.LeftWristPose = cvtOptional(input.LeftWristPose);
        sample.RightWristPose = cvtOptional(input.RightWristPose);
        if (input.LeftHandJointPoses.has_value())
        {
            for (size_t idx = 0; idx < carl_InputSample::HAND_JOINT::COUNT; ++idx)
            {
                auto found = openXrToCarlJointMap.find(static_cast<carl_InputSample::HAND_JOINT>(idx));
                if (found != openXrToCarlJointMap.end())
                {
                    cvt(input.LeftHandJointPoses->at(idx), sample.LeftHandJointPoses[static_cast<size_t>(found->second)]);
                }
            }
        }
        if (input.RightHandJointPoses.has_value())
        {
            for (size_t idx = 0; idx < carl_InputSample::HAND_JOINT::COUNT; ++idx)
            {
                auto found = openXrToCarlJointMap.find(static_cast<carl_InputSample::HAND_JOINT>(idx));
                if (found != openXrToCarlJointMap.end())
                {
                    cvt(input.RightHandJointPoses->at(idx), sample.RightHandJointPoses[static_cast<size_t>(found->second)]);
                }
            }
        }
        if (input.LeftControllerInput.has_value())
        {
            sample.LeftControllerState.Valid = true;
            sample.LeftControllerState.PrimaryClick = input.LeftControllerInput->at(static_cast<size_t>(carl::InputSample::ControllerInput::XR_KHR_generic_controller_PRIMARY_CLICK));
            sample.LeftControllerState.SecondaryClick = input.LeftControllerInput->at(static_cast<size_t>(carl::InputSample::ControllerInput::XR_KHR_generic_controller_SECONDARY_CLICK));
            sample.LeftControllerState.ThumbstickX = input.LeftControllerInput->at(static_cast<size_t>(carl::InputSample::ControllerInput::XR_KHR_generic_controller_THUMBSTICK_X));
            sample.LeftControllerState.ThumbstickY = input.LeftControllerInput->at(static_cast<size_t>(carl::InputSample::ControllerInput::XR_KHR_generic_controller_THUMBSTICK_Y));
            sample.LeftControllerState.ThumbstickClick = input.LeftControllerInput->at(static_cast<size_t>(carl::InputSample::ControllerInput::XR_KHR_generic_controller_THUMBSTICK_CLICK));
            sample.LeftControllerState.SqueezeValue = input.LeftControllerInput->at(static_cast<size_t>(carl::InputSample::ControllerInput::XR_KHR_generic_controller_SQUEEZE_VALUE));
            sample.LeftControllerState.TriggerValue = input.LeftControllerInput->at(static_cast<size_t>(carl::InputSample::ControllerInput::XR_KHR_generic_controller_TRIGGER_VALUE));
        }
        if (input.RightControllerInput.has_value())
        {
            sample.RightControllerState.Valid = true;
            sample.RightControllerState.PrimaryClick = input.RightControllerInput->at(static_cast<size_t>(carl::InputSample::ControllerInput::XR_KHR_generic_controller_PRIMARY_CLICK));
            sample.RightControllerState.SecondaryClick = input.RightControllerInput->at(static_cast<size_t>(carl::InputSample::ControllerInput::XR_KHR_generic_controller_SECONDARY_CLICK));
            sample.RightControllerState.ThumbstickX = input.RightControllerInput->at(static_cast<size_t>(carl::InputSample::ControllerInput::XR_KHR_generic_controller_THUMBSTICK_X));
            sample.RightControllerState.ThumbstickY = input.RightControllerInput->at(static_cast<size_t>(carl::InputSample::ControllerInput::XR_KHR_generic_controller_THUMBSTICK_Y));
            sample.RightControllerState.ThumbstickClick = input.RightControllerInput->at(static_cast<size_t>(carl::InputSample::ControllerInput::XR_KHR_generic_controller_THUMBSTICK_CLICK));
            sample.RightControllerState.SqueezeValue = input.RightControllerInput->at(static_cast<size_t>(carl::InputSample::ControllerInput::XR_KHR_generic_controller_SQUEEZE_VALUE));
            sample.RightControllerState.TriggerValue = input.RightControllerInput->at(static_cast<size_t>(carl::InputSample::ControllerInput::XR_KHR_generic_controller_TRIGGER_VALUE));
        }

        return sample;
    }
}

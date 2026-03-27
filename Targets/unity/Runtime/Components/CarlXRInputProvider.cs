/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

using Carl.Native;
using UnityEngine;
using UnityEngine.XR;
using UnityEngine.XR.Hands;
using System.Collections.Generic;

namespace Carl
{
    /// <summary>
    /// Reads XR input (hand tracking, HMD, controllers) each frame
    /// and feeds it to the CARL session as InputSamples.
    /// Requires a <see cref="CarlSessionManager"/> in the scene.
    /// </summary>
    [AddComponentMenu("CARL/XR Input Provider")]
    public class CarlXRInputProvider : MonoBehaviour
    {
        [Tooltip("Capture hand tracking joint data.")]
        [SerializeField] bool captureHands = true;

        [Tooltip("Capture VR controller state (buttons, triggers, thumbsticks).")]
        [SerializeField] bool captureControllers = true;

        [Tooltip("Capture head-mounted display (HMD) pose.")]
        [SerializeField] bool captureHmd = true;

        [Tooltip("Work around a Unity XRHandSubsystem bug where the Palm (index 0) and " +
                 "Wrist (index 1) joint positions are swapped. When enabled, the Palm joint " +
                 "is used as the wrist pose source. Disable this if your XR runtime reports " +
                 "these joints correctly.")]
        [SerializeField] bool palmWristSwapWorkaround = true;

        CarlInputSampleInterop _sample;
        readonly List<InputDevice> _deviceBuffer = new List<InputDevice>();
        XRHandSubsystem _handSubsystem;

        /// <summary>
        /// Provides public access to the most recently built sample,
        /// enabling other components (e.g., CarlGestureRecorder) to capture it.
        /// </summary>
        public ref CarlInputSampleInterop CurrentSample => ref _sample;

        void Update()
        {
            if (CarlSessionManager.Instance == null || CarlSessionManager.Instance.Session == null)
                return;

            _sample = default;
            _sample.Timestamp = Time.realtimeSinceStartupAsDouble;

            if (captureHmd)
                CaptureHmd();
            if (captureHands)
                CaptureHands();
            if (captureControllers)
                CaptureControllers();

            CarlSessionManager.Instance.Session.AddInputSample(ref _sample);
        }

        void CaptureHmd()
        {
            _deviceBuffer.Clear();
            InputDevices.GetDevicesAtXRNode(XRNode.Head, _deviceBuffer);
            if (_deviceBuffer.Count > 0 && _deviceBuffer[0].TryGetFeatureValue(CommonUsages.devicePosition, out Vector3 hmdPos)
                                      && _deviceBuffer[0].TryGetFeatureValue(CommonUsages.deviceRotation, out Quaternion hmdRot))
            {
                _sample.HmdPose = CarlCoordinateConversion.ToCarlTransform(hmdPos, hmdRot);
            }
        }

        void CaptureHands()
        {
            if (_handSubsystem == null || !_handSubsystem.running)
            {
                var subsystems = new List<XRHandSubsystem>();
                SubsystemManager.GetSubsystems(subsystems);
                _handSubsystem = subsystems.Count > 0 ? subsystems[0] : null;
            }

            if (_handSubsystem == null)
                return;

            CaptureHand(_handSubsystem.leftHand, ref _sample.LeftWristPose, ref _sample.LeftHandJointPoses, palmWristSwapWorkaround);
            CaptureHand(_handSubsystem.rightHand, ref _sample.RightWristPose, ref _sample.RightHandJointPoses, palmWristSwapWorkaround);
        }

        static void CaptureHand(XRHand hand, ref CarlOptionalTransform wristPose, ref CarlHandJointPoses jointPoses, bool palmWristSwap)
        {
            if (!hand.isTracked)
                return;

            // Which joint to use as the wrist pose source. Unity's XRHandSubsystem has a
            // known bug where Palm (index 0) and Wrist (index 1) positions are swapped:
            // on-device testing confirmed that a marker at Unity's "Palm" tracks with the
            // anatomical wrist, and vice versa. When palmWristSwap is enabled, we read
            // from Palm to get the correct anatomical wrist position.
            var wristJoint = palmWristSwap ? HandJoint.Palm : HandJoint.Wrist;

            for (int i = 0; i < CarlHandJointPoses.Count; i++)
            {
                var jointId = XRHandJointIDUtility.FromIndex(i);
                var joint = hand.GetJoint(jointId);
                if (joint.TryGetPose(out Pose pose))
                {
                    jointPoses[i] = CarlCoordinateConversion.ToCarlTransform(pose);

                    if (i == (int)wristJoint)
                    {
                        wristPose = CarlCoordinateConversion.ToCarlTransform(pose);
                    }
                }
            }
        }

        void CaptureControllers()
        {
            CaptureController(XRNode.LeftHand, ref _sample.LeftControllerState, ref _sample.LeftWristPose, _deviceBuffer);
            CaptureController(XRNode.RightHand, ref _sample.RightControllerState, ref _sample.RightWristPose, _deviceBuffer);
        }

        static void CaptureController(XRNode node, ref CarlOptionalControllerState state, ref CarlOptionalTransform wristPose, List<InputDevice> deviceBuffer)
        {
            deviceBuffer.Clear();
            InputDevices.GetDevicesAtXRNode(node, deviceBuffer);
            if (deviceBuffer.Count == 0)
                return;

            var device = deviceBuffer[0];
            if (!device.isValid)
                return;

            if (wristPose.Valid == 0
                && device.TryGetFeatureValue(CommonUsages.devicePosition, out Vector3 pos)
                && device.TryGetFeatureValue(CommonUsages.deviceRotation, out Quaternion rot))
            {
                wristPose = CarlCoordinateConversion.ToCarlTransform(pos, rot);
            }

            state.IsValid = true;

            if (device.TryGetFeatureValue(CommonUsages.primaryButton, out bool primaryBtn))
                state.PrimaryClick = primaryBtn ? 1.0 : 0.0;
            if (device.TryGetFeatureValue(CommonUsages.secondaryButton, out bool secondaryBtn))
                state.SecondaryClick = secondaryBtn ? 1.0 : 0.0;
            if (device.TryGetFeatureValue(CommonUsages.primary2DAxis, out Vector2 thumbstick))
            {
                state.ThumbstickX = thumbstick.x;
                state.ThumbstickY = thumbstick.y;
            }
            if (device.TryGetFeatureValue(CommonUsages.primary2DAxisClick, out bool thumbstickClick))
                state.ThumbstickClick = thumbstickClick ? 1.0 : 0.0;
            if (device.TryGetFeatureValue(CommonUsages.grip, out float grip))
                state.SqueezeValue = grip;
            if (device.TryGetFeatureValue(CommonUsages.trigger, out float trigger))
                state.TriggerValue = trigger;
        }
    }
}

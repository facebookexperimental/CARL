/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

using System;
using System.Runtime.InteropServices;

namespace Carl.Native
{
    [StructLayout(LayoutKind.Sequential)]
    public struct CarlVector3d
    {
        public double X;
        public double Y;
        public double Z;
    }

    [StructLayout(LayoutKind.Sequential)]
    public struct CarlQuaterniond
    {
        public double W;
        public double X;
        public double Y;
        public double Z;
    }

    [StructLayout(LayoutKind.Sequential)]
    public struct CarlOptionalTransform
    {
        public byte Valid;
        byte _pad0, _pad1, _pad2, _pad3, _pad4, _pad5, _pad6;
        public CarlVector3d Position;
        public CarlQuaterniond Orientation;

        public bool IsValid
        {
            get => Valid != 0;
            set => Valid = value ? (byte)1 : (byte)0;
        }
    }

    [StructLayout(LayoutKind.Sequential)]
    public struct CarlOptionalControllerState
    {
        public byte Valid;
        byte _pad0, _pad1, _pad2, _pad3, _pad4, _pad5, _pad6;
        public double PrimaryClick;
        public double SecondaryClick;
        public double ThumbstickX;
        public double ThumbstickY;
        public double ThumbstickClick;
        public double SqueezeValue;
        public double TriggerValue;

        public bool IsValid
        {
            get => Valid != 0;
            set => Valid = value ? (byte)1 : (byte)0;
        }
    }

    [StructLayout(LayoutKind.Sequential)]
    public struct CarlHandJointPoses
    {
        public CarlOptionalTransform Joint0;
        public CarlOptionalTransform Joint1;
        public CarlOptionalTransform Joint2;
        public CarlOptionalTransform Joint3;
        public CarlOptionalTransform Joint4;
        public CarlOptionalTransform Joint5;
        public CarlOptionalTransform Joint6;
        public CarlOptionalTransform Joint7;
        public CarlOptionalTransform Joint8;
        public CarlOptionalTransform Joint9;
        public CarlOptionalTransform Joint10;
        public CarlOptionalTransform Joint11;
        public CarlOptionalTransform Joint12;
        public CarlOptionalTransform Joint13;
        public CarlOptionalTransform Joint14;
        public CarlOptionalTransform Joint15;
        public CarlOptionalTransform Joint16;
        public CarlOptionalTransform Joint17;
        public CarlOptionalTransform Joint18;
        public CarlOptionalTransform Joint19;
        public CarlOptionalTransform Joint20;
        public CarlOptionalTransform Joint21;
        public CarlOptionalTransform Joint22;
        public CarlOptionalTransform Joint23;
        public CarlOptionalTransform Joint24;
        public CarlOptionalTransform Joint25;

        public unsafe CarlOptionalTransform this[int index]
        {
            get
            {
                if ((uint)index >= 26)
                    throw new ArgumentOutOfRangeException(nameof(index));
                fixed (CarlOptionalTransform* p = &Joint0)
                {
                    return p[index];
                }
            }
            set
            {
                if ((uint)index >= 26)
                    throw new ArgumentOutOfRangeException(nameof(index));
                fixed (CarlOptionalTransform* p = &Joint0)
                {
                    p[index] = value;
                }
            }
        }

        public const int Count = 26;
    }

    [StructLayout(LayoutKind.Sequential)]
    public struct CarlInputSampleInterop
    {
        public double Timestamp;
        public CarlOptionalTransform HmdPose;
        public CarlOptionalTransform LeftWristPose;
        public CarlOptionalTransform RightWristPose;
        public CarlHandJointPoses LeftHandJointPoses;
        public CarlHandJointPoses RightHandJointPoses;
        public CarlOptionalControllerState LeftControllerState;
        public CarlOptionalControllerState RightControllerState;

        public const int ExpectedSize = 3656;
    }

    internal static class InteropSizeValidation
    {
#if UNITY_5_3_OR_NEWER
        [UnityEngine.RuntimeInitializeOnLoadMethod(UnityEngine.RuntimeInitializeLoadType.SubsystemRegistration)]
#endif
        static void ValidateStructSizes()
        {
            int actualSize;
            unsafe { actualSize = sizeof(CarlInputSampleInterop); }
            if (actualSize != CarlInputSampleInterop.ExpectedSize)
            {
                throw new InvalidOperationException(
                    $"CARL interop struct size mismatch: expected {CarlInputSampleInterop.ExpectedSize} bytes, " +
                    $"got {actualSize} bytes. This indicates an ABI incompatibility between the managed and native layers.");
            }
        }
    }
}

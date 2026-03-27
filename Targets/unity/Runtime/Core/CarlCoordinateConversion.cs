/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

using Carl.Native;
using UnityEngine;

namespace Carl
{
    /// <summary>
    /// Converts between Unity's left-handed Y-up coordinate system
    /// and CARL's right-handed Y-up (OpenXR) coordinate system.
    /// Conversion: negate Z for positions, negate Z and W for quaternion components.
    /// This matches the convention established in the original CARL ActionStudio prototype.
    /// </summary>
    public static class CarlCoordinateConversion
    {
        public static CarlVector3d ToCarlPosition(Vector3 unityPos)
        {
            return new CarlVector3d
            {
                X = unityPos.x,
                Y = unityPos.y,
                Z = -unityPos.z,
            };
        }

        public static CarlQuaterniond ToCarlRotation(Quaternion unityRot)
        {
            return new CarlQuaterniond
            {
                W = -unityRot.w,
                X = unityRot.x,
                Y = unityRot.y,
                Z = -unityRot.z,
            };
        }

        public static CarlOptionalTransform ToCarlTransform(Pose unityPose)
        {
            return new CarlOptionalTransform
            {
                Valid = 1,
                Position = ToCarlPosition(unityPose.position),
                Orientation = ToCarlRotation(unityPose.rotation),
            };
        }

        public static CarlOptionalTransform ToCarlTransform(Vector3 position, Quaternion rotation)
        {
            return new CarlOptionalTransform
            {
                Valid = 1,
                Position = ToCarlPosition(position),
                Orientation = ToCarlRotation(rotation),
            };
        }

        public static CarlOptionalTransform InvalidTransform()
        {
            return default;
        }

        public static Vector3 ToUnityPosition(CarlVector3d carlPos)
        {
            return new Vector3(
                (float)carlPos.X,
                (float)carlPos.Y,
                (float)-carlPos.Z);
        }

        public static Quaternion ToUnityRotation(CarlQuaterniond carlRot)
        {
            return new Quaternion(
                (float)carlRot.X,
                (float)carlRot.Y,
                (float)-carlRot.Z,
                (float)-carlRot.W);
        }

        public static Pose ToUnityPose(CarlOptionalTransform carlTransform)
        {
            return new Pose(
                ToUnityPosition(carlTransform.Position),
                ToUnityRotation(carlTransform.Orientation));
        }
    }
}

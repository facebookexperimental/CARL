/**
 * Shared utility functions and constants for the ImmersiveExperience layer.
 *
 * Responsibilities:
 *  - OPENXR_JOINT_MAPPINGS: maps CARL joint indices to Babylon.js WebXRHandJoint values.
 *  - populateInputSample: fills an ICarlInputSample from live WebXR camera and hand data.
 *  - applyJointSampleToMeshes: applies a recorded/live input sample to arrays of meshes,
 *    re-orienting joint poses from the sample's coordinate frame into the point-of-view frame.
 */
import { WebXRHandJoint, Quaternion, Matrix, Vector3, WebXRCamera, WebXRHand, AbstractMesh } from "@babylonjs/core";
import { ICarlInputSample, ICarlOptionalTransform } from "./carlInterfaces";

export const OPENXR_JOINT_MAPPINGS = [
    WebXRHandJoint.WRIST,                               // XR_HAND_JOINT_PALM_EXT
    WebXRHandJoint.WRIST,                               // XR_HAND_JOINT_WRIST_EXT
    WebXRHandJoint.THUMB_METACARPAL,                    // XR_HAND_JOINT_THUMB_METACARPAL_EXT
    WebXRHandJoint.THUMB_PHALANX_PROXIMAL,              // XR_HAND_JOINT_THUMB_PROXIMAL_EXT
    WebXRHandJoint.THUMB_PHALANX_DISTAL,                // XR_HAND_JOINT_THUMB_DISTAL_EXT
    WebXRHandJoint.THUMB_TIP,                           // XR_HAND_JOINT_THUMB_TIP_EXT
    WebXRHandJoint.INDEX_FINGER_METACARPAL,             // XR_HAND_JOINT_INDEX_METACARPAL_EXT
    WebXRHandJoint.INDEX_FINGER_PHALANX_PROXIMAL,       // XR_HAND_JOINT_INDEX_PROXIMAL_EXT
    WebXRHandJoint.INDEX_FINGER_PHALANX_INTERMEDIATE,   // XR_HAND_JOINT_INDEX_INTERMEDIATE_EXT
    WebXRHandJoint.INDEX_FINGER_PHALANX_DISTAL,         // XR_HAND_JOINT_INDEX_DISTAL_EXT
    WebXRHandJoint.INDEX_FINGER_TIP,                    // XR_HAND_JOINT_INDEX_TIP_EXT
    WebXRHandJoint.MIDDLE_FINGER_METACARPAL,            // XR_HAND_JOINT_MIDDLE_METACARPAL_EXT
    WebXRHandJoint.MIDDLE_FINGER_PHALANX_PROXIMAL,      // XR_HAND_JOINT_MIDDLE_PROXIMAL_EXT
    WebXRHandJoint.MIDDLE_FINGER_PHALANX_INTERMEDIATE,  // XR_HAND_JOINT_MIDDLE_INTERMEDIATE_EXT
    WebXRHandJoint.MIDDLE_FINGER_PHALANX_DISTAL,        // XR_HAND_JOINT_MIDDLE_DISTAL_EXT
    WebXRHandJoint.MIDDLE_FINGER_TIP,                   // XR_HAND_JOINT_MIDDLE_TIP_EXT
    WebXRHandJoint.RING_FINGER_METACARPAL,              // XR_HAND_JOINT_RING_METACARPAL_EXT
    WebXRHandJoint.RING_FINGER_PHALANX_PROXIMAL,        // XR_HAND_JOINT_RING_PROXIMAL_EXT
    WebXRHandJoint.RING_FINGER_PHALANX_INTERMEDIATE,    // XR_HAND_JOINT_RING_INTERMEDIATE_EXT
    WebXRHandJoint.RING_FINGER_PHALANX_DISTAL,          // XR_HAND_JOINT_RING_DISTAL_EXT
    WebXRHandJoint.RING_FINGER_TIP,                     // XR_HAND_JOINT_RING_TIP_EXT
    WebXRHandJoint.PINKY_FINGER_METACARPAL,             // XR_HAND_JOINT_LITTLE_METACARPAL_EXT
    WebXRHandJoint.PINKY_FINGER_PHALANX_PROXIMAL,       // XR_HAND_JOINT_LITTLE_PROXIMAL_EXT
    WebXRHandJoint.PINKY_FINGER_PHALANX_INTERMEDIATE,   // XR_HAND_JOINT_LITTLE_INTERMEDIATE_EXT
    WebXRHandJoint.PINKY_FINGER_PHALANX_DISTAL,         // XR_HAND_JOINT_LITTLE_DISTAL_EXT
    WebXRHandJoint.PINKY_FINGER_TIP,                    // XR_HAND_JOINT_LITTLE_TIP_EXT
];

export const START_T: number  = Date.now() / 1000;
const WEBXR_TO_OPENXR_ROTATION_CONVERSION = Quaternion.Identity(); // TODO: I don't think this is necessary as the rotation conventions SEEM to be the same. Remove when confirmed.
const convertedQuaternion = Quaternion.Identity();

function populateHandPoses(
    hand: WebXRHand | undefined,
    wristPose: ICarlOptionalTransform,
    jointPoses: ICarlOptionalTransform[]
): void {
    if (!hand) {
        wristPose.valid = false;
        for (let idx = 0; idx < jointPoses.length; ++idx) {
            jointPoses[idx].valid = false;
        }
        return;
    }
    const wristJoint = hand.getJointMesh(WebXRHandJoint.WRIST);
    wristPose.valid = true;
    wristPose.position.x = wristJoint.position.x;
    wristPose.position.y = wristJoint.position.y;
    wristPose.position.z = wristJoint.position.z;
    wristJoint.rotationQuaternion!.multiplyToRef(WEBXR_TO_OPENXR_ROTATION_CONVERSION, convertedQuaternion);
    wristPose.orientation.w = convertedQuaternion.w;
    wristPose.orientation.x = convertedQuaternion.x;
    wristPose.orientation.y = convertedQuaternion.y;
    wristPose.orientation.z = convertedQuaternion.z;
    for (let idx = 0; idx < jointPoses.length; ++idx) {
        const joint = hand.getJointMesh(OPENXR_JOINT_MAPPINGS[idx]);
        jointPoses[idx].valid = true;
        jointPoses[idx].position.x = joint.position.x;
        jointPoses[idx].position.y = joint.position.y;
        jointPoses[idx].position.z = joint.position.z;
        joint.rotationQuaternion!.multiplyToRef(WEBXR_TO_OPENXR_ROTATION_CONVERSION, convertedQuaternion);
        jointPoses[idx].orientation.w = convertedQuaternion.w;
        jointPoses[idx].orientation.x = convertedQuaternion.x;
        jointPoses[idx].orientation.y = convertedQuaternion.y;
        jointPoses[idx].orientation.z = convertedQuaternion.z;
    }
}

export function populateInputSample(hmd: WebXRCamera, leftHand: WebXRHand | undefined, rightHand: WebXRHand | undefined, sample: ICarlInputSample): void {
    sample.timestamp = (Date.now() / 1000) - START_T;

    sample.hmdPose.valid = true;
    sample.hmdPose.position.x = hmd.position.x;
    sample.hmdPose.position.y = hmd.position.y;
    sample.hmdPose.position.z = hmd.position.z;
    hmd.rotationQuaternion!.multiplyToRef(WEBXR_TO_OPENXR_ROTATION_CONVERSION, convertedQuaternion);
    sample.hmdPose.orientation.w = convertedQuaternion.w;
    sample.hmdPose.orientation.x = convertedQuaternion.x;
    sample.hmdPose.orientation.y = convertedQuaternion.y;
    sample.hmdPose.orientation.z = convertedQuaternion.z;

    populateHandPoses(leftHand, sample.leftWristPose, sample.leftHandJointPoses);
    populateHandPoses(rightHand, sample.rightWristPose, sample.rightHandJointPoses);
}

// Module-level scratch objects reused by applyJointSampleToMeshes to avoid per-frame allocation.
const _scratch = {
    vec: new Vector3(),
    quat: new Quaternion(),
    povMat: new Matrix(),
    sampleMat: new Matrix(),
    sampleToPovMat: new Matrix(),
    jointMat: new Matrix(),
};

export function applyJointSampleToMeshes(
    sample: ICarlInputSample,
    leftMeshes: AbstractMesh[],
    rightMeshes: AbstractMesh[],
    samplePosition: Vector3,
    sampleForward: Vector3,
    povPosition: Vector3,
    povForward: Vector3,
): void {
    povPosition.addToRef(povForward, _scratch.vec);
    Matrix.LookAtRHToRef(povPosition, _scratch.vec, Vector3.UpReadOnly, _scratch.povMat);
    samplePosition.addToRef(sampleForward, _scratch.vec);
    Matrix.LookAtRHToRef(samplePosition, _scratch.vec, Vector3.UpReadOnly, _scratch.sampleMat);

    _scratch.povMat.invertToRef(_scratch.povMat);
    _scratch.sampleMat.multiplyToRef(_scratch.povMat, _scratch.sampleToPovMat);

    for (let idx = 0; idx < OPENXR_JOINT_MAPPINGS.length; ++idx) {
        let pose = sample.leftHandJointPoses[idx];
        _scratch.quat.copyFromFloats(pose.orientation.x, pose.orientation.y, pose.orientation.z, pose.orientation.w);
        Matrix.FromQuaternionToRef(_scratch.quat, _scratch.jointMat);
        _scratch.jointMat.setTranslationFromFloats(pose.position.x, pose.position.y, pose.position.z);
        _scratch.jointMat.multiplyToRef(_scratch.sampleToPovMat, _scratch.jointMat);
        _scratch.jointMat.decompose(undefined, leftMeshes[idx].rotationQuaternion!, leftMeshes[idx].position);

        pose = sample.rightHandJointPoses[idx];
        _scratch.quat.copyFromFloats(pose.orientation.x, pose.orientation.y, pose.orientation.z, pose.orientation.w);
        Matrix.FromQuaternionToRef(_scratch.quat, _scratch.jointMat);
        _scratch.jointMat.setTranslationFromFloats(pose.position.x, pose.position.y, pose.position.z);
        _scratch.jointMat.multiplyToRef(_scratch.sampleToPovMat, _scratch.jointMat);
        _scratch.jointMat.decompose(undefined, rightMeshes[idx].rotationQuaternion!, rightMeshes[idx].position);
    }
}

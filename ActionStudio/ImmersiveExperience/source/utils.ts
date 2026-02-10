import { WebXRHandJoint, Quaternion, WebXRCamera, WebXRHand, AbstractMesh } from "@babylonjs/core";
import { ICarlInputSample } from "./carlInterfaces";

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

    let joint: AbstractMesh;
    if (!leftHand) {
        sample.leftWristPose.valid = false;
        for (let idx = 0; idx < sample.leftHandJointPoses.length; ++idx) {
            sample.leftHandJointPoses[idx].valid = false;
        }
    } else {
        joint = leftHand.getJointMesh(WebXRHandJoint.WRIST);
        sample.leftWristPose.valid = true;
        sample.leftWristPose.position.x = joint.position.x;
        sample.leftWristPose.position.y = joint.position.y;
        sample.leftWristPose.position.z = joint.position.z;
        joint.rotationQuaternion!.multiplyToRef(WEBXR_TO_OPENXR_ROTATION_CONVERSION, convertedQuaternion);
        sample.leftWristPose.orientation.w = convertedQuaternion.w;
        sample.leftWristPose.orientation.x = convertedQuaternion.x;
        sample.leftWristPose.orientation.y = convertedQuaternion.y;
        sample.leftWristPose.orientation.z = convertedQuaternion.z;

        for (let idx = 0; idx < sample.leftHandJointPoses.length; ++idx) {
            joint = leftHand.getJointMesh(OPENXR_JOINT_MAPPINGS[idx]);
            sample.leftHandJointPoses[idx].valid = true;
            sample.leftHandJointPoses[idx].position.x = joint.position.x;
            sample.leftHandJointPoses[idx].position.y = joint.position.y;
            sample.leftHandJointPoses[idx].position.z = joint.position.z;
            joint.rotationQuaternion!.multiplyToRef(WEBXR_TO_OPENXR_ROTATION_CONVERSION, convertedQuaternion);
            sample.leftHandJointPoses[idx].orientation.w = convertedQuaternion.w;
            sample.leftHandJointPoses[idx].orientation.x = convertedQuaternion.x;
            sample.leftHandJointPoses[idx].orientation.y = convertedQuaternion.y;
            sample.leftHandJointPoses[idx].orientation.z = convertedQuaternion.z;
        }
    }

    if (!rightHand) {
        sample.rightWristPose.valid = false;
        for (let idx = 0; idx < sample.rightHandJointPoses.length; ++idx) {
            sample.rightHandJointPoses[idx].valid = false;
        }
    } else {
        joint = rightHand.getJointMesh(WebXRHandJoint.WRIST);
        sample.rightWristPose.valid = true;
        sample.rightWristPose.position.x = joint.position.x;
        sample.rightWristPose.position.y = joint.position.y;
        sample.rightWristPose.position.z = joint.position.z;
        joint.rotationQuaternion!.multiplyToRef(WEBXR_TO_OPENXR_ROTATION_CONVERSION, convertedQuaternion);
        sample.rightWristPose.orientation.w = convertedQuaternion.w;
        sample.rightWristPose.orientation.x = convertedQuaternion.x;
        sample.rightWristPose.orientation.y = convertedQuaternion.y;
        sample.rightWristPose.orientation.z = convertedQuaternion.z;

        for (let idx = 0; idx < sample.rightHandJointPoses.length; ++idx) {
            joint = rightHand.getJointMesh(OPENXR_JOINT_MAPPINGS[idx]);
            sample.rightHandJointPoses[idx].valid = true;
            sample.rightHandJointPoses[idx].position.x = joint.position.x;
            sample.rightHandJointPoses[idx].position.y = joint.position.y;
            sample.rightHandJointPoses[idx].position.z = joint.position.z;
            joint.rotationQuaternion!.multiplyToRef(WEBXR_TO_OPENXR_ROTATION_CONVERSION, convertedQuaternion);
            sample.rightHandJointPoses[idx].orientation.w = convertedQuaternion.w;
            sample.rightHandJointPoses[idx].orientation.x = convertedQuaternion.x;
            sample.rightHandJointPoses[idx].orientation.y = convertedQuaternion.y;
            sample.rightHandJointPoses[idx].orientation.z = convertedQuaternion.z;
        }
    }
}

import { WebXRHandJoint, WebXRHand, WebXRSessionManager } from "@babylonjs/core";
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

const WEBXR_JOINT_NAMES = [
    "wrist",
    "thumb-metacarpal",
    "thumb-phalanx-proximal",
    "thumb-phalanx-distal",
    "thumb-tip",
    "index-finger-metacarpal",
    "index-finger-phalanx-proximal",
    "index-finger-phalanx-intermediate",
    "index-finger-phalanx-distal",
    "index-finger-tip",
    "middle-finger-metacarpal",
    "middle-finger-phalanx-proximal",
    "middle-finger-phalanx-intermediate",
    "middle-finger-phalanx-distal",
    "middle-finger-tip",
    "ring-finger-metacarpal",
    "ring-finger-phalanx-proximal",
    "ring-finger-phalanx-intermediate",
    "ring-finger-phalanx-distal",
    "ring-finger-tip",
    "pinky-finger-metacarpal",
    "pinky-finger-phalanx-proximal",
    "pinky-finger-phalanx-intermediate",
    "pinky-finger-phalanx-distal",
    "pinky-finger-tip",
];

const OPENXR_TO_WEBXR_JOINT_INDICES = [
    0,  // XR_HAND_JOINT_PALM_EXT -> wrist (approximation, WebXR has no palm)
    0,  // XR_HAND_JOINT_WRIST_EXT -> wrist
    1,  // XR_HAND_JOINT_THUMB_METACARPAL_EXT -> thumb-metacarpal
    2,  // XR_HAND_JOINT_THUMB_PROXIMAL_EXT -> thumb-phalanx-proximal
    3,  // XR_HAND_JOINT_THUMB_DISTAL_EXT -> thumb-phalanx-distal
    4,  // XR_HAND_JOINT_THUMB_TIP_EXT -> thumb-tip
    5,  // XR_HAND_JOINT_INDEX_METACARPAL_EXT -> index-finger-metacarpal
    6,  // XR_HAND_JOINT_INDEX_PROXIMAL_EXT -> index-finger-phalanx-proximal
    7,  // XR_HAND_JOINT_INDEX_INTERMEDIATE_EXT -> index-finger-phalanx-intermediate
    8,  // XR_HAND_JOINT_INDEX_DISTAL_EXT -> index-finger-phalanx-distal
    9,  // XR_HAND_JOINT_INDEX_TIP_EXT -> index-finger-tip
    10, // XR_HAND_JOINT_MIDDLE_METACARPAL_EXT -> middle-finger-metacarpal
    11, // XR_HAND_JOINT_MIDDLE_PROXIMAL_EXT -> middle-finger-phalanx-proximal
    12, // XR_HAND_JOINT_MIDDLE_INTERMEDIATE_EXT -> middle-finger-phalanx-intermediate
    13, // XR_HAND_JOINT_MIDDLE_DISTAL_EXT -> middle-finger-phalanx-distal
    14, // XR_HAND_JOINT_MIDDLE_TIP_EXT -> middle-finger-tip
    15, // XR_HAND_JOINT_RING_METACARPAL_EXT -> ring-finger-metacarpal
    16, // XR_HAND_JOINT_RING_PROXIMAL_EXT -> ring-finger-phalanx-proximal
    17, // XR_HAND_JOINT_RING_INTERMEDIATE_EXT -> ring-finger-phalanx-intermediate
    18, // XR_HAND_JOINT_RING_DISTAL_EXT -> ring-finger-phalanx-distal
    19, // XR_HAND_JOINT_RING_TIP_EXT -> ring-finger-tip
    20, // XR_HAND_JOINT_LITTLE_METACARPAL_EXT -> pinky-finger-metacarpal
    21, // XR_HAND_JOINT_LITTLE_PROXIMAL_EXT -> pinky-finger-phalanx-proximal
    22, // XR_HAND_JOINT_LITTLE_INTERMEDIATE_EXT -> pinky-finger-phalanx-intermediate
    23, // XR_HAND_JOINT_LITTLE_DISTAL_EXT -> pinky-finger-phalanx-distal
    24, // XR_HAND_JOINT_LITTLE_TIP_EXT -> pinky-finger-tip
];

export const START_T: number  = Date.now() / 1000;

function populatePoseFromXRPose(pose: XRPose, target: ICarlOptionalTransform): void {
    target.valid = true;
    target.position.x = pose.transform.position.x;
    target.position.y = pose.transform.position.y;
    target.position.z = pose.transform.position.z;
    target.orientation.w = pose.transform.orientation.w;
    target.orientation.x = pose.transform.orientation.x;
    target.orientation.y = -pose.transform.orientation.y;
    target.orientation.z = -pose.transform.orientation.z;
}

function invalidatePose(target: ICarlOptionalTransform): void {
    target.valid = false;
}

export function populateInputSample(
    frame: XRFrame,
    sessionManager: WebXRSessionManager,
    leftHand: WebXRHand | undefined,
    rightHand: WebXRHand | undefined,
    sample: ICarlInputSample
): void {
    sample.timestamp = (Date.now() / 1000) - START_T;

    const referenceSpace = sessionManager.referenceSpace;

    const viewerPose = frame.getViewerPose(referenceSpace);
    if (viewerPose) {
        sample.hmdPose.valid = true;
        sample.hmdPose.position.x = viewerPose.transform.position.x;
        sample.hmdPose.position.y = viewerPose.transform.position.y;
        sample.hmdPose.position.z = viewerPose.transform.position.z;
        sample.hmdPose.orientation.w = viewerPose.transform.orientation.w;
        sample.hmdPose.orientation.x = viewerPose.transform.orientation.x;
        sample.hmdPose.orientation.y = viewerPose.transform.orientation.y;
        sample.hmdPose.orientation.z = viewerPose.transform.orientation.z;
    } else {
        invalidatePose(sample.hmdPose);
    }

    if (!leftHand) {
        invalidatePose(sample.leftWristPose);
        for (let idx = 0; idx < sample.leftHandJointPoses.length; ++idx) {
            invalidatePose(sample.leftHandJointPoses[idx]);
        }
    } else {
        const xrHand = leftHand.xrController.inputSource.hand;
        if (xrHand) {
            const wristJointSpace = xrHand.get("wrist");
            if (wristJointSpace) {
                const wristPose = frame.getJointPose(wristJointSpace, referenceSpace);
                if (wristPose) {
                    populatePoseFromXRPose(wristPose, sample.leftWristPose);
                } else {
                    invalidatePose(sample.leftWristPose);
                }
            } else {
                invalidatePose(sample.leftWristPose);
            }

            for (let idx = 0; idx < sample.leftHandJointPoses.length; ++idx) {
                const webxrJointIdx = OPENXR_TO_WEBXR_JOINT_INDICES[idx];
                const jointName = WEBXR_JOINT_NAMES[webxrJointIdx] as XRHandJoint;
                const jointSpace = xrHand.get(jointName);
                if (jointSpace) {
                    const jointPose = frame.getJointPose(jointSpace, referenceSpace);
                    if (jointPose) {
                        populatePoseFromXRPose(jointPose, sample.leftHandJointPoses[idx]);
                    } else {
                        invalidatePose(sample.leftHandJointPoses[idx]);
                    }
                } else {
                    invalidatePose(sample.leftHandJointPoses[idx]);
                }
            }
        } else {
            invalidatePose(sample.leftWristPose);
            for (let idx = 0; idx < sample.leftHandJointPoses.length; ++idx) {
                invalidatePose(sample.leftHandJointPoses[idx]);
            }
        }
    }

    if (!rightHand) {
        invalidatePose(sample.rightWristPose);
        for (let idx = 0; idx < sample.rightHandJointPoses.length; ++idx) {
            invalidatePose(sample.rightHandJointPoses[idx]);
        }
    } else {
        const xrHand = rightHand.xrController.inputSource.hand;
        if (xrHand) {
            const wristJointSpace = xrHand.get("wrist");
            if (wristJointSpace) {
                const wristPose = frame.getJointPose(wristJointSpace, referenceSpace);
                if (wristPose) {
                    populatePoseFromXRPose(wristPose, sample.rightWristPose);
                } else {
                    invalidatePose(sample.rightWristPose);
                }
            } else {
                invalidatePose(sample.rightWristPose);
            }

            for (let idx = 0; idx < sample.rightHandJointPoses.length; ++idx) {
                const webxrJointIdx = OPENXR_TO_WEBXR_JOINT_INDICES[idx];
                const jointName = WEBXR_JOINT_NAMES[webxrJointIdx] as XRHandJoint;
                const jointSpace = xrHand.get(jointName);
                if (jointSpace) {
                    const jointPose = frame.getJointPose(jointSpace, referenceSpace);
                    if (jointPose) {
                        populatePoseFromXRPose(jointPose, sample.rightHandJointPoses[idx]);
                    } else {
                        invalidatePose(sample.rightHandJointPoses[idx]);
                    }
                } else {
                    invalidatePose(sample.rightHandJointPoses[idx]);
                }
            }
        } else {
            invalidatePose(sample.rightWristPose);
            for (let idx = 0; idx < sample.rightHandJointPoses.length; ++idx) {
                invalidatePose(sample.rightHandJointPoses[idx]);
            }
        }
    }
}

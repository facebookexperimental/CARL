import { AbstractMesh, WebXRCamera, WebXRHand, WebXRHandJoint } from "@babylonjs/core";

export interface ICarlPosition {
    x: number;
    y: number;
    z: number;
}

export interface ICarlOrientation {
    w: number;
    x: number;
    y: number;
    z: number;
}

export interface ICarlOptionalTransform {
    valid: boolean;
    position: ICarlPosition;
    orientation: ICarlOrientation;
}

export interface ICarlOptionalControllerState {
    valid: boolean;
    primaryClick: number;
    secondaryClick: number;
    thumbstickX: number;
    thumbstickY: number;
    thumbstickClick: number;
    squeezeValue: number;
    triggerValue: number;
}

export interface ICarlInputSample {
    timestamp: number;
    hmdPose: ICarlOptionalTransform;
    leftWristPose: ICarlOptionalTransform;
    rightWristPose: ICarlOptionalTransform;
    leftHandJointPoses: ICarlOptionalTransform[];
    rightHandJointPoses: ICarlOptionalTransform[];
    leftControllerState: ICarlOptionalControllerState;
    rightControllerState: ICarlOptionalControllerState;
}

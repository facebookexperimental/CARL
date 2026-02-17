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

export interface ICarlRecordingInspector {
    getStartTimestamp(): number;
    getEndTimestamp(): number;
    inspect(timestamp: number): ICarlInputSample;
    dispose(): void;
}

export interface ICarlExample {
    getStartTimestamp(): number;
    setStartTimestamp(timestamp: number): number;
    getEndTimestamp(): number;
    setEndTimestamp(timestamp: number): number;
    getRecordingInspector(): ICarlRecordingInspector;
    dispose(): void;
}

export interface ICarlDefinition {
    getDefaultSensitivity(): number;
    setDefaultSensitivity(sensitivity: number): void;
    download(): void;
    dispose(): void;
}

export interface ICarlRecognizer {
    currentScore(): number;
    getSensitivity(): number;
    setSensitivity(sensitivity: number): void;
    dispose(): void;
}

export interface ICarl {
    startRecording(): number;
    stopRecording(recordingId: number): ICarlExample;
    getActionTypes(): {name: string, typeId: number}[];
    createDefinition(actionTypeId: number, examples: ICarlExample[], counterexamples: ICarlExample[]): ICarlDefinition;
    createRecognizer(definition: ICarlDefinition): ICarlRecognizer;
    createInputSample(): ICarlInputSample;
    handleInputSample(sample: ICarlInputSample): void;
}

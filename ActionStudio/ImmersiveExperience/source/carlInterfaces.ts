/**
 * TypeScript interfaces mirroring the CARL WASM bindings exposed by NativeIntegration.
 *
 * These types are the contract between the Babylon.js ImmersiveExperience layer and the
 * Emscripten-generated C++ library.  The shapes must stay in sync with the native API;
 * if the native API changes, update these interfaces and the NativeIntegration wrapper.
 */
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
    stopRecording(recordingId: number, metadata?: object): ICarlExample;
    getActionTypes(): {name: string, typeId: number}[];
    draftDefinition(actionTypeId: number, examples: ICarlExample[], counterexamples: ICarlExample[]): ICarlDefinition;
    finalizeDefinition(definition: ICarlDefinition, metadata?: object): void;
    createRecognizer(definition: ICarlDefinition): ICarlRecognizer;
    createInputSample(): ICarlInputSample;
    handleInputSample(sample: ICarlInputSample): void;
}

export interface IInitialBlock<T> {
    value: T;
    color: string;
    id?: number;   // DB record id — used to persist trim edits
}

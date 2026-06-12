/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * Strongly-typed wrapper classes that bridge the CARL WASM native bindings to idiomatic
 * TypeScript usage.  These classes own the lifecycle of their underlying native (Emscripten)
 * objects and must be disposed when no longer needed to avoid WASM memory leaks.
 *
 * @module carlIntegration
 */

import {
    ICarl,
    ICarlDefinition,
    ICarlExample,
    ICarlInputSample,
    ICarlRecognizer,
    ICarlRecordingInspector,
} from "./carlInterfaces";
import { CarlInitOptions, NativeCarlModule, initializeCarl } from "./wasmLoader";

// The native Emscripten objects are untyped; alias `any` for readability.
// eslint-disable-next-line @typescript-eslint/no-explicit-any
type Native = any;

function nativeBytesToUint8Array(bytes: Native): Uint8Array {
    const jsBytes = new Uint8Array(bytes.size());
    for (let idx = 0; idx < jsBytes.length; ++idx) {
        jsBytes[idx] = bytes.get(idx);
    }
    bytes.delete();
    return jsBytes;
}

/**
 * Wraps a native CARL RecordingInspector, allowing frame-by-frame inspection
 * of a recording by timestamp.
 */
class CarlRecordingInspector implements ICarlRecordingInspector {
    private _nativeInspector: Native;

    constructor(nativeInspector: Native) {
        this._nativeInspector = nativeInspector;
    }

    getStartTimestamp(): number {
        return this._nativeInspector.startTimestamp();
    }

    getEndTimestamp(): number {
        return this._nativeInspector.endTimestamp();
    }

    inspect(timestamp: number): ICarlInputSample {
        return this._nativeInspector.inspect(timestamp);
    }

    dispose(): void {
        this._nativeInspector.delete();
    }
}

/**
 * Wraps a native CARL Example.  Holds the raw recording data plus the
 * user-selected trim window (startTimestamp / endTimestamp).
 *
 * Call `dispose()` when done to free the underlying WASM allocation.
 * Set `onDisposed` to be notified when disposal occurs.
 */
export class CarlExample implements ICarlExample {
    _nativeExample: Native;
    onDisposed: (() => void) | undefined = undefined;

    constructor(nativeExample: Native) {
        this._nativeExample = nativeExample;
    }

    getStartTimestamp(): number {
        return this._nativeExample.getStartTimestamp();
    }

    setStartTimestamp(timestamp: number): void {
        this._nativeExample.setStartTimestamp(timestamp);
    }

    getEndTimestamp(): number {
        return this._nativeExample.getEndTimestamp();
    }

    setEndTimestamp(timestamp: number): void {
        this._nativeExample.setEndTimestamp(timestamp);
    }

    getRecordingInspector(): ICarlRecordingInspector {
        return new CarlRecordingInspector(this._nativeExample.getRecordingInspector());
    }

    serialize(): Uint8Array {
        return nativeBytesToUint8Array(this._nativeExample.serialize());
    }

    dispose(): void {
        if (this.onDisposed) this.onDisposed();
        this._nativeExample.delete();
    }
}

/**
 * Wraps a native CARL Definition — a trained gesture/pose definition built
 * from a set of examples and counterexamples.
 *
 * Call `dispose()` when done to free the underlying WASM allocation.
 * Set `onDisposed` to be notified when disposal occurs.
 */
export class CarlDefinition implements ICarlDefinition {
    _nativeDefinition: Native;
    onDisposed: (() => void) | undefined = undefined;

    constructor(nativeDefinition: Native) {
        this._nativeDefinition = nativeDefinition;
    }

    serialize(): Uint8Array {
        return nativeBytesToUint8Array(this._nativeDefinition.serialize());
    }

    getActionType(): number {
        return this._nativeDefinition.getActionType();
    }

    getExamplesCount(): number {
        return this._nativeDefinition.getExamplesCount();
    }

    getCounterexamplesCount(): number {
        return this._nativeDefinition.getCounterexamplesCount();
    }

    getExample(idx: number): CarlExample {
        return new CarlExample(this._nativeDefinition.getExample(idx));
    }

    getCounterexample(idx: number): CarlExample {
        return new CarlExample(this._nativeDefinition.getCounterexample(idx));
    }

    getDefaultSensitivity(): number {
        return this._nativeDefinition.getDefaultSensitivity();
    }

    setDefaultSensitivity(sensitivity: number): void {
        this._nativeDefinition.setDefaultSensitivity(sensitivity);
    }

    dispose(): void {
        if (this.onDisposed) this.onDisposed();
        this._nativeDefinition.delete();
    }
}

/**
 * Wraps a native CARL Recognizer.  Used during live XR sessions to score
 * incoming input samples against a definition.
 *
 * Call `dispose()` when done to free the underlying WASM allocation.
 */
export class CarlRecognizer implements ICarlRecognizer {
    private _nativeRecognizer: Native;

    constructor(nativeRecognizer: Native) {
        this._nativeRecognizer = nativeRecognizer;
    }

    currentScore(): number {
        return this._nativeRecognizer.currentScore();
    }

    getSensitivity(): number {
        return this._nativeRecognizer.getSensitivity();
    }

    setSensitivity(sensitivity: number): void {
        this._nativeRecognizer.setSensitivity(sensitivity);
    }

    dispose(): void {
        this._nativeRecognizer.delete();
    }
}

/**
 * High-level integration facade over the CARL WASM module.
 *
 * Owns the CARL Session, manages in-progress recordings, maps action-type
 * IDs to human-readable names, and provides serialization/deserialization
 * helpers for examples and definitions.
 *
 * Instantiate via `CarlIntegration.CreateAsync()`.
 */
export class CarlIntegration implements ICarl {
    private _carl: Native;
    private _session: Native;
    private _inProgressRecordings = new Map<number, Native>();

    private _nextId = 1;

    private _actionTypes: { name: string; typeId: number }[] = [];
    private _idToActionType = new Map<number, Native>();

    onExampleCreated: ((example: CarlExample, metadata?: object) => void) | undefined = undefined;
    onDefinitionCreated: ((definition: CarlDefinition, metadata?: object) => void) | undefined = undefined;

    constructor(carl: Native) {
        this._carl = carl;
        this._session = new this._carl.Session();

        // TODO: Find a better way to deal with this block than this manual approach.
        let actionTypeId = 0;
        this._idToActionType.set(actionTypeId, this._carl.ACTION_TYPE.LeftHandPose);
        this._actionTypes.push({ name: "Left Hand Pose", typeId: actionTypeId });
        ++actionTypeId;
        this._idToActionType.set(actionTypeId, this._carl.ACTION_TYPE.LeftHandGesture);
        this._actionTypes.push({ name: "Left Hand Gesture", typeId: actionTypeId });
        ++actionTypeId;
        this._idToActionType.set(actionTypeId, this._carl.ACTION_TYPE.RightHandPose);
        this._actionTypes.push({ name: "Right Hand Pose", typeId: actionTypeId });
        ++actionTypeId;
        this._idToActionType.set(actionTypeId, this._carl.ACTION_TYPE.RightHandGesture);
        this._actionTypes.push({ name: "Right Hand Gesture", typeId: actionTypeId });
        ++actionTypeId;
        this._idToActionType.set(actionTypeId, this._carl.ACTION_TYPE.TwoHandGesture);
        this._actionTypes.push({ name: "Two Hand Gesture", typeId: actionTypeId });
        ++actionTypeId;
        this._idToActionType.set(actionTypeId, this._carl.ACTION_TYPE.LeftControllerGesture);
        this._actionTypes.push({ name: "Left Controller Gesture", typeId: actionTypeId });
        ++actionTypeId;
        this._idToActionType.set(actionTypeId, this._carl.ACTION_TYPE.RightControllerGesture);
        this._actionTypes.push({ name: "Right Controller Gesture", typeId: actionTypeId });
        ++actionTypeId;
        this._idToActionType.set(actionTypeId, this._carl.ACTION_TYPE.TwoControllerGesture);
        this._actionTypes.push({ name: "Two Controller Gesture", typeId: actionTypeId });
        ++actionTypeId;
        this._idToActionType.set(actionTypeId, this._carl.ACTION_TYPE.LeftWristTrajectory);
        this._actionTypes.push({ name: "Left Wrist Trajectory", typeId: actionTypeId });
        ++actionTypeId;
        this._idToActionType.set(actionTypeId, this._carl.ACTION_TYPE.RightWristTrajectory);
        this._actionTypes.push({ name: "Right Wrist Trajectory", typeId: actionTypeId });
        ++actionTypeId;
        this._idToActionType.set(actionTypeId, this._carl.ACTION_TYPE.LeftHandShape);
        this._actionTypes.push({ name: "Left Hand Shape", typeId: actionTypeId });
        ++actionTypeId;
        this._idToActionType.set(actionTypeId, this._carl.ACTION_TYPE.RightHandShape);
        this._actionTypes.push({ name: "Right Hand Shape", typeId: actionTypeId });
    }

    static async CreateAsync(options?: CarlInitOptions): Promise<CarlIntegration> {
        return new CarlIntegration(await initializeCarl(options));
    }

    startRecording(): number {
        const id = this._nextId;
        this._nextId += 1;
        this._inProgressRecordings.set(id, new this._carl.InProgressRecording());
        return id;
    }

    stopRecording(id: number, metadata?: object): CarlExample {
        const ipr = this._inProgressRecordings.get(id);
        this._inProgressRecordings.delete(id);

        const recording = new this._carl.Recording(ipr);
        const inspector = new this._carl.RecordingInspector(recording);
        const startT = inspector.startTimestamp();
        const endT = inspector.endTimestamp();

        const example = new this._carl.Example(recording, startT, endT);

        inspector.delete();
        recording.delete();
        ipr.delete();

        const jsExample = new CarlExample(example);
        if (this.onExampleCreated) this.onExampleCreated(jsExample, metadata);
        return jsExample;
    }

    getActionTypesMap(): { name: string; typeId: number }[] {
        return this._actionTypes;
    }

    draftDefinition(actionTypeId: number, examples: CarlExample[], counterexamples: CarlExample[]): CarlDefinition {
        const definition = new this._carl.Definition(this._idToActionType.get(actionTypeId));
        examples.forEach(example => {
            definition.addExample(example._nativeExample);
        });
        counterexamples.forEach(example => {
            definition.addCounterexample(example._nativeExample);
        });
        return new CarlDefinition(definition);
    }

    finalizeDefinition(definition: CarlDefinition, metadata?: object): void {
        if (this.onDefinitionCreated) this.onDefinitionCreated(definition, metadata);
    }

    createRecognizer(definition: CarlDefinition): CarlRecognizer {
        return new CarlRecognizer(new this._carl.Recognizer(this._session, definition._nativeDefinition));
    }

    createInputSample(): ICarlInputSample {
        return this._carl.createInputSample();
    }

    handleInputSample(sample: ICarlInputSample): void {
        this._session.addInput(sample);
        this._inProgressRecordings.forEach(ipr => {
            ipr.addInput(sample);
        });
    }

    tryDeserializeExample(jsBytes: Uint8Array): CarlExample | null {
        const nativeBytes = new this._carl.SerializedBytes();
        nativeBytes.resize(jsBytes.length);
        for (let idx = 0; idx < jsBytes.length; ++idx) {
            nativeBytes.set(idx, jsBytes[idx]);
        }
        const nativeExample = this._carl.Example.tryDeserialize(nativeBytes);
        nativeBytes.delete();
        return nativeExample ? new CarlExample(nativeExample) : null;
    }

    tryDeserializeDefinition(jsBytes: Uint8Array): CarlDefinition | null {
        const nativeBytes = new this._carl.SerializedBytes();
        nativeBytes.resize(jsBytes.length);
        for (let idx = 0; idx < jsBytes.length; ++idx) {
            nativeBytes.set(idx, jsBytes[idx]);
        }
        const nativeDefinition = this._carl.Definition.tryDeserialize(nativeBytes);
        nativeBytes.delete();
        return nativeDefinition ? new CarlDefinition(nativeDefinition) : null;
    }

    testDefinition(definition: CarlDefinition, testExamples: CarlExample[], sensitivity: number): { time: number; score: number }[][] {
        return testExamples.map(example => {
            const tempSession = new this._carl.Session();
            tempSession.tickCallbacks();
            const tempRecognizer = new this._carl.Recognizer(tempSession, definition._nativeDefinition);
            tempRecognizer.setSensitivity(sensitivity);

            const inspector = example.getRecordingInspector();
            const recStart = inspector.getStartTimestamp();
            const recEnd = inspector.getEndTimestamp();
            const STEPS = 60 * (recEnd - recStart);
            const dataPoints: { time: number; score: number }[] = [];

            for (let i = 0; i <= STEPS; i++) {
                const t = recStart + (i / STEPS) * (recEnd - recStart);
                const sample = inspector.inspect(t);
                tempSession.addInput(sample);
                dataPoints.push({ time: t - recStart, score: tempRecognizer.currentScore() });
            }

            inspector.dispose();
            tempRecognizer.delete();
            tempSession.delete();
            return dataPoints;
        });
    }
}

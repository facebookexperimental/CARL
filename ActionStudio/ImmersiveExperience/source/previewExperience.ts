/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * Standalone 2D preview experience rendered into a canvas element.
 *
 * Creates a minimal Babylon.js scene (no physics, no WebXR) that plays back a recorded
 * ICarlExample as animated joint spheres.  Used by the `/preview/:id` React route so
 * users can inspect, trim, and scrub through a gesture recording before saving it.
 */
import { ArcRotateCamera, Color3, Engine, HemisphericLight, MeshBuilder, Quaternion, Scene, StandardMaterial, Vector3 } from "@babylonjs/core";
import { ICarlExample, ICarlInputSample } from "./carlInterfaces";
import { applyJointSampleToMeshes } from "./utils";

// 26 entries (one per OPENXR_JOINT_MAPPINGS index), all 0.005 m (= 1 cm diameter).
const JOINT_RADII: number[] = new Array(26).fill(0.01);

export interface IPreviewExperienceHandle {
    setTime(relativeTime: number): void;
    setPlaying(playing: boolean): void;
    setTrimBounds(start: number, end: number): void;
    dispose(): void;
    readonly duration: number;
    readonly trimStart: number;
    readonly trimEnd: number;
    onTimeUpdate: ((t: number) => void) | null;
    onPlaybackEnd: (() => void) | null;
}

export async function initializePreviewExperienceAsync(canvas: HTMLCanvasElement, example: ICarlExample): Promise<IPreviewExperienceHandle> {
    const engine = new Engine(canvas, true);
    const scene = new Scene(engine);
    scene.useRightHandedSystem = true;

    const light = new HemisphericLight("light", new Vector3(0, 1, 0), scene);
    light.intensity = 0.7;

    const IMMITATION_ORIGIN = Vector3.RightHandedForwardReadOnly.scale(2);
    const camera = new ArcRotateCamera("camera", -Math.PI / 2, Math.PI / 3, 3, IMMITATION_ORIGIN.add(Vector3.UpReadOnly), scene);
    camera.upperRadiusLimit = 3;
    camera.lowerRadiusLimit = 1;
    camera.wheelDeltaPercentage = 0.03;
    camera.maxZ = 10;
    camera.minZ = 0.05;
    camera.lowerBetaLimit = Math.PI / 4;
    camera.upperBetaLimit = 3 * Math.PI / 4;
    camera.attachControl(canvas, true);

    const leftMat = new StandardMaterial("leftMat", scene);
    leftMat.diffuseColor = Color3.Blue();
    const rightMat = new StandardMaterial("rightMat", scene);
    rightMat.diffuseColor = Color3.Red();

    const leftMeshes = JOINT_RADII.map((r, i) => {
        const mesh = MeshBuilder.CreateSphere(`left_joint_${i}`, { diameter: r * 2 }, scene);
        mesh.material = leftMat;
        mesh.rotationQuaternion = Quaternion.Identity();
        return mesh;
    });
    const rightMeshes = JOINT_RADII.map((r, i) => {
        const mesh = MeshBuilder.CreateSphere(`right_joint_${i}`, { diameter: r * 2 }, scene);
        mesh.material = rightMat;
        mesh.rotationQuaternion = Quaternion.Identity();
        return mesh;
    });

    const inspector = example.getRecordingInspector();
    const recordingStart = inspector.getStartTimestamp();
    const duration = inspector.getEndTimestamp() - recordingStart;
    const trimStart = example.getStartTimestamp() - recordingStart;
    const trimEnd = example.getEndTimestamp() - recordingStart;

    // Spatial transform constants matching ExamplePreviewer / InputPuppet
    const samplePosition = Vector3.ZeroReadOnly;
    const sampleForward = Vector3.RightHandedForwardReadOnly;
    const povPosition = IMMITATION_ORIGIN;
    const povForward = Vector3.RightHandedBackwardReadOnly;

    function positionJointMeshes(sample: ICarlInputSample): void {
        applyJointSampleToMeshes(sample, leftMeshes, rightMeshes,
            samplePosition, sampleForward, povPosition, povForward);
    }

    let currentT = trimStart;
    let isPlaying = false;
    let activeTrimStart = trimStart;
    let activeTrimEnd = trimEnd;
    let onTimeUpdate: ((t: number) => void) | null = null;
    let onPlaybackEnd: (() => void) | null = null;

    // Render initial pose at trimStart
    positionJointMeshes(inspector.inspect(currentT + recordingStart));

    scene.onBeforeRenderObservable.add(() => {
        if (!isPlaying) return;
        currentT += scene.deltaTime / 1000;
        if (currentT >= activeTrimEnd) {
            currentT = activeTrimEnd;
            isPlaying = false;
            positionJointMeshes(inspector.inspect(currentT + recordingStart));
            onTimeUpdate?.(currentT);
            onPlaybackEnd?.();
            return;
        }
        positionJointMeshes(inspector.inspect(currentT + recordingStart));
        onTimeUpdate?.(currentT);
    });

    engine.runRenderLoop(() => scene.render());

    const handle: IPreviewExperienceHandle = {
        get duration() { return duration; },
        get trimStart() { return trimStart; },
        get trimEnd() { return trimEnd; },
        get onTimeUpdate() { return onTimeUpdate; },
        set onTimeUpdate(v) { onTimeUpdate = v; },
        get onPlaybackEnd() { return onPlaybackEnd; },
        set onPlaybackEnd(v) { onPlaybackEnd = v; },
        setTime(relativeTime: number) {
            currentT = relativeTime;
            positionJointMeshes(inspector.inspect(relativeTime + recordingStart));
        },
        setPlaying(playing: boolean) {
            if (playing) {
                currentT = activeTrimStart;
            }
            isPlaying = playing;
        },
        setTrimBounds(start: number, end: number) {
            activeTrimStart = start;
            activeTrimEnd = end;
        },
        dispose() {
            engine.stopRenderLoop();
            inspector.dispose();
            scene.dispose();
            engine.dispose();
        },
    };

    return handle;
}

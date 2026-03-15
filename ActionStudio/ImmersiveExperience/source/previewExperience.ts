import { ArcRotateCamera, Color3, Engine, HemisphericLight, Matrix, MeshBuilder, Quaternion, Scene, StandardMaterial, Vector3 } from "@babylonjs/core";
import { ICarlExample, ICarlInputSample } from "./carlInterfaces";

// 26 entries (one per OPENXR_JOINT_MAPPINGS index), all 0.005 m (= 1 cm diameter).
const JOINT_RADII: number[] = new Array(26).fill(0.005);

export interface IPreviewExperienceHandle {
    setTime(relativeTime: number): void;
    setPlaying(playing: boolean): void;
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

    const light = new HemisphericLight("light", new Vector3(0, 1, 0), scene);
    light.intensity = 0.7;

    const IMMITATION_ORIGIN = Vector3.RightHandedForwardReadOnly.scale(2);
    const camera = new ArcRotateCamera("camera", -Math.PI / 2, Math.PI / 3, 3, IMMITATION_ORIGIN, scene);
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

    const scratchVec = new Vector3();
    const scratchQuat = new Quaternion();
    const povMat = new Matrix();
    const sampleMat = new Matrix();
    const sampleToPovMat = new Matrix();
    const jointMat = new Matrix();

    function positionJointMeshes(sample: ICarlInputSample): void {
        povPosition.addToRef(povForward, scratchVec);
        Matrix.LookAtRHToRef(povPosition, scratchVec, Vector3.UpReadOnly, povMat);
        samplePosition.addToRef(sampleForward, scratchVec);
        Matrix.LookAtRHToRef(samplePosition, scratchVec, Vector3.UpReadOnly, sampleMat);

        povMat.invertToRef(povMat);
        sampleMat.multiplyToRef(povMat, sampleToPovMat);

        for (let idx = 0; idx < 26; ++idx) {
            let pose = sample.leftHandJointPoses[idx];
            scratchQuat.copyFromFloats(pose.orientation.x, pose.orientation.y, pose.orientation.z, pose.orientation.w);
            Matrix.FromQuaternionToRef(scratchQuat, jointMat);
            jointMat.setTranslationFromFloats(pose.position.x, pose.position.y, pose.position.z);
            jointMat.multiplyToRef(sampleToPovMat, jointMat);
            jointMat.decompose(undefined, leftMeshes[idx].rotationQuaternion!, leftMeshes[idx].position);

            pose = sample.rightHandJointPoses[idx];
            scratchQuat.copyFromFloats(pose.orientation.x, pose.orientation.y, pose.orientation.z, pose.orientation.w);
            Matrix.FromQuaternionToRef(scratchQuat, jointMat);
            jointMat.setTranslationFromFloats(pose.position.x, pose.position.y, pose.position.z);
            jointMat.multiplyToRef(sampleToPovMat, jointMat);
            jointMat.decompose(undefined, rightMeshes[idx].rotationQuaternion!, rightMeshes[idx].position);
        }
    }

    let currentT = trimStart;
    let isPlaying = false;
    let onTimeUpdate: ((t: number) => void) | null = null;
    let onPlaybackEnd: (() => void) | null = null;

    // Render initial pose at trimStart
    positionJointMeshes(inspector.inspect(currentT + recordingStart));

    scene.onBeforeRenderObservable.add(() => {
        if (!isPlaying) return;
        currentT += scene.deltaTime / 1000;
        if (currentT >= trimEnd) {
            currentT = trimEnd;
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
            if (playing && currentT >= trimEnd) {
                currentT = trimStart;
            }
            isPlaying = playing;
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

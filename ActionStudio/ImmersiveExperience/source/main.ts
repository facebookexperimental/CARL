import { _InstancesBatch, AbstractMesh, Color3, Engine, FreeCamera, HavokPlugin, HemisphericLight, IDisposable, Matrix, Mesh, MeshBuilder, Observable, Observer, PBRMaterial, PhysicsAggregate, PhysicsShape, PhysicsShapeType, PhysicsVortexEventData, Quaternion, ReciprocalBlock, Scene, Tools, Vector3, WebXRCamera, WebXRControllerPointerSelection, WebXRFeatureName, WebXRHand, WebXRHandJoint, WebXRState } from "@babylonjs/core";
import HavokPhysics from "@babylonjs/havok";
import "@babylonjs/loaders/glTF/2.0";
import { ICarl, ICarlExample, ICarlInputSample, ICarlRecordingInspector } from "./carlInterfaces";
import { HandPinchGrabber } from "./handPinchGrabber";
import { PhysicsEnabledScene } from "./physicsEnabledScene";
import { PhysicsGrabBehavior } from "./physicsGrabBehavior";

const OPENXR_JOINT_MAPPINGS = [
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

const START_T: number  = Date.now() / 1000;
const WEBXR_TO_OPENXR_ROTATION_CONVERSION = Quaternion.Identity(); // TODO: I don't think this is necessary as the rotation conventions SEEM to be the same. Remove when confirmed.
const convertedQuaternion = Quaternion.Identity();
function populateInputSample(hmd: WebXRCamera, leftHand: WebXRHand | undefined, rightHand: WebXRHand | undefined, sample: ICarlInputSample): void {
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

export class InputPuppet {
    private _leftHandMeshes: AbstractMesh[];
    private _rightHandMeshes: AbstractMesh[];
    private _scratchVec: Vector3 = new Vector3();
    private _scratchQuat: Quaternion = new Quaternion();
    private _povMat: Matrix = new Matrix();
    private _sampleMat: Matrix = new Matrix();
    private _sampleToPovMat: Matrix = new Matrix();
    private _jointMat: Matrix = new Matrix();

    public constructor(leftHand: WebXRHand, rightHand: WebXRHand) {
        this._leftHandMeshes = OPENXR_JOINT_MAPPINGS.map(joint => {
            const mesh = leftHand.getJointMesh(joint);
            return mesh.clone(joint, null)!;
        });
        this._rightHandMeshes = OPENXR_JOINT_MAPPINGS.map(joint => {
            const mesh = rightHand.getJointMesh(joint);
            return mesh.clone(joint, null)!;
        });

        this.setEnabled(false);
    }

    public immitateInputSample(sample: ICarlInputSample, samplePosition: Vector3, sampleForward: Vector3, povPosition: Vector3, povForward: Vector3) {
        povPosition.addToRef(povForward, this._scratchVec);
        Matrix.LookAtRHToRef(povPosition, this._scratchVec, Vector3.UpReadOnly, this._povMat);
        samplePosition.addToRef(sampleForward, this._scratchVec);
        Matrix.LookAtRHToRef(samplePosition, this._scratchVec, Vector3.UpReadOnly, this._sampleMat);

        this._povMat.invertToRef(this._povMat);
        this._sampleMat.multiplyToRef(this._povMat, this._sampleToPovMat);

        for (let idx = 0; idx < OPENXR_JOINT_MAPPINGS.length; ++idx) {
            let pose = sample.leftHandJointPoses[idx];
            this._scratchQuat.copyFromFloats(pose.orientation.x, pose.orientation.y, pose.orientation.z, pose.orientation.w);
            Matrix.FromQuaternionToRef(this._scratchQuat, this._jointMat);
            this._jointMat.setTranslationFromFloats(pose.position.x, pose.position.y, pose.position.z);
            this._jointMat.multiplyToRef(this._sampleToPovMat, this._jointMat);
            this._jointMat.decompose(undefined, this._leftHandMeshes[idx].rotationQuaternion!, this._leftHandMeshes[idx].position);
            
            pose = sample.rightHandJointPoses[idx];
            this._scratchQuat.copyFromFloats(pose.orientation.x, pose.orientation.y, pose.orientation.z, pose.orientation.w);
            Matrix.FromQuaternionToRef(this._scratchQuat, this._jointMat);
            this._jointMat.setTranslationFromFloats(pose.position.x, pose.position.y, pose.position.z);
            this._jointMat.multiplyToRef(this._sampleToPovMat, this._jointMat);
            this._jointMat.decompose(undefined, this._rightHandMeshes[idx].rotationQuaternion!, this._rightHandMeshes[idx].position);
        }
    }

    public setEnabled(enabled: boolean): void {
        this._leftHandMeshes.map(mesh => mesh.setEnabled(enabled));
        this._rightHandMeshes.map(mesh => mesh.setEnabled(enabled));
    }
}

class ExampleBlockSpawner {
    private _scene: PhysicsEnabledScene;
    private _nextId: number;
    private _exampleBlock: AbstractMesh;

    constructor (scene: PhysicsEnabledScene) {
        this._scene = scene;
        this._nextId = 0;
        this._exampleBlock = scene.getMeshByName("example_block")!;
        this._exampleBlock.setEnabled(false);
    }

    spawnNewExampleBlock(): AbstractMesh {
            const clone = this._exampleBlock.clone(`example_block_${++this._nextId}`, null)!;
            const aggregate = new PhysicsAggregate(clone, PhysicsShapeType.BOX, { mass: 1 }, this._scene);
            aggregate.transformNode.rotationQuaternion = Quaternion.FromEulerAngles(
                aggregate.transformNode.rotation.x, 
                aggregate.transformNode.rotation.y, 
                aggregate.transformNode.rotation.z);
            this._scene.grabbables.push(PhysicsGrabBehavior.attach(clone));
            clone.setEnabled(true);
            return clone;
    }
}

enum PokeButtonState {
    Default,
    Poked,
    Cooldown,
}

class PokeButton implements IDisposable {
    private scene: PhysicsEnabledScene;
    private onButton: AbstractMesh;
    private offButton: AbstractMesh | undefined;

    private beforeRenderObserver: Observer<Scene>;
    private state: PokeButtonState = PokeButtonState.Default;
    private cooldown: number = 0;

    private static readonly POKE_ENTER_THRESHOLD = 0.015;
    private static readonly POKE_EXIT_THRESHOLD = 0.02;

    public onPokeObservable: Observable<boolean> = new Observable<boolean>();

    public constructor(scene: PhysicsEnabledScene, onButtonName: string, offButtonName: string | undefined) {
        this.scene = scene;
        this.onButton = this.scene.getMeshByName(onButtonName)!;
        if (offButtonName) {
            this.offButton = this.scene.getMeshByName(offButtonName)!;
        }
        this.offButton?.setEnabled(false);
        this.onPokeObservable.add(isPoked => {
            if (isPoked && this.offButton) {
                this.onButton.setEnabled(!this.onButton.isEnabled());
                this.offButton.setEnabled(!this.offButton.isEnabled());
            }
        });

        this.beforeRenderObserver = this.scene.onBeforeRenderObservable.add(() => this.update());
        this.scene.onDisposeObservable.addOnce(() => this.dispose());
    }

    private update(): void {
        this.cooldown -= this.scene.deltaTime;

        let minDist: number = 10;
        if (this.scene.leftHand) {
            const tipPosition = this.scene.leftHand.getJointMesh(WebXRHandJoint.INDEX_FINGER_TIP).position;
            minDist = Math.min(Vector3.Distance(tipPosition, this.onButton.position), minDist);
        }
        if (this.scene.rightHand) {
            const tipPosition = this.scene.rightHand.getJointMesh(WebXRHandJoint.INDEX_FINGER_TIP).position;
            minDist = Math.min(Vector3.Distance(tipPosition, this.onButton.position), minDist);
        }

        switch (this.state) {
            case PokeButtonState.Default:
                if (minDist < PokeButton.POKE_ENTER_THRESHOLD) {
                    this.state = PokeButtonState.Poked;
                    this.onPokeObservable.notifyObservers(true);
                }
                break;
            case PokeButtonState.Poked:
                if (minDist > PokeButton.POKE_EXIT_THRESHOLD) {
                    this.cooldown = 500;
                    this.state = PokeButtonState.Cooldown;
                }
                break;
            case PokeButtonState.Cooldown:
                if (minDist < PokeButton.POKE_EXIT_THRESHOLD) {
                    this.state = PokeButtonState.Poked;
                } else if (this.cooldown < 0) {
                    this.state = PokeButtonState.Default;
                    this.onPokeObservable.notifyObservers(false);
                }
                break;
        }
    }

    public dispose(): void {
        this.onPokeObservable.clear();
        this.beforeRenderObserver.remove();
    }
}

class ExamplePreviewer {
    private scene: PhysicsEnabledScene;
    private previewing: boolean = false;
    private shouldPlay: boolean = false;

    public constructor(scene: PhysicsEnabledScene) {
        this.scene = scene;
    }

    public previewExample(example: ICarlExample): () => void {
        this.previewing = true;
        this.scene.onBeforeRenderObservable.runCoroutineAsync(this.previewExampleCoroutine(example));
        return () => {
            this.previewing = false;
        }
    }

    public play(): void {
        this.shouldPlay = true;
    }

    private *previewExampleCoroutine(example: ICarlExample) {
        const inspector = example.getRecordingInspector();
        const duration = inspector.getEndTimestamp() - inspector.getStartTimestamp();
        this.scene.inputPuppet?.immitateInputSample(inspector.inspect(inspector.getStartTimestamp()), Vector3.ZeroReadOnly, Vector3.RightHandedBackwardReadOnly, Vector3.RightHandedBackwardReadOnly, Vector3.RightHandedForwardReadOnly);
        this.scene.inputPuppet?.setEnabled(true);

        let currentT: number | undefined = undefined;

        let minT = example.getStartTimestamp();
        let maxT = example.getEndTimestamp();

        const sliderObservers: Observer<void>[] = [];
        for (let idx = 0; idx < this.scene.sliders.length; ++idx) {
            const slider = this.scene.sliders[idx];
            if (idx == 0) {
                slider.value = (minT - inspector.getStartTimestamp()) / duration;
            } else {
                slider.value = (maxT - inspector.getStartTimestamp()) / duration;
            }

            const otherSlider = this.scene.sliders[(idx + 1) % 2];
            sliderObservers.push(slider.onUpdatedObservable.add(() => {
                currentT = undefined;

                let t = slider.value;
                t = t * duration + inspector.getStartTimestamp();

                if (slider.value < otherSlider.value) {
                    minT = t;
                } else {
                    maxT = t;
                }

                const sample = inspector.inspect(t);
                this.scene.inputPuppet?.immitateInputSample(sample, Vector3.ZeroReadOnly, Vector3.RightHandedBackwardReadOnly, Vector3.RightHandedBackwardReadOnly, Vector3.RightHandedForwardReadOnly);
            }));
        }
        
        while (this.previewing) {
            if (this.shouldPlay) {
                currentT = minT;
                this.shouldPlay = false;
            }

            if (currentT) {
                const sample = inspector.inspect(currentT);
                this.scene.inputPuppet?.immitateInputSample(sample, Vector3.ZeroReadOnly, Vector3.RightHandedBackwardReadOnly, Vector3.RightHandedBackwardReadOnly, Vector3.RightHandedForwardReadOnly);

                currentT += (this.scene.deltaTime / 1000);
                if (currentT > maxT) {
                    currentT = undefined;
                }
            }
            yield;
        }

        sliderObservers.forEach(observer => {
            observer.remove();
        });
        example.setStartTimestamp(minT);
        example.setEndTimestamp(maxT);
        inspector.dispose();
        this.scene.inputPuppet?.setEnabled(false);
    }
}

export async function initializeImmersiveExperienceAsync(canvas: HTMLCanvasElement, carl: ICarl): Promise<void> {
    const engine = new Engine(canvas, true);

    const resizeHandler = () => {
        engine.resize();
    };
    window.addEventListener("resize", resizeHandler);
    engine.onDisposeObservable.add(() => {
        window.removeEventListener("resize", resizeHandler);
    });

    const havok = await HavokPhysics();
    const hk = new HavokPlugin(true, havok);

    const scene = await PhysicsEnabledScene.loadAsync("./assets/action_studio.glb", engine, hk);
    scene.useRightHandedSystem = true;
    engine.runRenderLoop(() => {
        scene.render();
    });

    scene.enablePhysics(new Vector3(0, -9.8, 0), hk);

    const camera = new FreeCamera("camera1", new Vector3(0, 1, 0), scene);
    camera.attachControl(canvas, true);

    const light = new HemisphericLight("light", new Vector3(0, 1, 0), scene);
    light.intensity = 0.7;
    
    const spawner = new ExampleBlockSpawner(scene);

    const xr = await scene.createDefaultXRExperienceAsync();

    xr.baseExperience.onStateChangedObservable.add((state) => {
        if (state === WebXRState.IN_XR) {
            canvas.style.display = "none";
        } else {
            canvas.style.display = "block";
        }
    });

    xr.baseExperience.featuresManager.disableFeature(WebXRControllerPointerSelection.Name);

    let meshJoint = MeshBuilder.CreateBox("jointParent", { size: 1 });
    meshJoint.isVisible = false;

    const handTracking = xr.baseExperience.featuresManager.enableFeature(WebXRFeatureName.HAND_TRACKING, "latest", {
        xrInput: xr.input,
        jointMeshes: {
            sourceMesh: meshJoint,
            keepOriginalVisible: true,
        },
        handMeshes: {
            disableDefaultMeshes: true
        }
    });

    handTracking.onHandAddedObservable.add((hand) => {
        switch (hand.xrController.inputSource.handedness) {
            case "left":
                scene.leftHand = hand;
                break;
            case "right":
                scene.rightHand = hand;
                break;
        }

        const grabber = new HandPinchGrabber("handGrabber", hand, xr.input.xrSessionManager, scene);
        grabber.onGrabChanged.add(() => {
            PhysicsGrabBehavior.handleGrab(grabber);
        });
        scene.grabbers.set(hand.xrController.uniqueId, grabber);
        hand.xrController.inputSource.handedness
    });

    handTracking.onHandRemovedObservable.add((hand) => {
        if (scene.grabbers.has(hand.xrController.uniqueId)) {
            scene.grabbers.get(hand.xrController.uniqueId)?.dispose();
            scene.grabbers.delete(hand.xrController.uniqueId);
        }

        switch (hand.xrController.inputSource.handedness) {
            case "left":
                if (scene.leftHand?.xrController.uniqueId === hand.xrController.uniqueId) {
                    scene.leftHand = undefined;
                }
                break;
            case "right":
                if (scene.rightHand?.xrController.uniqueId === hand.xrController.uniqueId) {
                    scene.rightHand = undefined;
                }
        }
    });

    const cachedInputSample = carl.createInputSample();
    function createInputSample(): ICarlInputSample {
        return {...cachedInputSample};
    }

    const previewer = new ExamplePreviewer(scene);

    // TODO: Testing the recorder.
    const recordingStartPoke = new PokeButton(scene, "recording_start", "recording_stop");
    let recording: number | undefined = undefined;
    let example: ICarlExample | undefined = undefined;
    recordingStartPoke.onPokeObservable.add(poked => {
        if (poked) {
            if (recording) {
                example = carl.stopRecording(recording);
                recording = undefined;
                spawner.spawnNewExampleBlock();

                // TODO: Testing the previewer.
                Tools.DelayAsync(2000).then(() => {
                    const stopper = previewer.previewExample(example!);
                    Tools.DelayAsync(2000).then(() => {
                        previewer.play();
                    });
                    Tools.DelayAsync(10000).then(() => {
                        stopper();
                    });
                });
            } else {
                recording = carl.startRecording();
            }
        }
    });

    xr.input.xrSessionManager.onXRFrameObservable.add((frame) => {
        let sample = createInputSample();
        populateInputSample(xr.input.xrCamera, scene.leftHand, scene.rightHand, sample);
        carl.handleInputSample(sample);

        if (scene.inputPuppet === undefined && scene.leftHand && scene.rightHand && scene.leftHand.getJointMesh(OPENXR_JOINT_MAPPINGS[0]).scaling.x < 1) {
            scene.inputPuppet = new InputPuppet(scene.leftHand, scene.rightHand);
        }
    });

    // Testing CARL functionality sans interaction (yet)
    xr.input.xrSessionManager.onXRFrameObservable.addOnce(async () => {
        await Tools.DelayAsync(2000);
        const recId = carl.startRecording();
        const recId2 = carl.startRecording();
        await Tools.DelayAsync(500);
        const example = carl.stopRecording(recId);
        const definition = carl.createDefinition(2, [example], []);
        const recognizer = carl.createRecognizer(definition);
        
        example.dispose();
        definition.dispose();
    });

    

        // if (scene.inputPuppet) {
        //     scene.inputPuppet.immitateInputSample(sample, Vector3.ZeroReadOnly, Vector3.RightHandedBackwardReadOnly, Vector3.RightHandedBackwardReadOnly.scale(0.1), Vector3.RightHandedBackwardReadOnly);
        // }
}

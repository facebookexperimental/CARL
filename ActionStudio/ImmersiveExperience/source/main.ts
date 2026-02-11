import { _InstancesBatch, Color4, Engine, FreeCamera, HavokPlugin, HemisphericLight, MeshBuilder, TransformNode as AbstractMesh, Vector3, WebXRControllerPointerSelection, WebXRFeatureName, WebXRState } from "@babylonjs/core";
import HavokPhysics from "@babylonjs/havok";
import "@babylonjs/loaders/glTF/2.0";
import { ICarl, ICarlExample, ICarlInputSample } from "./carlInterfaces";
import { HandPinchGrabber } from "./handPinchGrabber";
import { PhysicsEnabledScene } from "./physicsEnabledScene";
import { PhysicsGrabBehavior } from "./physicsGrabBehavior";
import { OPENXR_JOINT_MAPPINGS, populateInputSample } from "./utils";
import { ExampleBlockSpawner } from "./exampleBlockSpawner";
import { PokeButton } from "./pokeButton";
import { ExamplePreviewer } from "./examplePreviewer";
import { InputPuppet } from "./inputPuppet";
import { RecognitionGraph } from "./recognitionGraph";

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

    scene.transformNodes.map(node => {
        if (node.name.startsWith("invisible")) {
            node.parent!.setEnabled(false);
        }
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
    let previewStopper: any = null;
    let currentlyPreviewed: AbstractMesh | null = null;

    const editorMesh = scene.getMeshByName("editor");

    const addPreviewingSupportToExampleBlock = (block: AbstractMesh, blockGrabbable: PhysicsGrabBehavior) => {
        const tryStopPreviewing = () => {
            if (previewStopper) previewStopper();
            previewStopper = null;
            currentlyPreviewed = null;
        };
        blockGrabbable.onGrabStartedObservable.add(() =>  {
            if (currentlyPreviewed === block) {
                tryStopPreviewing();
            }
        });
        blockGrabbable.onGrabEndedObservable.add(() => {
            if (Vector3.Distance(block.position, editorMesh!.position) < 0.08) {
                tryStopPreviewing();
                currentlyPreviewed = block;
                previewStopper = previewer.previewExample(spawner.getExampleFromBlock(block)!);
            }
        });
    };

    const currentGraph = new RecognitionGraph(scene, Color4.FromInts(255, 255, 255, 255));
    scene.onDisposeObservable.add(() => {
        currentGraph.dispose();
    });
    const isExampleMatrix = scene.getMeshByName("examples_volume")!.computeWorldMatrix(true).invert();
    const scratchVec = new Vector3();
    const addCurrentDefinitionSupportToExampleBlock = (block: AbstractMesh, blockGrabbable: PhysicsGrabBehavior) => {
        blockGrabbable.onGrabStartedObservable.add(() =>  {
            // TODO: Do this for real.
            currentGraph.recognizer = undefined;
        });
        blockGrabbable.onGrabEndedObservable.add(() => {
            Vector3.TransformCoordinatesToRef(block.position, isExampleMatrix, scratchVec);
            if (Math.abs(scratchVec.x) < 1 && Math.abs(scratchVec.y) < 1 && Math.abs(scratchVec.z) < 1) {
                const example = spawner.getExampleFromBlock(block)!;
                const definition = carl.createDefinition(2, [example], []);
                currentGraph.recognizer = carl.createRecognizer(definition);
                definition.dispose();
            }
        });
    };

    // TODO: Testing the recorder.
    const recordingPokeButton = new PokeButton(scene, "recording_start", "recording_stop");
    let recording: number | undefined = undefined;
    let example: ICarlExample | undefined = undefined;
    recordingPokeButton.onPokeObservable.add(poked => {
        if (poked) {
            if (recording) {
                example = carl.stopRecording(recording);
                recording = undefined;
                
                const block = spawner.spawnNewExampleBlock(example);
                const blockGrabbable = PhysicsGrabBehavior.get(block);
                addPreviewingSupportToExampleBlock(block, blockGrabbable);
                addCurrentDefinitionSupportToExampleBlock(block, blockGrabbable);
            } else {
                recording = carl.startRecording();
                console.log("Recording! " + recording);
            }
        }
    });

    const playPoke = new PokeButton(scene, "example_play");
    playPoke.onPokeObservable.add(poked => {
        if (poked && currentlyPreviewed) {
            previewer.play();
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
}


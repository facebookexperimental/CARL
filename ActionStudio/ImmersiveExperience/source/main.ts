import { _InstancesBatch, Color4, Engine, FreeCamera, HavokPlugin, HemisphericLight, MeshBuilder, TransformNode as AbstractMesh, Vector3, WebXRControllerPointerSelection, WebXRFeatureName, WebXRState, TransformNode } from "@babylonjs/core";
import HavokPhysics from "@babylonjs/havok";
import "@babylonjs/loaders/glTF/2.0";
import { ICarl, ICarlDefinition, ICarlExample, ICarlInputSample } from "./carlInterfaces";
import { HandPinchGrabber } from "./handPinchGrabber";
import { PhysicsEnabledScene } from "./physicsEnabledScene";
import { PhysicsGrabBehavior } from "./physicsGrabBehavior";
import { OPENXR_JOINT_MAPPINGS, populateInputSample } from "./utils";
import { ExampleBlockSpawner } from "./exampleBlockSpawner";
import { PokeButton } from "./pokeButton";
import { ExamplePreviewer } from "./examplePreviewer";
import { InputPuppet } from "./inputPuppet";
import { RecognitionGraph } from "./recognitionGraph";
import { SliderBehavior } from "./slider";

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

    const xr = await scene.createDefaultXRExperienceAsync({
        uiOptions: {
            sessionMode: "immersive-ar",
        }
    });

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

    const exampleBlocks = new Set<TransformNode>();
    const counterexampleBlocks = new Set<TransformNode>();

    const currentGraph = new RecognitionGraph(scene, Color4.FromInts(255, 255, 255, 255));
    scene.onDisposeObservable.add(() => {
        currentGraph.dispose();
    });

    let currentActionType = 10;
    const actionTypeNameMeshes: AbstractMesh[] = [];
    for (let id = 0; id < 12; ++id) {
        actionTypeNameMeshes.push(scene.getMeshByName(`text_action_type_${id}`)!);
    }
    const actionTypePoke = new PokeButton(scene, "action_type_button");
    const actionTypeChangedHandler = () => {
        currentActionType = (currentActionType + 1) % actionTypeNameMeshes.length;
        actionTypeNameMeshes.map((mesh, idx) => mesh.isVisible = (idx === currentActionType));
    };
    actionTypeChangedHandler();
    actionTypePoke.onPokeObservable.add(poked => {
        if (poked) {
            actionTypeChangedHandler();
        }
    });

    const sensitivityMesh = scene.getMeshByName("sensitivity_slider")!;
    const sensitivitySlider = SliderBehavior.GetForNode(sensitivityMesh)!;
    const sensitivityGrabbable = PhysicsGrabBehavior.get(sensitivityMesh)!;

    let currentDefinition: ICarlDefinition | undefined = undefined;
    let currentSensitivity = 5;
    const regenerateDefinition = () => {
        currentDefinition?.dispose();
        if (exampleBlocks.size < 1) {
            currentDefinition = undefined;
            currentGraph.recognizer = undefined;
            return;
        }

        const examples: ICarlExample[] = [];
        exampleBlocks.forEach(block => {
            const example = spawner.getExampleFromBlock(block);
            if (example) {
                examples.push(example);
            }
        });
        const counterexamples: ICarlExample[] = [];
        counterexampleBlocks.forEach(block => {
            const counterexample = spawner.getExampleFromBlock(block);
            if (counterexample) {
                counterexamples.push(counterexample);
            }
        });
        currentDefinition = carl.createDefinition(currentActionType, examples, counterexamples);
        currentDefinition.setDefaultSensitivity(currentSensitivity);
        currentGraph.recognizer = carl.createRecognizer(currentDefinition);
        sensitivitySlider.value = currentSensitivity / 10;
    };

    const inEditorMatrix = scene.getMeshByName("editor_volume")!.computeWorldMatrix(true).clone().invert();
    const isExampleMatrix = scene.getMeshByName("examples_volume")!.computeWorldMatrix(true).clone().invert();
    const isCounterexampleMatrix = scene.getMeshByName("counterexamples_volume")!.computeWorldMatrix(true).clone().invert();

    const scratchVec = new Vector3();
    const addExampleBlockPlacementDetection = (block: AbstractMesh) => {
        const blockState = {
            inEditor: false,
            isExample: false,
            isCounterexample: false,
        };
        scene.onBeforeRenderObservable.add(() => {
            Vector3.TransformCoordinatesToRef(block.position, inEditorMatrix, scratchVec);
            let inBounds = Math.abs(scratchVec.x) < 1 && Math.abs(scratchVec.y) < 1 && Math.abs(scratchVec.z) < 1;
            if (!blockState.inEditor && inBounds) {
                previewStopper = previewer.previewExample(spawner.getExampleFromBlock(block)!);
                blockState.inEditor = true;
                return;
            } else if (blockState.inEditor && !inBounds) {
                previewStopper();
                previewStopper = null;
                blockState.inEditor = false;
            }
            
            Vector3.TransformCoordinatesToRef(block.position, isExampleMatrix, scratchVec);
            inBounds = Math.abs(scratchVec.x) < 1 && Math.abs(scratchVec.y) < 1 && Math.abs(scratchVec.z) < 1;
            if (!blockState.isExample && inBounds) {
                exampleBlocks.add(block);
                regenerateDefinition();
                blockState.isExample = true;
                return;
            } else if (blockState.isExample && !inBounds) {
                exampleBlocks.delete(block);
                regenerateDefinition();
                blockState.isExample = false;
            }
            
            Vector3.TransformCoordinatesToRef(block.position, isCounterexampleMatrix, scratchVec);
            inBounds = Math.abs(scratchVec.x) < 1 && Math.abs(scratchVec.y) < 1 && Math.abs(scratchVec.z) < 1;
            if (!blockState.isCounterexample && inBounds) {
                counterexampleBlocks.add(block);
                regenerateDefinition();
                blockState.isCounterexample = true;
                return;
            } else if (blockState.isCounterexample && !inBounds) {
                counterexampleBlocks.delete(block);
                regenerateDefinition();
                blockState.isCounterexample = false;
            }
        });
    }

    sensitivityGrabbable.onGrabEndedObservable.add(() => {
        currentSensitivity = 10 * sensitivitySlider.value;
        regenerateDefinition();
    });

    const recordingPokeButton = new PokeButton(scene, "recording_start", "recording_stop");
    let recording: number | undefined = undefined;
    let example: ICarlExample | undefined = undefined;
    recordingPokeButton.onPokeObservable.add(poked => {
        if (poked) {
            if (recording) {
                example = carl.stopRecording(recording);
                recording = undefined;
                
                const block = spawner.spawnNewExampleBlock(example);
                addExampleBlockPlacementDetection(block);
            } else {
                recording = carl.startRecording();
            }
        }
    });

    const playPoke = new PokeButton(scene, "example_play");
    playPoke.onPokeObservable.add(poked => {
        if (poked) {
            previewer.play();
        }
    });

    const downloadPoke = new PokeButton(scene, "download_button");
    downloadPoke.onPokeObservable.add(poked => {
        if (poked && currentDefinition) {
            currentDefinition.download();
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

    const exitPoke = new PokeButton(scene, "exit_xr_button");
    exitPoke.onPokeObservable.addOnce(() => {
        xr.input.xrSessionManager.exitXRAsync();
        // TODO: Dispose scene? Or reset?
    });
}

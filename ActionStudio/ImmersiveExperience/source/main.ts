/**
 * Entry point for the 3D immersive XR experience.
 *
 * `initializeImmersiveExperienceAsync` bootstraps the full Babylon.js / WebXR session:
 * loads the glTF scene, enables Havok physics, sets up hand tracking, wires up all
 * interactive UI elements (sliders, poke buttons, block spawners, recognition graph),
 * and starts the per-frame input-sample loop that feeds hand poses into CARL.
 */
import { _InstancesBatch, Color4, Engine, FreeCamera, HavokPlugin, HemisphericLight, MeshBuilder, TransformNode as AbstractMesh, Vector3, WebXRControllerPointerSelection, WebXRFeatureName, TransformNode, Color3, StandardMaterial, Mesh } from "@babylonjs/core";
import HavokPhysics from "@babylonjs/havok";
import "@babylonjs/loaders/glTF/2.0";
import { ICarl, ICarlDefinition, ICarlExample, ICarlInputSample, IInitialBlock } from "./carlInterfaces";
import { HandPinchGrabber } from "./handPinchGrabber";
import { PhysicsEnabledScene } from "./physicsEnabledScene";
import { PhysicsGrabBehavior } from "./physicsGrabBehavior";
import { OPENXR_JOINT_MAPPINGS, populateInputSample } from "./utils";
import { PokeButton } from "./pokeButton";
import { ExamplePreviewer } from "./examplePreviewer";
import { InputPuppet } from "./inputPuppet";
import { RecognitionGraph } from "./recognitionGraph";
import { SliderBehavior } from "./slider";
import { BlockSpawner } from "./blockSpawner";

export interface IImmersiveExperience {
    dispose(): void;
}

export interface IImmersiveExperienceOptions {
    initialExamples?: IInitialBlock<ICarlExample>[];
    initialDefinitions?: IInitialBlock<ICarlDefinition>[];
    onSessionEnded?: () => void;
}

export async function initializeImmersiveExperienceAsync(canvas: HTMLCanvasElement, carl: ICarl, options?: IImmersiveExperienceOptions): Promise<IImmersiveExperience> {
    const initialExamples = options?.initialExamples ?? [];
    const initialDefinitions = options?.initialDefinitions ?? [];
    const onSessionEnded = options?.onSessionEnded;

    // --- Engine and scene ---
    const engine = new Engine(canvas, true);
    engine.displayLoadingUI = () => {};
    engine.hideLoadingUI();

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

    // --- UI state ---
    const exampleSpawner = new BlockSpawner<ICarlExample>(scene, "example_block");
    const definitionSpawner = new BlockSpawner<ICarlDefinition>(scene, "definition_block");

    // --- WebXR session ---
    const xr = await scene.createDefaultXRExperienceAsync({
        uiOptions: {
            sessionMode: "immersive-ar",
        },
        disableDefaultUI: true,
        disablePointerSelection: true,
        disableTeleportation: true,
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

    // --- Hand tracking ---
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

    // --- Preview and recognition state ---
    const previewer = new ExamplePreviewer(scene);
    let previewStopper: any = null;

    const exampleBlocks = new Set<TransformNode>();
    const counterexampleBlocks = new Set<TransformNode>();

    const currentGraph = new RecognitionGraph(scene, Color4.FromInts(255, 255, 255, 255));

    let currentActionType = 10;

    // --- Color and sensitivity sliders ---
    const rMesh = scene.getMeshByName("r_slider")!;
    const rSlider = SliderBehavior.GetForNode(rMesh)!;
    const gMesh = scene.getMeshByName("g_slider")!;
    const gSlider = SliderBehavior.GetForNode(gMesh)!;
    const bMesh = scene.getMeshByName("b_slider")!;
    const bSlider = SliderBehavior.GetForNode(bMesh)!;

    const colorMat = new StandardMaterial("colorMat", scene);
    const handleColorChanged = () => {
        colorMat.diffuseColor.set(rSlider.value, gSlider.value, bSlider.value);
    };
    rSlider.onUpdatedObservable.add(handleColorChanged);
    gSlider.onUpdatedObservable.add(handleColorChanged);
    bSlider.onUpdatedObservable.add(handleColorChanged);
    scene.getMeshByName("color_swatch")!.material = colorMat;
    handleColorChanged();

    const sensitivityMesh = scene.getMeshByName("sensitivity_slider")!;
    const sensitivitySlider = SliderBehavior.GetForNode(sensitivityMesh)!;
    const sensitivityGrabbable = PhysicsGrabBehavior.get(sensitivityMesh)!;

    // --- Definition regeneration ---
    let draftDefinition: ICarlDefinition | undefined = undefined;
    let currentSensitivity = 5;
    const regenerateDefinition = () => {
        draftDefinition?.dispose();
        if (exampleBlocks.size < 1) {
            draftDefinition = undefined;
            currentGraph.recognizer = undefined;
            return;
        }

        const examples: ICarlExample[] = [];
        exampleBlocks.forEach(block => {
            const example = exampleSpawner.getValueFromBlock(block);
            if (example) {
                examples.push(example);
            }
        });
        const counterexamples: ICarlExample[] = [];
        counterexampleBlocks.forEach(block => {
            const counterexample = exampleSpawner.getValueFromBlock(block);
            if (counterexample) {
                counterexamples.push(counterexample);
            }
        });
        draftDefinition = carl.draftDefinition(currentActionType, examples, counterexamples);
        draftDefinition.setDefaultSensitivity(currentSensitivity);
        currentGraph.recognizer = carl.createRecognizer(draftDefinition);
        sensitivitySlider.value = currentSensitivity / 10;
    };

    // --- Spatial volume matrices ---
    const isDiscardedMatrix = scene.getMeshByName("discard_volume")!.computeWorldMatrix(true).clone().invert();
    const inEditorMatrix = scene.getMeshByName("editor_volume")!.computeWorldMatrix(true).clone().invert();
    const isExampleMatrix = scene.getMeshByName("examples_volume")!.computeWorldMatrix(true).clone().invert();
    const isCounterexampleMatrix = scene.getMeshByName("counterexamples_volume")!.computeWorldMatrix(true).clone().invert();
    const isDemoMatrix = scene.getMeshByName("demo_volume")!.computeWorldMatrix(true).clone().invert();

    // --- Block placement detection ---
    const scratchVec = new Vector3();
    const addExampleBlockPlacementDetection = (block: AbstractMesh) => {
        const blockState = {
            inEditor: false,
            isExample: false,
            isCounterexample: false,
        };
        scene.onBeforeRenderObservable.add(() => {
            Vector3.TransformCoordinatesToRef(block.position, isDiscardedMatrix, scratchVec);
            let inBounds = Math.abs(scratchVec.x) < 1 && Math.abs(scratchVec.y) < 1 && Math.abs(scratchVec.z) < 1;
            if (inBounds) {
                exampleSpawner.getValueFromBlock(block)?.dispose();
                block.dispose();
                return;
            }

            Vector3.TransformCoordinatesToRef(block.position, inEditorMatrix, scratchVec);
            inBounds = Math.abs(scratchVec.x) < 1 && Math.abs(scratchVec.y) < 1 && Math.abs(scratchVec.z) < 1;
            if (!blockState.inEditor && inBounds) {
                previewStopper = previewer.previewExample(exampleSpawner.getValueFromBlock(block)!);
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

    const addDefinitionBlockPlacementDetection = (block: AbstractMesh) => {
        interface IBlockState {
            graph: RecognitionGraph | undefined;
        }
        const blockState: IBlockState = {
            graph: undefined,
        };
        scene.onBeforeRenderObservable.add(() => {
            Vector3.TransformCoordinatesToRef(block.position, isDiscardedMatrix, scratchVec);
            let inBounds = Math.abs(scratchVec.x) < 1 && Math.abs(scratchVec.y) < 1 && Math.abs(scratchVec.z) < 1;
            if (inBounds) {
                definitionSpawner.getValueFromBlock(block)?.dispose();
                block.dispose();
                return;
            }

            Vector3.TransformCoordinatesToRef(block.position, isDemoMatrix, scratchVec);
            inBounds = Math.abs(scratchVec.x) < 1 && Math.abs(scratchVec.y) < 1 && Math.abs(scratchVec.z) < 1;
            if (!blockState.graph && inBounds) {
                const definition = definitionSpawner.getValueFromBlock(block)!;
                blockState.graph = new RecognitionGraph(scene, Color4.FromColor3(((block as Mesh).material as StandardMaterial).diffuseColor));
                blockState.graph.recognizer = carl.createRecognizer(definition);
                return;
            } else if (blockState.graph && !inBounds) {
                blockState.graph.dispose();
                blockState.graph = undefined;
            }
        });
    }

    // --- Poke buttons ---
    sensitivityGrabbable.onGrabEndedObservable.add(() => {
        currentSensitivity = 10 * sensitivitySlider.value;
        regenerateDefinition();
    });

    const actionTypeNameMeshes: AbstractMesh[] = [];
    for (let id = 0; id < 12; ++id) {
        actionTypeNameMeshes.push(scene.getMeshByName(`text_action_type_${id}`)!);
    }
    const actionTypePoke = new PokeButton(scene, "action_type_button");
    const actionTypeChangedHandler = () => {
        currentActionType = (currentActionType + 1) % actionTypeNameMeshes.length;
        actionTypeNameMeshes.map((mesh, idx) => mesh.isVisible = (idx === currentActionType));
        regenerateDefinition();
    };
    actionTypeChangedHandler();
    actionTypePoke.onPokeObservable.add(poked => {
        if (poked) {
            actionTypeChangedHandler();
        }
    });

    const recordingPokeButton = new PokeButton(scene, "recording_start", "recording_stop");
    let recording: number | undefined = undefined;
    let example: ICarlExample | undefined = undefined;
    recordingPokeButton.onPokeObservable.add(poked => {
        if (poked) {
            if (recording) {
                const colorString = colorMat.diffuseColor.toHexString();
                example = carl.stopRecording(recording, { color: colorString });
                recording = undefined;

                const block = exampleSpawner.spawnNewBlock(example);
                block.material = colorMat.clone(`${block.name}_mat`);
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
        if (poked && draftDefinition) {
            draftDefinition.download();
        }
    });

    const createDefinitionPoke = new PokeButton(scene, "definition_button");
    createDefinitionPoke.onPokeObservable.add(poked => {
        if (poked && draftDefinition) {
            const colorString = colorMat.diffuseColor.toHexString();
            carl.finalizeDefinition(draftDefinition, { color: colorString });
            const definitionBlock = definitionSpawner.spawnNewBlock(draftDefinition);
            definitionBlock.material = colorMat.clone(`${definitionBlock.name}_mat`);
            addDefinitionBlockPlacementDetection(definitionBlock);

            draftDefinition = undefined;
            regenerateDefinition();
        }
    });

    // --- Input sample loop ---
    xr.input.xrSessionManager.onXRFrameObservable.add((frame) => {
        let sample = createInputSample();
        populateInputSample(xr.input.xrCamera, scene.leftHand, scene.rightHand, sample);
        carl.handleInputSample(sample);

        if (scene.inputPuppet === undefined && scene.leftHand && scene.rightHand && scene.leftHand.getJointMesh(OPENXR_JOINT_MAPPINGS[0]).scaling.x < 1) {
            scene.inputPuppet = new InputPuppet(scene.leftHand, scene.rightHand);
        }
    });

    // --- Disposal ---
    let disposed = false;
    const shutdown = async () => {
        if (disposed) return;
        disposed = true;

        try {
            await xr.baseExperience.exitXRAsync();
        } catch (e) {
            // XR session may already be ended
        }

        engine.stopRenderLoop();
        scene.dispose();
        engine.dispose();
        onSessionEnded?.();
    };

    const exitPoke = new PokeButton(scene, "exit_xr_button");
    exitPoke.onPokeObservable.add(() => {
        shutdown();
    });

    // --- Auto-spawn initial example blocks ---
    for (const initial of initialExamples) {
        const mat = new StandardMaterial(`initial_example_mat`, scene);
        mat.diffuseColor = Color3.FromHexString(initial.color);
        const block = exampleSpawner.spawnNewBlock(initial.value);
        block.material = mat;
        addExampleBlockPlacementDetection(block);
    }

    // --- Auto-spawn initial definition blocks ---
    for (const initial of initialDefinitions) {
        const mat = new StandardMaterial(`initial_definition_mat`, scene);
        mat.diffuseColor = Color3.FromHexString(initial.color);
        const block = definitionSpawner.spawnNewBlock(initial.value);
        block.material = mat;
        addDefinitionBlockPlacementDetection(block);
    }

    // --- Enter XR automatically ---
    xr.baseExperience.enterXRAsync("immersive-ar", "local-floor");

    return {
        dispose: () => { shutdown(); },
    };
}

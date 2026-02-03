import { _InstancesBatch, AbstractMesh, Color3, Engine, FreeCamera, HavokPlugin, HemisphericLight, Mesh, MeshBuilder, PBRMaterial, Quaternion, Scene, Tools, Vector3, WebXRCamera, WebXRControllerPointerSelection, WebXRFeatureName, WebXRHand, WebXRHandJoint, WebXRState } from "@babylonjs/core";
import HavokPhysics from "@babylonjs/havok";
import "@babylonjs/loaders/glTF/2.0";
import { ICarl, ICarlInputSample } from "./carlInterfaces";
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
function populateInputSample(hmd: WebXRCamera, leftHand: WebXRHand | null, rightHand: WebXRHand | null, sample: ICarlInputSample): void {
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

    const hands = new Map<XRHandedness, WebXRHand>();
    handTracking.onHandAddedObservable.add((hand) => {
        hands.set(hand.xrController.inputSource.handedness, hand);

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

        const handedness = hand.xrController.inputSource.handedness;
        if (hands.has(handedness) && hands.get(handedness)?.xrController.uniqueId == hand.xrController.uniqueId) {
            hands.delete(handedness);
        }
    });

    const cachedInputSample = carl.createInputSample();
    function createInputSample(): ICarlInputSample {
        return {...cachedInputSample};
    }

    xr.input.xrSessionManager.onXRFrameObservable.add((frame) => {
        let sample = createInputSample();
        populateInputSample(xr.input.xrCamera, hands.has("left") ? hands.get("left")! : null, hands.has("right") ? hands.get("right")! : null, sample);
        carl.handleInputSample(sample);
    });

    // Testing CARL functionality sans interaction (yet)
    xr.input.xrSessionManager.onXRFrameObservable.addOnce(async () => {
        await Tools.DelayAsync(2000);
        const recId = carl.startRecording();
        await Tools.DelayAsync(500);
        const example = carl.stopRecording(recId);
        const definition = carl.createDefinition(2, [example], []);
        const recognizer = carl.createRecognizer(definition);
        
        example.dispose();
        definition.dispose();
    
        scene.onBeforeRenderObservable.add(() => {
            console.log("Current score: " + recognizer.currentScore());
        });
    });
}

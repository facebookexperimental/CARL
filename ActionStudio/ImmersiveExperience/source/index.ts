import { Clamp, Engine, FreeCamera, HavokPlugin, HemisphericLight, Matrix, MeshBuilder, Nullable, Observable, Observer, PhysicsAggregate, PhysicsMotionType, PhysicsShapeType, Quaternion, Scene, TransformNode, Vector3, WebXRControllerPointerSelection, WebXRFeatureName, WebXRHand, WebXRHandJoint, WebXRSessionManager } from "@babylonjs/core";
import HavokPhysics from "@babylonjs/havok";

class HandPinchGrabber extends TransformNode {
    public onGrabChanged: Observable<boolean> = new Observable<boolean>();

    private _isGrabbing: boolean = false;
    private _grabFrameCount: number = 0;
    private static readonly GRAB_FRAMES_TO_CONFIRM: number = 4;
    private _pinchThreshold: number = 0.03;
    private _frameObserver: Observer<XRFrame>;
    private _sessionManager: WebXRSessionManager;

    constructor(name: string, hand: WebXRHand, sessionManager: WebXRSessionManager, scene: Scene) {
        super(name, scene);
        this.rotationQuaternion = Quaternion.Identity();

        this._sessionManager = sessionManager;

        const indexPosition = new Vector3();
        const thumbPosition = new Vector3();
        this._frameObserver = this._sessionManager.onXRFrameObservable.add(() => {
            let index = hand.getJointMesh(WebXRHandJoint.INDEX_FINGER_TIP)!;
            let thumb = hand.getJointMesh(WebXRHandJoint.THUMB_TIP)!;
            indexPosition.copyFrom(index.absolutePosition);
            thumbPosition.copyFrom(thumb.absolutePosition);

            this.position.copyFrom(indexPosition);
            this.position.addInPlace(thumbPosition);
            this.position.scaleInPlace(0.5);

            Quaternion.SlerpToRef(index.rotationQuaternion!, thumb.rotationQuaternion!, 0.5, this.rotationQuaternion!);

            if (Vector3.Distance(indexPosition, thumbPosition) > this._pinchThreshold * 1.2) {
                this._grabFrameCount = Math.max(-HandPinchGrabber.GRAB_FRAMES_TO_CONFIRM, Math.min(this._grabFrameCount - 1, -1));
            } else if (!this._isGrabbing && Vector3.Distance(indexPosition, thumbPosition) <= this._pinchThreshold) {
                this._grabFrameCount = Math.max(1, Math.min(this._grabFrameCount + 1, HandPinchGrabber.GRAB_FRAMES_TO_CONFIRM));
            }

            if (this._isGrabbing && this._grabFrameCount == -HandPinchGrabber.GRAB_FRAMES_TO_CONFIRM) {
                this._isGrabbing = false;
                this.onGrabChanged.notifyObservers(this._isGrabbing);
            } else if (!this._isGrabbing && this._grabFrameCount == HandPinchGrabber.GRAB_FRAMES_TO_CONFIRM) {
                this._isGrabbing = true;
                this.onGrabChanged.notifyObservers(this._isGrabbing);
            }
        });
    }

    public dispose(doNotRecurse?: boolean, disposeMaterialAndTextures?: boolean): void {
        super.dispose(doNotRecurse, disposeMaterialAndTextures);
        this._sessionManager.onXRFrameObservable.remove(this._frameObserver);
    }
}

export function initializeImmersiveExperience(canvas: HTMLCanvasElement): void {
    const engine = new Engine(canvas, true);

    const resizeHandler = () => {
        engine.resize();
    };
    window.addEventListener("resize", resizeHandler);
    engine.onDisposeObservable.add(() => {
        window.removeEventListener("resize", resizeHandler);
    });

    const scene = new Scene(engine);
    const camera = new FreeCamera("camera1", new Vector3(0, 1, 0), scene);
    camera.attachControl(canvas, true);

    const light = new HemisphericLight("light", new Vector3(0, 1, 0), scene);
    light.intensity = 0.7;

    const ground = MeshBuilder.CreateGround("ground", {width: 10, height: 10}, scene);
    const table = MeshBuilder.CreateBox("table", {width: 2, depth: 1, height: 0.02}, scene);
    table.position.set(0, 1, 1);
    const recording = MeshBuilder.CreateBox("recording", {width: 0.04, depth: 0.04, height: 0.04}, scene);
    recording.position.set(0, 1.1, 1);

    const recording2 = MeshBuilder.CreateBox("recording2", {width: 0.04, depth: 0.04, height: 0.04}, scene);
    recording2.position.set(0.5, 1.1, 1);
    
    HavokPhysics().then((havok) => {
        const hk = new HavokPlugin(true, havok);
        scene.enablePhysics(new Vector3(0, -9.8, 0), hk);
        
        const groundAggregate = new PhysicsAggregate(ground, PhysicsShapeType.BOX, { mass: 0 }, scene);
        const tableAggregate = new PhysicsAggregate(table, PhysicsShapeType.BOX, { mass: 0 }, scene);
        const recordingAggregate = new PhysicsAggregate(recording, PhysicsShapeType.BOX, { mass: 1 }, scene);
        recordingAggregate.transformNode.rotationQuaternion = Quaternion.Identity();
        const recording2Aggregate = new PhysicsAggregate(recording2, PhysicsShapeType.BOX, { mass: 1 }, scene);
    });
    
    scene.createDefaultXRExperienceAsync().then((xr) => {
        xr.baseExperience.featuresManager.disableFeature(WebXRControllerPointerSelection.Name);

        const handTracking = xr.baseExperience.featuresManager.enableFeature(WebXRFeatureName.HAND_TRACKING, "latest", {
            xrInput: xr.input,
            jointMeshes: {
                sourceMesh: MeshBuilder.CreateBox("jointParent", { size: 1 }),
                keepOriginalVisible: true,
            },
        });

        handTracking.onHandAddedObservable.add((hand) => {
            const grabber = new HandPinchGrabber("handGrabber", hand, xr.input.xrSessionManager, scene);
            let observer: Nullable<Observer<XRFrame>> = null;
            const offsetMatrix = new Matrix();
            const worldMatrix = new Matrix();
            grabber.onGrabChanged.add((isGrabbing) => {
                if (isGrabbing && Vector3.Distance(grabber.position, recording.position) < 0.04) {
                    recording.physicsBody!.setMotionType(PhysicsMotionType.ANIMATED);
                    recording.physicsBody!.disablePreStep = false;
                    recording.getWorldMatrix().multiplyToRef(grabber.getWorldMatrix().invert(), offsetMatrix);
                    observer = xr.input.xrSessionManager.onXRFrameObservable.add(() => {
                        offsetMatrix.multiplyToRef(grabber.getWorldMatrix(), worldMatrix);
                        worldMatrix.decompose(undefined, recording.physicsBody!.transformNode.rotationQuaternion!, recording.physicsBody!.transformNode.position);
                    });
                } else if (!isGrabbing) {
                    recording.physicsBody!.disablePreStep = true;
                    recording.physicsBody!.setMotionType(PhysicsMotionType.DYNAMIC);
                    if (observer) {
                        xr.input.xrSessionManager.onXRFrameObservable.remove(observer);
                        observer = null;
                    }
                }
            });
        });

        handTracking.onHandRemovedObservable.add((hand) => {
            console.log("Hand removed:", hand);
        });
    });

    engine.runRenderLoop(() => {
        scene.render();
    });
}

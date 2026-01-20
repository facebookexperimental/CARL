import { _InstancesBatch, Engine, FreeCamera, HavokPlugin, HemisphericLight, IDisposable, Matrix, MeshBuilder, Observable, Observer, PhysicsAggregate, PhysicsMotionType, PhysicsShapeType, Quaternion, Scene, TransformNode, Vector3, WebXRControllerPointerSelection, WebXRFeatureName, WebXRHand, WebXRHandJoint, WebXRSessionManager, WebXRState } from "@babylonjs/core";
import HavokPhysics from "@babylonjs/havok";

interface IGrabber extends TransformNode {
    onGrabChanged: Observable<boolean>;
    get isGrabbing(): boolean;
}

interface IProposal {
    score: number;
    accept: () => void;
}

class HandPinchGrabber extends TransformNode implements IGrabber {
    public onGrabChanged: Observable<boolean> = new Observable<boolean>();

    private _isGrabbing: boolean = false;
    public get isGrabbing(): boolean {
      return this._isGrabbing;  
    }
    private _grabFrameCount: number = 0;
    private static readonly GRAB_FRAMES_TO_CONFIRM: number = 4;
    private _pinchThreshold: number = 0.02;
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
            thumbPosition.scaleAndAddToRef(3, this.position);
            this.position.scaleInPlace(0.25);
            this.computeWorldMatrix();

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

class PhysicsGrabBehavior implements IDisposable {
    private readonly _node: TransformNode;
    private readonly _offsetMatrix = Matrix.Identity();
    private readonly _worldMatrix = Matrix.Identity();
    private _observer: Observer<TransformNode> | null = null;
    private _currentGrabber: IGrabber | null = null;

    private constructor(node: TransformNode) {
        this._node = node;
        node.onDisposeObservable.addOnce(() => this.dispose());
    }

    private static _behaviors = new Map<number, PhysicsGrabBehavior>();

    public static attach(node: TransformNode): PhysicsGrabBehavior {
        if (!PhysicsGrabBehavior._behaviors.has(node.uniqueId)) {
            PhysicsGrabBehavior._behaviors.set(node.uniqueId, new PhysicsGrabBehavior(node));
        }
        return PhysicsGrabBehavior._behaviors.get(node.uniqueId)!;
    }

    public static handleGrab(grabber: IGrabber) {
        let bestProposal: IProposal | null = null;
        PhysicsGrabBehavior._behaviors.forEach((grabbable) => {
            const p = grabbable.proposeGrab(grabber);
            if (bestProposal === null || p.score > bestProposal.score) {
                bestProposal = p;
            }
        });
        if (bestProposal && (bestProposal as IProposal).score > 0) {
            (bestProposal as IProposal).accept();
        }
    }

    public proposeGrab(grabber: IGrabber): IProposal {
        if (grabber.isGrabbing) {
            const DISTANCE = 0.04;
            const score = Math.min(Math.max((DISTANCE - Vector3.Distance(grabber.position, this._node.position)) / DISTANCE, 0), 1);
            return {
                score: score,
                accept: () => {
                    if (this._currentGrabber !== null) {
                        // If we're already being grabbed, we need to stop observing that grab, but we're already set up for animated physics.
                        this._observer?.remove();
                    } else {
                        // Otherwise, we need to change our physics mode from to permit animation.
                        this._node.physicsBody!.setMotionType(PhysicsMotionType.ANIMATED);
                        this._node.physicsBody!.disablePreStep = false;
                    }
                    this._currentGrabber = grabber;
                    this._node.getWorldMatrix().multiplyToRef(grabber.getWorldMatrix().invert(), this._offsetMatrix);
                    this._observer = this._currentGrabber.onAfterWorldMatrixUpdateObservable.add(() => {
                        this._offsetMatrix.multiplyToRef(grabber.getWorldMatrix(), this._worldMatrix);
                        this._worldMatrix.decompose(undefined, this._node.physicsBody!.transformNode.rotationQuaternion!, this._node.physicsBody!.transformNode.position);
                    });
                    // TODO: Handle what happens if the grabber is disposed.
                }
            }
        } else if (this._currentGrabber === grabber && !grabber.isGrabbing) {
            return {
                score: 1,
                accept: () => {
                    this._currentGrabber = null;
                    this._node.physicsBody!.disablePreStep = true;
                    this._node.physicsBody!.setMotionType(PhysicsMotionType.DYNAMIC);
                    this._observer?.remove();
                    this._observer = null;
                }
            }
        } else {
            return {score: 0, accept: () => {}};
        }
    }

    public dispose(): void {
        PhysicsGrabBehavior._behaviors.delete(this._node.uniqueId);
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
    scene.useRightHandedSystem = true;
    const camera = new FreeCamera("camera1", new Vector3(0, 1, 0), scene);
    camera.attachControl(canvas, true);

    const light = new HemisphericLight("light", new Vector3(0, 1, 0), scene);
    light.intensity = 0.7;

    const ground = MeshBuilder.CreateGround("ground", {width: 10, height: 10}, scene);
    const table = MeshBuilder.CreateBox("table", {width: 2, depth: 1, height: 0.02}, scene);
    table.position.set(0, 1, 1);
    const recording = MeshBuilder.CreateBox("recording", {width: 0.04, depth: 0.04, height: 0.04}, scene);
    recording.position.set(0, 1.1, 1);
    const grabbables: PhysicsGrabBehavior[] = [];
    const grabbers = new Map<string, HandPinchGrabber>();

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
        grabbables.push(PhysicsGrabBehavior.attach(recording));
        grabbables.push(PhysicsGrabBehavior.attach(recording2));
    });
    
    scene.createDefaultXRExperienceAsync().then((xr) => {
        xr.baseExperience.onStateChangedObservable.add((state) => {
            if (state === WebXRState.IN_XR) {
                canvas.style.display = "none";
            } else {
                canvas.style.display = "block";
            }
        });

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
            grabber.onGrabChanged.add(() => {
                PhysicsGrabBehavior.handleGrab(grabber);
            });
            grabbers.set(hand.xrController.uniqueId, grabber);
        });

        handTracking.onHandRemovedObservable.add((hand) => {
            if (grabbers.has(hand.xrController.uniqueId)) {
                grabbers.get(hand.xrController.uniqueId)?.dispose();
                grabbers.delete(hand.xrController.uniqueId);
            }
        });
    });

    engine.runRenderLoop(() => {
        scene.render();
    });
}

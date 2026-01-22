import { _InstancesBatch, AbstractMesh, AppendSceneAsync, Engine, FreeCamera, HavokPlugin, HemisphericLight, IDisposable, Matrix, MeshBuilder, Node, Observable, Observer, PhysicsAggregate, PhysicsMotionType, PhysicsShape, PhysicsShapeType, Quaternion, Scene, TransformNode, Vector3, WebXRControllerPointerSelection, WebXRFeatureName, WebXRHand, WebXRHandJoint, WebXRSessionManager, WebXRState } from "@babylonjs/core";
import HavokPhysics from "@babylonjs/havok";
import "@babylonjs/loaders/glTF/2.0";

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
    private _grabberObserver: Observer<TransformNode> | null = null;
    private _currentGrabber: IGrabber | null = null;
    private _currentGrabberDisposeObserver: Observer<Node> | null = null;
    private _disposeObserver: Observer<Node>;

    private constructor(node: TransformNode) {
        this._node = node;
        this._disposeObserver = node.onDisposeObservable.add(() => this.dispose());
        node.onDisposeObservable.addOnce(() => this.dispose());
    }

    private static _behaviors = new Map<number, PhysicsGrabBehavior>();

    public static attach(node: TransformNode): PhysicsGrabBehavior {
        if (!PhysicsGrabBehavior._behaviors.has(node.uniqueId)) {
            const behavior = new PhysicsGrabBehavior(node);
            PhysicsGrabBehavior._behaviors.set(node.uniqueId, behavior);
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

    private _pausePhysics(): void {
        this._node.physicsBody!.setMotionType(PhysicsMotionType.ANIMATED);
        this._node.physicsBody!.disablePreStep = false;
    }

    private _resumePhysics(): void {
        this._currentGrabber = null;
        this._node.physicsBody!.disablePreStep = true;
        this._node.physicsBody!.setMotionType(PhysicsMotionType.DYNAMIC);
        this._grabberObserver?.remove();
        this._grabberObserver = null;
        this._currentGrabberDisposeObserver?.remove();
        this._currentGrabberDisposeObserver = null;
    }

    public proposeGrab(grabber: IGrabber): IProposal {
        if (grabber.isGrabbing) {
            const k = 1.2 * this._node.physicsBody!.getBoundingBox().extendSizeWorld!.length();
            const score = Math.min(Math.max((k - Vector3.Distance(grabber.position, this._node.position)) / k, 0), 1);
            return {
                score: score,
                accept: () => {
                    if (this._currentGrabber !== null) {
                        // If we're already being grabbed, we need to stop observing that grab, but we're already set up for animated physics.
                        this._grabberObserver?.remove();
                    } else {
                        // Otherwise, we need to change our physics mode from to permit animation.
                        this._pausePhysics();
                    }
                    this._currentGrabber = grabber;
                    this._node.getWorldMatrix().multiplyToRef(grabber.getWorldMatrix().invert(), this._offsetMatrix);
                    this._grabberObserver = this._currentGrabber.onAfterWorldMatrixUpdateObservable.add(() => {
                        this._offsetMatrix.multiplyToRef(grabber.getWorldMatrix(), this._worldMatrix);
                        this._worldMatrix.decompose(undefined, this._node.physicsBody!.transformNode.rotationQuaternion!, this._node.physicsBody!.transformNode.position);
                    });
                    const grabberDisposedObserver = this._currentGrabber.onDisposeObservable.add(() => {
                        this._resumePhysics();
                    });
                    // TODO: Handle what happens if the grabber is disposed.
                }
            }
        } else if (this._currentGrabber === grabber && !grabber.isGrabbing) {
            return {
                score: 1,
                accept: () => this._resumePhysics()
            };
        } else {
            return { score: 0, accept: () => { } };
        }
    }

    public dispose(): void {
        PhysicsGrabBehavior._behaviors.delete(this._node.uniqueId);
        this._disposeObserver?.remove();
    }
}

class PhysicsEnabledScene extends Scene {
    public grabbables: PhysicsGrabBehavior[] = [];
    public grabbers = new Map<string, HandPinchGrabber>();

    private constructor(engine: Engine) {
        super(engine);
    }

    public static async loadAsync(url: string, engine: Engine, havok: HavokPlugin): Promise<PhysicsEnabledScene> {
        const scene = new PhysicsEnabledScene(engine);
        scene.useRightHandedSystem = true;
        await AppendSceneAsync(url, scene);
        scene.enablePhysics(Vector3.DownReadOnly.scale(9.81), havok);

        const physicsNulls: TransformNode[] =  [];

        scene.transformNodes.forEach(node => {
            if (node.name.startsWith("physics")) {
                physicsNulls.push(node);
            }
        });

        physicsNulls.forEach(physicsNull => {
            const mesh = physicsNull.parent! as AbstractMesh;
            const shapeType = physicsNull.name.split(';')[0];
            let shape: PhysicsShapeType;
            switch (shapeType) {
                case "physics_sphere":
                    shape = PhysicsShapeType.SPHERE;
                    break;
                case "physics_mesh":
                    shape = PhysicsShapeType.MESH;
                    break;
                default:
                    shape = PhysicsShapeType.BOX;
                    break;
            }
            physicsNull.setParent(null);
            physicsNull.dispose();

            let mass = 0;
            const massNulls = mesh.getChildTransformNodes(true, child => child.name.startsWith("mass:"));
            if (massNulls.length > 0) {
                const massNull = massNulls[0];
                mass = parseFloat(massNull.name.split(';')[0].split(':')[1]);
                massNull.setParent(null);
                massNull.dispose();
            }
            
            let grabbable = false;
            const grabbableNulls = mesh.getChildTransformNodes(true, child => child.name.startsWith("grabbable"));
            if (grabbableNulls.length > 0) {
                grabbable = true;
                const grabbableNull = grabbableNulls[0];
                grabbableNull.setParent(null);
                grabbableNull.dispose();
            }

            const aggregate = new PhysicsAggregate(mesh, shape, { mass: mass }, scene);
            aggregate.transformNode.rotationQuaternion = Quaternion.Identity();
            if (grabbable) {
                scene.grabbables.push(PhysicsGrabBehavior.attach(mesh));
            }
        });

        return scene;
    }
}

export async function initializeImmersiveExperienceAsync(canvas: HTMLCanvasElement): Promise<void> {
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

    const handTracking = xr.baseExperience.featuresManager.enableFeature(WebXRFeatureName.HAND_TRACKING, "latest", {
        xrInput: xr.input,
        jointMeshes: {
            sourceMesh: MeshBuilder.CreateBox("jointParent", { size: 1 }),
            //keepOriginalVisible: true,
        },
    });

    handTracking.onHandAddedObservable.add((hand) => {
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
    });

    engine.runRenderLoop(() => {
        scene.render();
    });
}

import { IDisposable, Matrix, Node, Observable, Observer, PhysicsMotionType, TransformNode, Vector3 } from "@babylonjs/core";
import { IGrabber, IProposal } from "./handPinchGrabber";

export class PhysicsGrabBehavior implements IDisposable {
    private readonly _node: TransformNode;
    private readonly _offsetMatrix = Matrix.Identity();
    private readonly _worldMatrix = Matrix.Identity();
    private _grabberObserver: Observer<TransformNode> | null = null;
    private _currentGrabber: IGrabber | null = null;
    private _currentGrabberDisposeObserver: Observer<Node> | null = null;
    private _disposeObserver: Observer<Node>;

    public onGrabMovedObservable: Observable<TransformNode> = new Observable<TransformNode>();

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
                        this.onGrabMovedObservable.notifyObservers(this._node.physicsBody!.transformNode);
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

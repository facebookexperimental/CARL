/**
 * HandPinchGrabber — detects a pinch gesture on a WebXR hand and emits grab events.
 *
 * Tracks the distance between the index fingertip and thumb tip each XR frame.
 * Uses hysteresis (separate enter/exit thresholds) and a 4-frame confirmation window
 * to avoid spurious grab transitions caused by tracking jitter.
 */
import { Observable, Observer, Quaternion, Scene, TransformNode, Vector3, WebXRHand, WebXRHandJoint, WebXRSessionManager } from "@babylonjs/core";
import { PhysicsEnabledScene } from "./physicsEnabledScene";


export interface IGrabber extends TransformNode {
    onGrabChanged: Observable<boolean>;
    get isGrabbing(): boolean;
}

export interface IProposal {
    score: number;
    accept: () => void;
}

export class HandPinchGrabber extends TransformNode implements IGrabber {
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

    constructor(name: string, hand: WebXRHand, sessionManager: WebXRSessionManager, scene: PhysicsEnabledScene) {
        super(name, scene);
        scene.grabbers.set(hand.xrController.inputSource.handedness, this);
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

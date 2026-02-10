import { IDisposable, AbstractMesh, Observer, Scene, Observable, WebXRHandJoint, Vector3 } from "@babylonjs/core";
import { PhysicsEnabledScene } from "./physicsEnabledScene";

enum PokeButtonState {
    Default,
    Poked,
    Cooldown,
}

export class PokeButton implements IDisposable {
    private scene: PhysicsEnabledScene;
    private onButton: AbstractMesh;
    private offButton: AbstractMesh | undefined;

    private beforeRenderObserver: Observer<Scene>;
    private state: PokeButtonState = PokeButtonState.Default;
    private cooldown: number = 0;

    private static readonly POKE_ENTER_THRESHOLD = 0.015;
    private static readonly POKE_EXIT_THRESHOLD = 0.02;

    public onPokeObservable: Observable<boolean> = new Observable<boolean>();

    public constructor(scene: PhysicsEnabledScene, onButtonName: string, offButtonName: string | undefined = undefined) {
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

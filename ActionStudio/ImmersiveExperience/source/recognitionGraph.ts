import { IDisposable, AbstractMesh, ParticleSystem, DynamicTexture, Observer, Scene, Color4, MeshBuilder, Vector3 } from "@babylonjs/core";
import { ICarlRecognizer } from "./carlInterfaces";
import { PhysicsEnabledScene } from "./physicsEnabledScene";

export class RecognitionGraph implements IDisposable {
    private _recognizer: ICarlRecognizer | undefined = undefined;
    private _emitter: AbstractMesh;
    private _particleSystem: ParticleSystem;
    private _texture: DynamicTexture;
    private _observer: Observer<Scene>;

    public get recognizer(): ICarlRecognizer | undefined {
        return this._recognizer;
    }

    public set recognizer(value: ICarlRecognizer | undefined) {
        this._recognizer?.dispose();
        if (this._recognizer && !value) {
            this._particleSystem.stop();
        } else if (!this._recognizer && value) {
            this._particleSystem.start();
        }
        this._recognizer = value;
    }

    public constructor(scene: PhysicsEnabledScene, color: Color4) {
        this._emitter = MeshBuilder.CreateBox("recognizerEmitter", {size: 0.01});
        this._emitter.position.set(1, 1, 2);
        this._emitter.isVisible = false;

        this._particleSystem = new ParticleSystem("particles", 2000, scene);

        // TODO: Pull this out so that all the graphs can use the same texture.
        const size = 1;
        this._texture = new DynamicTexture("whiteTexture", { width: size, height: size }, scene, false);
        const ctx = this._texture.getContext();
        ctx.fillStyle = "#FFFFFF";
        ctx.fillRect(0, 0, size, size);
        this._texture.update();
        this._particleSystem.particleTexture = this._texture;

        this._particleSystem.color1 = color;
        this._particleSystem.color2 = color;
        this._particleSystem.colorDead = color;
        this._particleSystem.blendMode = ParticleSystem.BLENDMODE_ONEONE;

        // Position where the particles are emitted from
        this._particleSystem.emitter = this._emitter;
        this._particleSystem.minEmitBox = Vector3.Zero();
        this._particleSystem.maxEmitBox = Vector3.Zero();
        
        this._particleSystem.minSize = 0.01;
        this._particleSystem.maxSize = 0.01;
        
        this._particleSystem.emitRate = 50;
        this._particleSystem.minLifeTime = 4;
        this._particleSystem.maxLifeTime = 4;
        
        this._particleSystem.minEmitPower = 0.5;
        this._particleSystem.maxEmitPower = 0.5;

        this._particleSystem.direction1 = Vector3.Left();
        this._particleSystem.direction2 = Vector3.Left();

        this._observer = scene.onBeforeRenderObservable.add(() => {
            if (this._recognizer) {
                this._emitter.position.y = 1 + 0.7 * this._recognizer.currentScore();
            }
        });
    }

    public dispose(): void {
        this._recognizer?.dispose();
        this._observer.remove();
        this._particleSystem.dispose();
        this._texture.dispose();
        this._emitter.dispose();
    }
}

/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * RecognitionGraph — visualises a CARL recognizer's current score as a particle stream.
 *
 * A particle emitter floats at a height proportional to `recognizer.currentScore()`, painting
 * a scrolling bar graph in 3D space.  Multiple RecognitionGraph instances share a single
 * white DynamicTexture (ref-counted) to avoid redundant GPU allocations.
 */
import { IDisposable, AbstractMesh, ParticleSystem, DynamicTexture, Observer, Scene, Color4, MeshBuilder, Vector3 } from "@babylonjs/core";
import { ICarlRecognizer } from "./carlInterfaces";
import { PhysicsEnabledScene } from "./physicsEnabledScene";

export class RecognitionGraph implements IDisposable {
    private static _sharedTexture: DynamicTexture | null = null;
    private static _sharedTextureRefCount: number = 0;

    private _recognizer: ICarlRecognizer | undefined = undefined;
    private _emitter: AbstractMesh;
    private _particleSystem: ParticleSystem;
    private _texture: DynamicTexture;
    private _onBeforeRenderObserver: Observer<Scene>;
    private _onSceneDisposeObserver: Observer<Scene>;

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
        this._emitter.position.set(1, 1, 0).addInPlace(Vector3.RightHandedForwardReadOnly).addInPlace(Vector3.RightHandedForwardReadOnly);
        this._emitter.isVisible = false;

        this._particleSystem = new ParticleSystem("particles", 2000, scene);

        if (RecognitionGraph._sharedTexture === null) {
            const size = 1;
            RecognitionGraph._sharedTexture = new DynamicTexture("whiteTexture", { width: size, height: size }, scene, false);
            const ctx = RecognitionGraph._sharedTexture.getContext();
            ctx.fillStyle = "#FFFFFF";
            ctx.fillRect(0, 0, size, size);
            RecognitionGraph._sharedTexture.update();
        }
        RecognitionGraph._sharedTextureRefCount++;
        this._texture = RecognitionGraph._sharedTexture;
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

        this._onBeforeRenderObserver = scene.onBeforeRenderObservable.add(() => {
            if (this._recognizer) {
                this._emitter.position.y = 1 + 0.7 * this._recognizer.currentScore();
            }
        });

        this._onSceneDisposeObserver = scene.onDisposeObservable.add(() => {
            this.dispose();
        });
    }

    public dispose(): void {
        this._recognizer = undefined;
        this._onBeforeRenderObserver.remove();
        this._onSceneDisposeObserver.remove();
        this._particleSystem.dispose();
        RecognitionGraph._sharedTextureRefCount--;
        if (RecognitionGraph._sharedTextureRefCount === 0) {
            RecognitionGraph._sharedTexture!.dispose();
            RecognitionGraph._sharedTexture = null;
        }
        this._emitter.dispose();
    }
}

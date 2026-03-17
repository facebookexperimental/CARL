/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * ExamplePreviewer — drives the InputPuppet to animate a stored example inside the XR scene.
 *
 * When an example block is dropped into the editor volume, `previewExample` runs a coroutine
 * that scrubs through the recording in sync with two in-scene slider meshes (trim start/end).
 * Calling `play()` triggers a real-time playback pass from minT to maxT.  On exit the updated
 * trim bounds are written back to the example via `setStartTimestamp` / `setEndTimestamp`.
 */
import { Vector3, Observer } from "@babylonjs/core";
import { ICarlExample } from "./carlInterfaces";
import { PhysicsEnabledScene } from "./physicsEnabledScene";
import { SliderBehavior } from "./slider";

export class ExamplePreviewer {
    private scene: PhysicsEnabledScene;
    private sliders: SliderBehavior[];
    private previewing: boolean = false;
    private shouldPlay: boolean = false;

    public constructor(scene: PhysicsEnabledScene) {
        this.scene = scene;
        this.sliders = [];

        scene.meshes.map(mesh => {
            if (mesh.name.startsWith("editor_slider_")) {
                this.sliders.push(SliderBehavior.GetForNode(mesh)!);
            }
        });
    }

    public previewExample(example: ICarlExample, onStopped?: (example: ICarlExample) => void): () => void {
        this.previewing = true;
        this.scene.onBeforeRenderObservable.runCoroutineAsync(this.previewExampleCoroutine(example, onStopped));
        return () => {
            this.previewing = false;
        }
    }

    public play(): void {
        this.shouldPlay = this.previewing;
    }

    private *previewExampleCoroutine(example: ICarlExample, onStopped?: (example: ICarlExample) => void) {
        const IMMITATION_ORIGIN = Vector3.RightHandedForwardReadOnly.scale(2);
        const inspector = example.getRecordingInspector();
        const duration = inspector.getEndTimestamp() - inspector.getStartTimestamp();
        this.scene.inputPuppet?.immitateInputSample(inspector.inspect(inspector.getStartTimestamp()), Vector3.ZeroReadOnly, Vector3.RightHandedForwardReadOnly, IMMITATION_ORIGIN, Vector3.RightHandedBackwardReadOnly);
        this.scene.inputPuppet?.setEnabled(true);

        let currentT: number | undefined = undefined;

        let minT = example.getStartTimestamp();
        let maxT = example.getEndTimestamp();

        this.sliders[0].value = 0;
        this.sliders[1].value = 1;

        const sliderObservers: Observer<void>[] = [];
        for (let idx = 0; idx < this.sliders.length; ++idx) {
            const slider = this.sliders[idx];
            const otherSlider = this.sliders[(idx + 1) % 2];
            sliderObservers.push(slider.onUpdatedObservable.add(() => {
                currentT = undefined;

                let t = slider.value;
                t = t * duration + inspector.getStartTimestamp();

                if (slider.value < otherSlider.value) {
                    minT = t;
                } else {
                    maxT = t;
                }

                const sample = inspector.inspect(t);
                this.scene.inputPuppet?.immitateInputSample(sample, Vector3.ZeroReadOnly, Vector3.RightHandedForwardReadOnly, IMMITATION_ORIGIN, Vector3.RightHandedBackwardReadOnly);
            }));
        }

        this.sliders[1].value = (maxT - inspector.getStartTimestamp()) / duration;
        this.sliders[0].value = (minT - inspector.getStartTimestamp()) / duration;
        
        while (this.previewing) {
            if (this.shouldPlay) {
                currentT = minT;
                this.shouldPlay = false;
            }

            if (currentT) {
                const sample = inspector.inspect(currentT);
                this.scene.inputPuppet?.immitateInputSample(sample, Vector3.ZeroReadOnly, Vector3.RightHandedForwardReadOnly, IMMITATION_ORIGIN, Vector3.RightHandedBackwardReadOnly);

                currentT += (this.scene.deltaTime / 1000);
                if (currentT > maxT) {
                    currentT = undefined;
                }
            }
            yield;
        }

        sliderObservers.forEach(observer => {
            observer.remove();
        });

        example.setStartTimestamp(minT);
        example.setEndTimestamp(maxT);
        inspector.dispose();
        onStopped?.(example);
        this.scene.inputPuppet?.setEnabled(false);
    }
}

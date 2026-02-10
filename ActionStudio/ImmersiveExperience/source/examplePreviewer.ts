import { Vector3, Observer } from "@babylonjs/core";
import { ICarlExample } from "./carlInterfaces";
import { PhysicsEnabledScene } from "./physicsEnabledScene";

export class ExamplePreviewer {
    private scene: PhysicsEnabledScene;
    private previewing: boolean = false;
    private shouldPlay: boolean = false;

    public constructor(scene: PhysicsEnabledScene) {
        this.scene = scene;
    }

    public previewExample(example: ICarlExample): () => void {
        this.previewing = true;
        this.scene.onBeforeRenderObservable.runCoroutineAsync(this.previewExampleCoroutine(example));
        return () => {
            this.previewing = false;
        }
    }

    public play(): void {
        this.shouldPlay = true;
    }

    private *previewExampleCoroutine(example: ICarlExample) {
        const IMMITATION_ORIGIN = Vector3.RightHandedBackwardReadOnly.scale(2);
        const inspector = example.getRecordingInspector();
        const duration = inspector.getEndTimestamp() - inspector.getStartTimestamp();
        this.scene.inputPuppet?.immitateInputSample(inspector.inspect(inspector.getStartTimestamp()), Vector3.ZeroReadOnly, Vector3.RightHandedBackwardReadOnly, IMMITATION_ORIGIN, Vector3.RightHandedForwardReadOnly);
        this.scene.inputPuppet?.setEnabled(true);

        let currentT: number | undefined = undefined;

        let minT = example.getStartTimestamp();
        let maxT = example.getEndTimestamp();

        this.scene.sliders[0].value = 0;
        this.scene.sliders[1].value = 1;

        const sliderObservers: Observer<void>[] = [];
        for (let idx = 0; idx < this.scene.sliders.length; ++idx) {
            const slider = this.scene.sliders[idx];
            const otherSlider = this.scene.sliders[(idx + 1) % 2];
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
                this.scene.inputPuppet?.immitateInputSample(sample, Vector3.ZeroReadOnly, Vector3.RightHandedBackwardReadOnly, IMMITATION_ORIGIN, Vector3.RightHandedForwardReadOnly);
            }));
        }

        this.scene.sliders[1].value = (maxT - inspector.getStartTimestamp()) / duration;
        this.scene.sliders[0].value = (minT - inspector.getStartTimestamp()) / duration;
        
        while (this.previewing) {
            if (this.shouldPlay) {
                currentT = minT;
                this.shouldPlay = false;
            }

            if (currentT) {
                const sample = inspector.inspect(currentT);
                this.scene.inputPuppet?.immitateInputSample(sample, Vector3.ZeroReadOnly, Vector3.RightHandedBackwardReadOnly, IMMITATION_ORIGIN, Vector3.RightHandedForwardReadOnly);

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
        this.scene.inputPuppet?.setEnabled(false);
    }
}

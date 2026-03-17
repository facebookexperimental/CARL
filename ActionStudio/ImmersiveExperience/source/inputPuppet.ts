/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * InputPuppet renders a ghost copy of recorded hand-joint poses alongside the live hands.
 *
 * It clones the live hand joint meshes at construction time and re-positions them each
 * frame by delegating to `applyJointSampleToMeshes` in `utils.ts`.  The puppet is used
 * by ExamplePreviewer to visualise a stored example while the user edits its trim bounds.
 */
import { AbstractMesh, Vector3, WebXRHand } from "@babylonjs/core";
import { ICarlInputSample } from "./carlInterfaces";
import { OPENXR_JOINT_MAPPINGS, applyJointSampleToMeshes } from "./utils";

export class InputPuppet {
    private _leftHandMeshes: AbstractMesh[];
    private _rightHandMeshes: AbstractMesh[];

    public constructor(leftHand: WebXRHand, rightHand: WebXRHand) {
        this._leftHandMeshes = OPENXR_JOINT_MAPPINGS.map(joint => {
            const mesh = leftHand.getJointMesh(joint);
            return mesh.clone(joint, null)!;
        });
        this._rightHandMeshes = OPENXR_JOINT_MAPPINGS.map(joint => {
            const mesh = rightHand.getJointMesh(joint);
            return mesh.clone(joint, null)!;
        });

        this.setEnabled(false);
    }

    public immitateInputSample(sample: ICarlInputSample, samplePosition: Vector3, sampleForward: Vector3, povPosition: Vector3, povForward: Vector3) {
        applyJointSampleToMeshes(sample, this._leftHandMeshes, this._rightHandMeshes,
            samplePosition, sampleForward, povPosition, povForward);
    }

    public setEnabled(enabled: boolean): void {
        this._leftHandMeshes.map(mesh => mesh.setEnabled(enabled));
        this._rightHandMeshes.map(mesh => mesh.setEnabled(enabled));
    }
}

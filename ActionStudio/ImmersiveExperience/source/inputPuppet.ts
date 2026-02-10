import { AbstractMesh, Vector3, Quaternion, Matrix, WebXRHand } from "@babylonjs/core";
import { ICarlInputSample } from "./carlInterfaces";
import { OPENXR_JOINT_MAPPINGS } from "./utils";

export class InputPuppet {
    private _leftHandMeshes: AbstractMesh[];
    private _rightHandMeshes: AbstractMesh[];
    private _scratchVec: Vector3 = new Vector3();
    private _scratchQuat: Quaternion = new Quaternion();
    private _povMat: Matrix = new Matrix();
    private _sampleMat: Matrix = new Matrix();
    private _sampleToPovMat: Matrix = new Matrix();
    private _jointMat: Matrix = new Matrix();

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
        povPosition.addToRef(povForward, this._scratchVec);
        Matrix.LookAtRHToRef(povPosition, this._scratchVec, Vector3.UpReadOnly, this._povMat);
        samplePosition.addToRef(sampleForward, this._scratchVec);
        Matrix.LookAtRHToRef(samplePosition, this._scratchVec, Vector3.UpReadOnly, this._sampleMat);

        this._povMat.invertToRef(this._povMat);
        this._sampleMat.multiplyToRef(this._povMat, this._sampleToPovMat);

        for (let idx = 0; idx < OPENXR_JOINT_MAPPINGS.length; ++idx) {
            let pose = sample.leftHandJointPoses[idx];
            this._scratchQuat.copyFromFloats(pose.orientation.x, pose.orientation.y, pose.orientation.z, pose.orientation.w);
            Matrix.FromQuaternionToRef(this._scratchQuat, this._jointMat);
            this._jointMat.setTranslationFromFloats(pose.position.x, pose.position.y, pose.position.z);
            this._jointMat.multiplyToRef(this._sampleToPovMat, this._jointMat);
            this._jointMat.decompose(undefined, this._leftHandMeshes[idx].rotationQuaternion!, this._leftHandMeshes[idx].position);
            
            pose = sample.rightHandJointPoses[idx];
            this._scratchQuat.copyFromFloats(pose.orientation.x, pose.orientation.y, pose.orientation.z, pose.orientation.w);
            Matrix.FromQuaternionToRef(this._scratchQuat, this._jointMat);
            this._jointMat.setTranslationFromFloats(pose.position.x, pose.position.y, pose.position.z);
            this._jointMat.multiplyToRef(this._sampleToPovMat, this._jointMat);
            this._jointMat.decompose(undefined, this._rightHandMeshes[idx].rotationQuaternion!, this._rightHandMeshes[idx].position);
        }
    }

    public setEnabled(enabled: boolean): void {
        this._leftHandMeshes.map(mesh => mesh.setEnabled(enabled));
        this._rightHandMeshes.map(mesh => mesh.setEnabled(enabled));
    }
}

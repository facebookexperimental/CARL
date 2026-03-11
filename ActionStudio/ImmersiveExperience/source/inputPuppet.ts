import { AbstractMesh, Vector3, Quaternion, Matrix, WebXRHand } from "@babylonjs/core";
import { ICarlInputSample } from "./carlInterfaces";
import { convertOpenXRSampleToBabylonConventions, OPENXR_JOINT_MAPPINGS } from "./utils";

export class InputPuppet {
    private _leftHandMeshes: AbstractMesh[];
    private _rightHandMeshes: AbstractMesh[];
    private _scratchVec: Vector3 = new Vector3();
    private _scratchQuat: Quaternion = new Quaternion();
    private _povMat: Matrix = new Matrix();
    private _sampleMat: Matrix = new Matrix();
    private _sampleToPovMat: Matrix = new Matrix();
    private _jointMat: Matrix = new Matrix();
    private _convertedSample: ICarlInputSample | undefined;

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
        if (!this._convertedSample) {
            this._convertedSample = JSON.parse(JSON.stringify(sample));
        } else {
            this._deepCopySample(sample, this._convertedSample);
        }
        convertOpenXRSampleToBabylonConventions(this._convertedSample!);

        povPosition.addToRef(povForward, this._scratchVec);
        Matrix.LookAtLHToRef(povPosition, this._scratchVec, Vector3.UpReadOnly, this._povMat);
        samplePosition.addToRef(sampleForward, this._scratchVec);
        Matrix.LookAtLHToRef(samplePosition, this._scratchVec, Vector3.UpReadOnly, this._sampleMat);

        this._povMat.invertToRef(this._povMat);
        this._sampleMat.multiplyToRef(this._povMat, this._sampleToPovMat);

        for (let idx = 0; idx < OPENXR_JOINT_MAPPINGS.length; ++idx) {
            let pose = this._convertedSample!.leftHandJointPoses[idx];
            this._scratchQuat.copyFromFloats(pose.orientation.x, pose.orientation.y, pose.orientation.z, pose.orientation.w);
            Matrix.FromQuaternionToRef(this._scratchQuat, this._jointMat);
            this._jointMat.setTranslationFromFloats(pose.position.x, pose.position.y, pose.position.z);
            this._jointMat.multiplyToRef(this._sampleToPovMat, this._jointMat);
            this._jointMat.decompose(undefined, this._leftHandMeshes[idx].rotationQuaternion!, this._leftHandMeshes[idx].position);

            pose = this._convertedSample!.rightHandJointPoses[idx];
            this._scratchQuat.copyFromFloats(pose.orientation.x, pose.orientation.y, pose.orientation.z, pose.orientation.w);
            Matrix.FromQuaternionToRef(this._scratchQuat, this._jointMat);
            this._jointMat.setTranslationFromFloats(pose.position.x, pose.position.y, pose.position.z);
            this._jointMat.multiplyToRef(this._sampleToPovMat, this._jointMat);
            this._jointMat.decompose(undefined, this._rightHandMeshes[idx].rotationQuaternion!, this._rightHandMeshes[idx].position);
        }
    }

    private _deepCopySample(src: ICarlInputSample, dst: ICarlInputSample): void {
        dst.timestamp = src.timestamp;
        this._copyPose(src.hmdPose, dst.hmdPose);
        this._copyPose(src.leftWristPose, dst.leftWristPose);
        this._copyPose(src.rightWristPose, dst.rightWristPose);
        for (let i = 0; i < src.leftHandJointPoses.length; i++) {
            this._copyPose(src.leftHandJointPoses[i], dst.leftHandJointPoses[i]);
        }
        for (let i = 0; i < src.rightHandJointPoses.length; i++) {
            this._copyPose(src.rightHandJointPoses[i], dst.rightHandJointPoses[i]);
        }
    }

    private _copyPose(src: ICarlInputSample["hmdPose"], dst: ICarlInputSample["hmdPose"]): void {
        dst.valid = src.valid;
        dst.position.x = src.position.x;
        dst.position.y = src.position.y;
        dst.position.z = src.position.z;
        dst.orientation.w = src.orientation.w;
        dst.orientation.x = src.orientation.x;
        dst.orientation.y = src.orientation.y;
        dst.orientation.z = src.orientation.z;
    }

    public setEnabled(enabled: boolean): void {
        this._leftHandMeshes.map(mesh => mesh.setEnabled(enabled));
        this._rightHandMeshes.map(mesh => mesh.setEnabled(enabled));
    }
}

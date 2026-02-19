import { AbstractMesh, TransformNode, PhysicsAggregate, PhysicsShapeType, Quaternion, IDisposable } from "@babylonjs/core";
import { PhysicsEnabledScene } from "./physicsEnabledScene";
import { PhysicsGrabBehavior } from "./physicsGrabBehavior";

export class BlockSpawner<T> {
    private _scene: PhysicsEnabledScene;
    private _nextId: number;
    private _block: AbstractMesh;
    private _nodeToT = new Map<TransformNode, T>();

    constructor (scene: PhysicsEnabledScene, name: string) {
        this._scene = scene;
        this._nextId = 0;
        this._block = scene.getMeshByName(name)!;
        this._block.setEnabled(false);
    }

    // Note: blocks do NOT take ownership of the value they represent and WILL NOT dispose it on their dispose.
    // This must be handled separately.
    spawnNewBlock(value: T): AbstractMesh {
            const clonedExampleBlock = this._block.clone(`${this._block.name}_${++this._nextId}`, null)!;
            const aggregate = new PhysicsAggregate(clonedExampleBlock, PhysicsShapeType.BOX, { mass: 1 }, this._scene);
            aggregate.transformNode.rotationQuaternion = Quaternion.FromEulerAngles(
                aggregate.transformNode.rotation.x, 
                aggregate.transformNode.rotation.y, 
                aggregate.transformNode.rotation.z);
            this._scene.grabbables.push(PhysicsGrabBehavior.get(clonedExampleBlock));
            clonedExampleBlock.setEnabled(true);
            
            this._nodeToT.set(clonedExampleBlock, value);
            clonedExampleBlock.onDisposeObservable.addOnce(() => {
                this._nodeToT.delete(clonedExampleBlock);
            });

            return clonedExampleBlock;
    }

    getValueFromBlock(block: TransformNode): T | undefined {
        return this._nodeToT.get(block);
    }
}

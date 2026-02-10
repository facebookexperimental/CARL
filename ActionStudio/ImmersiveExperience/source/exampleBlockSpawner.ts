import { AbstractMesh, TransformNode, PhysicsAggregate, PhysicsShapeType, Quaternion } from "@babylonjs/core";
import { ICarlExample } from "./carlInterfaces";
import { PhysicsEnabledScene } from "./physicsEnabledScene";
import { PhysicsGrabBehavior } from "./physicsGrabBehavior";

export class ExampleBlockSpawner {
    private _scene: PhysicsEnabledScene;
    private _nextId: number;
    private _exampleBlock: AbstractMesh;
    private _nodeToExample = new Map<TransformNode, ICarlExample>();

    constructor (scene: PhysicsEnabledScene) {
        this._scene = scene;
        this._nextId = 0;
        this._exampleBlock = scene.getMeshByName("example_block")!;
        this._exampleBlock.setEnabled(false);
    }

    spawnNewExampleBlock(example: ICarlExample): AbstractMesh {
            const clonedExampleBlock = this._exampleBlock.clone(`example_block_${++this._nextId}`, null)!;
            const aggregate = new PhysicsAggregate(clonedExampleBlock, PhysicsShapeType.BOX, { mass: 1 }, this._scene);
            aggregate.transformNode.rotationQuaternion = Quaternion.FromEulerAngles(
                aggregate.transformNode.rotation.x, 
                aggregate.transformNode.rotation.y, 
                aggregate.transformNode.rotation.z);
            this._scene.grabbables.push(PhysicsGrabBehavior.get(clonedExampleBlock));
            clonedExampleBlock.setEnabled(true);
            
            this._nodeToExample.set(clonedExampleBlock, example);
            clonedExampleBlock.onDisposeObservable.addOnce(() => {
                this._nodeToExample.delete(clonedExampleBlock);
            });

            return clonedExampleBlock;
    }

    getExampleFromBlock(block: TransformNode): ICarlExample | undefined {
        return this._nodeToExample.get(block);
    }
}

import { PhysicsAggregate, PhysicsBody, TransformNode, Vector3 } from "@babylonjs/core";
import { PhysicsGrabBehavior } from "./physicsGrabBehavior";

export class SliderConstraints {
    public a: Vector3;
    public b: Vector3;

    private _ab = new Vector3();
    private _ap = new Vector3();

    constructor(a: Vector3, b: Vector3) {
        this.a = a.clone();
        this.b = b.clone();
    }

    public handleTransformMoved(p: Vector3): number {
        this.b.subtractToRef(this.a, this._ab); // Vector from A to B
        p.subtractToRef(this.a, this._ap); // Vector from A to P
        const ab2 = Vector3.Dot(this._ab, this._ab); // |AB|^2
        if (ab2 === 0) {
            p.copyFrom(this.a); // A and B are the same point
        }

        // Projection factor t of AP onto AB
        let t = Vector3.Dot(this._ap, this._ab) / ab2;

        // Clamp t to [0, 1] for segment; remove clamp for infinite line
        t = Math.max(0, Math.min(1, t));

        // Nearest point = A + t * AB
        this._ab.scaleInPlace(t);
        this.a.addToRef(this._ab, p);

        return t;
    }
}

export class SliderBehavior {
    private static _idToValue = new Map<number, number>();

    static attach(node: TransformNode, grabBehavior: PhysicsGrabBehavior, constraints: SliderConstraints): void {
        node.physicsBody!.setMassProperties({
            mass: node.physicsBody!.getMassProperties().mass,
            inertia: Vector3.Zero(),
        });
        const initialOrientation = node.rotationQuaternion!.clone();
        grabBehavior.onGrabMovedObservable.add(() => {
            node.rotationQuaternion?.copyFrom(initialOrientation);
            const t = constraints.handleTransformMoved(node.position);
            SliderBehavior._idToValue.set(node.uniqueId, t);
        });

        SliderBehavior._idToValue.set(node.uniqueId, 0);
        node.onDisposeObservable.add(() => {
            SliderBehavior._idToValue.delete(node.uniqueId);
        });
    }

    static GetValue(node: TransformNode): number {
        return SliderBehavior._idToValue.get(node.uniqueId) ?? -1;
    }
}

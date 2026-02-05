import { Observable, Observer, TransformNode, Vector3 } from "@babylonjs/core";
import { PhysicsGrabBehavior } from "./physicsGrabBehavior";

export class SliderConstraints {
    public a: Vector3;
    public b: Vector3;
    constructor(a: Vector3, b: Vector3) {
        this.a = a.clone();
        this.b = b.clone();
    }
}

export class SliderBehavior {
    private _node: TransformNode;
    private _constraints: SliderConstraints;
    private _value: number = -1;
    private _onGrabMoveObserver: Observer<TransformNode>;

    private _ab = new Vector3();
    private _ap = new Vector3();

    public get value() {
        return this._value;
    }
    public set value(newValue: number) {
        this._value = newValue;
        this._constraints.b.subtractToRef(this._constraints.a, this._ab);
        this._ab.scaleInPlace(this._value);
        this._constraints.a.addToRef(this._ab, this._node.position);
        this.onUpdatedObservable.notifyObservers();
    }
    public onUpdatedObservable: Observable<void> = new Observable<void>();

    private static _idToSlider = new Map<number, SliderBehavior>();

    private constructor(node: TransformNode, grabBehavior: PhysicsGrabBehavior, constraints: SliderConstraints, value: number = 0) {
        this._node = node;
        this._constraints = constraints;
        this.value = value;

        node.physicsBody!.setMassProperties({
            mass: node.physicsBody!.getMassProperties().mass,
            inertia: Vector3.Zero(),
        });

        const initialOrientation = this._node.rotationQuaternion!.clone();
        this._onGrabMoveObserver = grabBehavior.onGrabMovedObservable.add(() => {
            const t = this._calculateTAfterMove(this._node.position);
            this._node.rotationQuaternion?.copyFrom(initialOrientation);
            this.value = t;
        });

        this._node.onDisposeObservable.addOnce(() => {
            this.dispose();
        })

        SliderBehavior._idToSlider.set(node.uniqueId, this);
    }

    public dispose(): void {
        this._onGrabMoveObserver.remove();
        SliderBehavior._idToSlider.delete(this._node.uniqueId);
    }

    static attach(node: TransformNode, grabBehavior: PhysicsGrabBehavior, constraints: SliderConstraints): SliderBehavior {
        const initialValue = SliderBehavior._idToSlider.size > 0 ? 1 : 0;
        return new SliderBehavior(node, grabBehavior, constraints, initialValue);
    }

    static GetForNode(node: TransformNode): SliderBehavior | undefined {
        return SliderBehavior._idToSlider.get(node.uniqueId);
    }

    public _calculateTAfterMove(p: Vector3): number {
        this._constraints.b.subtractToRef(this._constraints.a, this._ab);
        p.subtractToRef(this._constraints.a, this._ap);
        const ab2 = Vector3.Dot(this._ab, this._ab);
        if (ab2 === 0) {
            p.copyFrom(this._constraints.a);
        }
        let t = Vector3.Dot(this._ap, this._ab) / ab2;
        return Math.max(0, Math.min(1, t));
    }
}

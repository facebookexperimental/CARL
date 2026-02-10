import { AbstractMesh, AppendSceneAsync, Engine, HavokPlugin, PhysicsAggregate, PhysicsShapeType, Quaternion, Scene, TransformNode, Vector3, WebXRHand } from "@babylonjs/core";
import { PhysicsGrabBehavior } from "./physicsGrabBehavior";
import { HandPinchGrabber } from "./handPinchGrabber";
import { SliderBehavior, SliderConstraints } from "./slider";
import { InputPuppet } from "./main";

export class PhysicsEnabledScene extends Scene {
    public grabbables: PhysicsGrabBehavior[] = [];
    public grabbers = new Map<string, HandPinchGrabber>();
    public sliders: SliderBehavior[] = [];

    public leftHand: WebXRHand | undefined;
    public rightHand: WebXRHand | undefined;

    public inputPuppet: InputPuppet | undefined;

    private constructor(engine: Engine) {
        super(engine);
    }

    public static async loadAsync(url: string, engine: Engine, havok: HavokPlugin): Promise<PhysicsEnabledScene> {
        const scene = new PhysicsEnabledScene(engine);
        scene.useRightHandedSystem = true;
        await AppendSceneAsync(url, scene);
        scene.enablePhysics(Vector3.DownReadOnly.scale(9.81), havok);

        const physicsNulls: TransformNode[] =  [];

        scene.transformNodes.forEach(node => {
            if (node.name.startsWith("physics")) {
                physicsNulls.push(node);
            }
        });

        physicsNulls.forEach(physicsNull => {
            const mesh = physicsNull.parent! as AbstractMesh;
            const shapeType = physicsNull.name.split(';')[0];
            let shape: PhysicsShapeType;
            switch (shapeType) {
                case "physics_sphere":
                    shape = PhysicsShapeType.SPHERE;
                    break;
                case "physics_mesh":
                    shape = PhysicsShapeType.MESH;
                    break;
                default:
                    shape = PhysicsShapeType.BOX;
                    break;
            }
            physicsNull.setParent(null);
            physicsNull.dispose();

            let mass = 0;
            const massNulls = mesh.getChildTransformNodes(true, child => child.name.startsWith("mass:"));
            if (massNulls.length > 0) {
                const massNull = massNulls[0];
                mass = parseFloat(massNull.name.split(';')[0].split(':')[1]);
                massNull.setParent(null);
                massNull.dispose();
            }
            
            let grabbable = false;
            let sliderConstraints: SliderConstraints | undefined = undefined;
            const grabbableNulls = mesh.getChildTransformNodes(true, child => child.name.startsWith("grabbable"));
            if (grabbableNulls.length > 0) {
                grabbable = true;
                const grabbableNull = grabbableNulls[0];
                grabbableNull.setParent(null);
                grabbableNull.dispose();

                const sliderNulls = mesh.getChildTransformNodes(true, child => child.name.startsWith("slider"));
                if (sliderNulls.length === 2) {
                    const ab = sliderNulls[0].name.split(';')[1].startsWith("a");
                    const a = ab ? sliderNulls[0] : sliderNulls[1];
                    const b = ab ? sliderNulls[1] : sliderNulls[0];
                    a.setParent(null);
                    b.setParent(null);
                    sliderConstraints = new SliderConstraints(a.position, b.position);
                    a.dispose();
                    b.dispose();
                }
            }

            const aggregate = new PhysicsAggregate(mesh, shape, { mass: mass }, scene);
            aggregate.transformNode.rotationQuaternion = Quaternion.FromEulerAngles(
                aggregate.transformNode.rotation.x, 
                aggregate.transformNode.rotation.y, 
                aggregate.transformNode.rotation.z);
            if (grabbable) {
                const behavior = PhysicsGrabBehavior.get(mesh);
                scene.grabbables.push(behavior);

                if (sliderConstraints) {
                    const slider = SliderBehavior.attach(aggregate.transformNode, behavior, sliderConstraints);
                    scene.sliders.push(slider);
                }
            }
        });

        return scene;
    }
}

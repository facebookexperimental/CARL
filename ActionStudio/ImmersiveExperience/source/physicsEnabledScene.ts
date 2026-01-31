import { AbstractMesh, AppendSceneAsync, Engine, HavokPlugin, PhysicsAggregate, PhysicsShapeType, Quaternion, Scene, TransformNode, Vector3 } from "@babylonjs/core";
import { PhysicsGrabBehavior } from "./physicsGrabBehavior";
import { HandPinchGrabber } from "./handPinchGrabber";

export class PhysicsEnabledScene extends Scene {
    public grabbables: PhysicsGrabBehavior[] = [];
    public grabbers = new Map<string, HandPinchGrabber>();

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
            const grabbableNulls = mesh.getChildTransformNodes(true, child => child.name.startsWith("grabbable"));
            if (grabbableNulls.length > 0) {
                grabbable = true;
                const grabbableNull = grabbableNulls[0];
                grabbableNull.setParent(null);
                grabbableNull.dispose();
            }

            const aggregate = new PhysicsAggregate(mesh, shape, { mass: mass }, scene);
            aggregate.transformNode.rotationQuaternion = Quaternion.Identity();
            if (grabbable) {
                scene.grabbables.push(PhysicsGrabBehavior.attach(mesh));
            }
        });

        return scene;
    }
}

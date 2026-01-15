import { Engine, FreeCamera, HavokPlugin, HemisphericLight, MeshBuilder, PhysicsAggregate, PhysicsShapeType, Scene, Vector3 } from "@babylonjs/core";
import HavokPhysics from "@babylonjs/havok";

export function initializeImmersiveExperience(canvas: HTMLCanvasElement): void {
    const engine = new Engine(canvas, true);

    const resizeHandler = () => {
        engine.resize();
    };
    window.addEventListener("resize", resizeHandler);
    engine.onDisposeObservable.add(() => {
        window.removeEventListener("resize", resizeHandler);
    });

    const scene = new Scene(engine);
    var camera = new FreeCamera("camera1", new Vector3(0, 5, -10), scene);

    camera.setTarget(Vector3.Zero());
    camera.attachControl(canvas, true);

    var light = new HemisphericLight("light", new Vector3(0, 1, 0), scene);
    light.intensity = 0.7;

    var sphere = MeshBuilder.CreateSphere("sphere", {diameter: 2, segments: 32}, scene);
    sphere.position.y = 4;
    var ground = MeshBuilder.CreateGround("ground", {width: 10, height: 10}, scene);

    HavokPhysics().then((havok) => {
        var hk = new HavokPlugin(true, havok);
        scene.enablePhysics(new Vector3(0, -9.8, 0), hk);
        var sphereAggregate = new PhysicsAggregate(sphere, PhysicsShapeType.SPHERE, { mass: 1, restitution:0.75}, scene);
        var groundAggregate = new PhysicsAggregate(ground, PhysicsShapeType.BOX, { mass: 0 }, scene);
    });

    engine.runRenderLoop(() => {
        scene.render();
    });
}

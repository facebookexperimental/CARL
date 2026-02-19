import * as React from 'react';
import * as ReactDOM from 'react-dom/client';
import { initializeImmersiveExperienceAsync } from 'carl-actionstudio-immersiveexperience';
import { initializeNativeIntegrationAsync } from 'carl-actionstudio-nativeintegration';

const downloadSerializedBytes = (jsBytes, filename) => {
    const buffer = jsBytes.buffer;
    const blob = new Blob([buffer], { type: "application/octet-stream" });
    const link = document.createElement('a');
    link.href = URL.createObjectURL(blob);
    link.download = filename;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    URL.revokeObjectURL(link.href);
};

class CarlRecordingInspector {
    _nativeInspector;

    constructor(nativeInspector) {
        this._nativeInspector= nativeInspector;
    }

    getStartTimestamp() {
        return this._nativeInspector.startTimestamp();
    }

    getEndTimestamp() {
        return this._nativeInspector.endTimestamp();
    }

    inspect(timestamp) {
        return this._nativeInspector.inspect(timestamp);
    }

    dispose() {
        this._nativeInspector.delete();
    }
}

class CarlExample {
    _nativeExample;

    constructor(nativeExample) {
        this._nativeExample = nativeExample;
    }

    getStartTimestamp() {
        return this._nativeExample.getStartTimestamp();
    }

    setStartTimestamp(timestamp) {
        this._nativeExample.setStartTimestamp(timestamp);
    }

    getEndTimestamp() {
        return this._nativeExample.getEndTimestamp();
    }

    setEndTimestamp(timestamp) {
        this._nativeExample.setEndTimestamp(timestamp);
    }

    getRecordingInspector() {
        const inspector = this._nativeExample.getRecordingInspector();
        return new CarlRecordingInspector(inspector);
    }

    dispose() {
        this._nativeExample.delete();
    }
}

class CarlDefinition {
    _nativeDefinition;

    constructor(nativeDefinition) {
        this._nativeDefinition = nativeDefinition;
    }

    download() {
        const bytes = this.serialize();
        downloadSerializedBytes(bytes, "definition_" + (Date.now() / 1000) + ".bin");
        bytes.delete();
    }

    serialize() {
        const bytes = this._nativeDefinition.serialize();
        const jsBytes = new Uint8Array(bytes.size());
        for (let idx = 0; idx < jsBytes.length; ++idx) {
            jsBytes[idx] = bytes.get(idx);
        }
        bytes.delete();
        return jsBytes;
    }

    getDefaultSensitivity() {
        return this._nativeDefinition.getDefaultSensitivity();
    }

    setDefaultSensitivity(sensitivity) {
        this._nativeDefinition.setDefaultSensitivity(sensitivity);
    }

    dispose() {
        this._nativeDefinition.delete();
    }
}

class CarlRecognizer {
    _nativeRecognizer;

    constructor(nativeRecognizer) {
        this._nativeRecognizer = nativeRecognizer;
    }

    currentScore() {
        return this._nativeRecognizer.currentScore();
    }

    getSensitivity() {
        return this._nativeRecognizer.getSensitivity();
    };

    setSensitivity(sensitivity) {
        this._nativeRecognizer.setSensitivity(sensitivity);
    }

    dispose() {
        this._nativeRecognizer.delete();
    }
}

class CarlIntegration {
    _carl;
    _session;
    _inProgressRecordings = new Map();

    _nextId = 1;

    _actionTypes = [];
    _idToActionType = new Map();

    constructor(carl) {
        this._carl = carl;
        this._session = new this._carl.Session();
        
        // TODO: Find a better way to deal with this block than this manual approach.
        let actionTypeId = 0;
        this._idToActionType.set(actionTypeId, this._carl.ACTION_TYPE.LeftHandPose);
        this._actionTypes.push({ name: "Left Hand Pose", typeId: actionTypeId });
        ++actionTypeId;
        this._idToActionType.set(actionTypeId, this._carl.ACTION_TYPE.LeftHandGesture);
        this._actionTypes.push({ name: "Left Hand Gesture", typeId: actionTypeId });
        ++actionTypeId;
        this._idToActionType.set(actionTypeId, this._carl.ACTION_TYPE.RightHandPose);
        this._actionTypes.push({ name: "Right Hand Pose", typeId: actionTypeId });
        ++actionTypeId;
        this._idToActionType.set(actionTypeId, this._carl.ACTION_TYPE.RightHandGesture);
        this._actionTypes.push({ name: "Right Hand Gesture", typeId: actionTypeId });
        ++actionTypeId;
        this._idToActionType.set(actionTypeId, this._carl.ACTION_TYPE.TwoHandGesture);
        this._actionTypes.push({ name: "Two Hand Gesture", typeId: actionTypeId });
        ++actionTypeId;
        this._idToActionType.set(actionTypeId, this._carl.ACTION_TYPE.LeftControllerGesture);
        this._actionTypes.push({ name: "Left Controller Gesture", typeId: actionTypeId });
        ++actionTypeId;
        this._idToActionType.set(actionTypeId, this._carl.ACTION_TYPE.RightControllerGesture);
        this._actionTypes.push({ name: "Right Controller Gesture", typeId: actionTypeId });
        ++actionTypeId;
        this._idToActionType.set(actionTypeId, this._carl.ACTION_TYPE.TwoControllerGesture);
        this._actionTypes.push({ name: "Two Controller Gesture", typeId: actionTypeId });
        ++actionTypeId;
        this._idToActionType.set(actionTypeId, this._carl.ACTION_TYPE.LeftWristTrajectory);
        this._actionTypes.push({ name: "Left Wrist Trajectory", typeId: actionTypeId });
        ++actionTypeId;
        this._idToActionType.set(actionTypeId, this._carl.ACTION_TYPE.RightWristTrajectory);
        this._actionTypes.push({ name: "Right Wrist Trajectory", typeId: actionTypeId });
        ++actionTypeId;
        this._idToActionType.set(actionTypeId, this._carl.ACTION_TYPE.LeftHandShape);
        this._actionTypes.push({ name: "Left Hand Shape", typeId: actionTypeId });
        ++actionTypeId;
        this._idToActionType.set(actionTypeId, this._carl.ACTION_TYPE.RightHandShape);
        this._actionTypes.push({ name: "Right Hand Shape", typeId: actionTypeId });
    }

    static async CreateAsync() {
        return new CarlIntegration(await initializeNativeIntegrationAsync());
    }

    startRecording() {
        const id = this._nextId;
        this._nextId += 1;
        this._inProgressRecordings.set(id, new this._carl.InProgressRecording());
        return id;
    }

    stopRecording(id) {
        const ipr = this._inProgressRecordings.get(id);
        this._inProgressRecordings.delete(id);

        const recording = new this._carl.Recording(ipr);
        const inspector = new this._carl.RecordingInspector(recording);
        const startT = inspector.startTimestamp();
        const endT = inspector.endTimestamp();

        const example = new this._carl.Example(recording, startT, endT);

        inspector.delete();
        recording.delete();
        ipr.delete();

        return new CarlExample(example)
    }

    getActionTypesMap() {
        return this._actionTypes;
    }

    createDefinition(actionTypeId, examples, counterexamples) {
        const definition = new this._carl.Definition(this._idToActionType.get(actionTypeId));
        examples.forEach(example => {
            definition.addExample(example._nativeExample);
        });
        counterexamples.forEach(example => {
            definition.addCounterexample(example._nativeExample);
        });
        return new CarlDefinition(definition);
    }

    createRecognizer(definition) {
        return new CarlRecognizer(new this._carl.Recognizer(this._session, definition._nativeDefinition));
    }

    createInputSample() {
        return this._carl.createInputSample();
    }

    handleInputSample(sample) {
        this._session.addInput(sample);
        this._inProgressRecordings.forEach(ipr => {
            ipr.addInput(sample);
        });
    }
}

class SerializationsDB {
    static _DB_NAME = "carl_files";
    static _OBJECT_STORE_NAME = "serializations";
    _db;

    static async loadAsync() {
        const db = new SerializationsDB();
        db._db = await new Promise((resolve, reject) => {
            const request = window.indexedDB.open(SerializationsDB._DB_NAME, 3);
            request.onsuccess = evt => resolve(evt.target.result);
            request.onerror = () => {
                reject("Error fetching database: recordings will not be persisted in-browser");
            };
            request.onupgradeneeded = (evt) => {
                const db = evt.target.result;

                db.createObjectStore(SerializationsDB._OBJECT_STORE_NAME, { keyPath: "id", autoIncrement: true });
            };
        });

        db._db.onerror = (evt) => {
            console.error(`Database error: ${evt.target.error?.message}`);
        };

        return db;
    }

    static async resetAsync() {
        await new Promise((resolve, reject) => {
            const request = window.indexedDB.deleteDatabase(SerializationsDB._DB_NAME);

            request.onerror = (evt) => {
                reject("Error deleting database: we are in an unknown state");
            };
            request.onsuccess = (evt) => {
                resolve();
            };
        });
    }

    async addAsync(carlObject, type) {
        const transaction = this._db.transaction(SerializationsDB._OBJECT_STORE_NAME, "readwrite")
        const os = transaction.objectStore(SerializationsDB._OBJECT_STORE_NAME);

        os.add({
            type: type,
            bytes: carlObject.serialize(),
        });
        
        return await new Promise((resolve, reject) => {
            transaction.onerror = () => reject("Error storing CARL object");
            transaction.onsuccess = resolve();
        });
    }

    async fetchAllAsync() {
        const transaction = this._db.transaction(SerializationsDB._OBJECT_STORE_NAME, "readonly")
        const os = transaction.objectStore(SerializationsDB._OBJECT_STORE_NAME);

        const carlObjects = await new Promise((resolve, reject) => {
            const request = os.getAll();
            request.onerror = () => reject("Error in getAll() request");
            request.onsuccess = (evt) => {
                resolve(evt.target.result);
            };
        });
        
        return await new Promise((resolve, reject) => {
            transaction.onerror = () => reject("Error loading CARL objects");
            transaction.onsuccess = resolve(carlObjects);
        });
    }

    async updateAsync(id, carlObject) {
        const transaction = this._db.transaction(SerializationsDB._OBJECT_STORE_NAME, "readwrite")
        const os = transaction.objectStore(SerializationsDB._OBJECT_STORE_NAME);

        const fetched = await new Promise((resolve, reject) => {
            const request = os.get(id);
            request.onerror = () => reject("Error in get() request");
            request.onsuccess = (evt) => {
                resolve(evt.target.result);
            };
        });
        fetched.bytes = carlObject.serialize();

        await new Promise((resolve, reject) => {
            const request = os.put(fetched);
            request.onerror = () => reject("Error in put() request");
            request.onsuccess = () => resolve();
        });
        
        return await new Promise((resolve, reject) => {
            transaction.onerror = () => reject("Error loading CARL objects");
            transaction.onsuccess = resolve();
        });
    }

    async deleteAsync(id) {
        const transaction = this._db.transaction(SerializationsDB._OBJECT_STORE_NAME, "readwrite")
        const os = transaction.objectStore(SerializationsDB._OBJECT_STORE_NAME);
        
        os.delete(id);
        
        return await new Promise((resolve, reject) => {
            transaction.onerror = () => reject("Error deleting CARL objects");
            transaction.onsuccess = resolve();
        });
    }
}

const Greet = () => {
    const canvasRef = React.useRef(null);
    const buttonRef = React.useRef(null);
    async function bootstrapAsync() {
        const carl = await CarlIntegration.CreateAsync();
        const immersiveExperience = await initializeImmersiveExperienceAsync(canvasRef.current, carl);
        buttonRef.current.addEventListener("click", () => {
            immersiveExperience.enterImmersiveMode();
        });
    }
    React.useEffect(() => {
        bootstrapAsync();
    }, []);

    return <>
        <h1>Hello, world!</h1>
        <p>How's life?</p>
        <button ref={buttonRef}>Enter XR</button>
        <canvas ref={canvasRef} width={16} height={16} style={{ display: 'none' }}></canvas>
    </>;
};

const root = ReactDOM.createRoot(document.getElementById('root'));
root.render(<Greet />);

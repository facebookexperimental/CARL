import React, { useState } from 'react';
import { BrowserRouter as Router, Routes, Route, Navigate } from 'react-router-dom';
import Navigation from './components/Navigation';
import LibraryView from './components/LibraryView';
import PreviewMode from './components/PreviewMode';
import DefinitionBuilder from './components/DefinitionBuilder';
import { mockExamples, mockDefinitions } from './data/mockData';
import { initializeImmersiveExperienceAsync } from 'carl-actionstudio-immersiveexperience';
import { initializeNativeIntegrationAsync } from 'carl-actionstudio-nativeintegration';
import './styles/App.css';

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
    onDisposed = undefined;

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
        if (this.onDisposed) this.onDisposed();
        this._nativeExample.delete();
    }
}

class CarlDefinition {
    _nativeDefinition;
    onDisposed = undefined;

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
        if (this.onDisposed) this.onDisposed();
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

    onExampleCreated = undefined;
    onDefinitionCreated = undefined;

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

        const jsExample = new CarlExample(example);
        if (this.onExampleCreated) this.onExampleCreated(jsExample);
        return jsExample;
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
        const jsDefinition = new CarlDefinition(definition);
        if (this.onDefinitionCreated) this.onDefinitionCreated(jsDefinition);
        return jsDefinition;
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

function App() {
  const canvasRef = React.useRef(null);
  const [examples, setExamples] = useState(mockExamples);
  const [definitions, setDefinitions] = useState(mockDefinitions);
  let immersiveExperience = undefined;
  
  async function bootstrapAsync() {
    const carl = await CarlIntegration.CreateAsync();
    carl.onExampleCreated = example => {
      // examples.add(example);
      // setExamplesCount(examples.size);
      // example.onDisposed = () => {
      //   examples.delete(example);
      //   setExamplesCount(examples.size);
      // };
    };
    carl.onDefinitionCreated = definition => {
      // definitions.add(definition);
      // setDefinitionsCount(definitions.size);
      // definition.onDisposed = () => {
      //   definitions.delete(definition);
      //   setDefinitionsCount(definitions.size);
      // };
    };
    immersiveExperience = await initializeImmersiveExperienceAsync(canvasRef.current, carl);
  }
  React.useEffect(() => {
    bootstrapAsync();
  }, []);

  // Handler for recording new actions (stub for integration)
  const handleRecordNewActions = () => {
    immersiveExperience?.enterImmersiveMode();
  };

  // Handler for updating an example
  const updateExample = (id, updates) => {
    setExamples(examples.map(ex => ex.id === id ? { ...ex, ...updates } : ex));
  };

  // Handler for deleting an example
  const deleteExample = (id) => {
    setExamples(examples.filter(ex => ex.id !== id));
  };

  // Handler for updating a definition
  const updateDefinition = (id, updates) => {
    setDefinitions(definitions.map(def => def.id === id ? { ...def, ...updates } : def));
  };

  // Handler for deleting a definition
  const deleteDefinition = (id) => {
    setDefinitions(definitions.filter(def => def.id !== id));
  };

  // Handler for unpacking a definition
  const unpackDefinition = (id) => {
    const definition = definitions.find(def => def.id === id);
    if (definition) {
      console.log(`Unpacking definition: ${definition.name}`);
      // Examples are already in the examples array, just remove the definition
      deleteDefinition(id);
    }
  };

  // Handler for creating a new definition
  const createDefinition = (newDefinition) => {
    const id = `${Date.now()}`; // TODO: Get this from the database.
    setDefinitions([...definitions, { ...newDefinition, id }]);
  };

  return (
    <Router>
      <div className="app">
        <Navigation />
        <Routes>
          <Route 
            path="/library" 
            element={
              <LibraryView 
                examples={examples}
                definitions={definitions}
                onRecordNewActions={handleRecordNewActions}
                onUpdateExample={updateExample}
                onDeleteExample={deleteExample}
                onUpdateDefinition={updateDefinition}
                onDeleteDefinition={deleteDefinition}
                onUnpackDefinition={unpackDefinition}
              />
            } 
          />
          <Route 
            path="/preview/:exampleId?" 
            element={
              <PreviewMode 
                examples={examples}
                onUpdateExample={updateExample}
                onDeleteExample={deleteExample}
              />
            } 
          />
          <Route 
            path="/definition-builder" 
            element={
              <DefinitionBuilder 
                examples={examples}
                definitions={definitions}
                onCreateDefinition={createDefinition}
              />
            } 
          />
          <Route path="/" element={<Navigate to="/library" replace />} />
        </Routes>
      </div>
      <canvas ref={canvasRef} width={16} height={16} style={{ display: 'none' }}></canvas>
    </Router>
  );
}

export default App;

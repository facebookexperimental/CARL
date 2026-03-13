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
        this._nativeInspector = nativeInspector;
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

    serialize() {
        const bytes = this._nativeExample.serialize();
        const jsBytes = new Uint8Array(bytes.size());
        for (let idx = 0; idx < jsBytes.length; ++idx) {
            jsBytes[idx] = bytes.get(idx);
        }
        bytes.delete();
        return jsBytes;
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

    getActionType() {
        return this._nativeDefinition.getActionType();
    }

    getExamplesCount() {
        return this._nativeDefinition.getExamplesCount();
    }

    getCounterexamplesCount() {
        return this._nativeDefinition.getCounterexamplesCount();
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

    stopRecording(id, metadata) {
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
        if (this.onExampleCreated) this.onExampleCreated(jsExample, metadata);
        return jsExample;
    }

    getActionTypesMap() {
        return this._actionTypes;
    }

    draftDefinition(actionTypeId, examples, counterexamples) {
        const definition = new this._carl.Definition(this._idToActionType.get(actionTypeId));
        examples.forEach(example => {
            definition.addExample(example._nativeExample);
        });
        counterexamples.forEach(example => {
            definition.addCounterexample(example._nativeExample);
        });
        return new CarlDefinition(definition);
    }

    finalizeDefinition(definition, metadata) {
        if (this.onDefinitionCreated) this.onDefinitionCreated(definition, metadata);
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

    tryDeserializeExample(jsBytes) {
        const nativeBytes = new this._carl.SerializedBytes();
        nativeBytes.resize(jsBytes.length);
        for (let idx = 0; idx < jsBytes.length; ++idx) {
            nativeBytes.set(idx, jsBytes[idx]);
        }
        const nativeExample = this._carl.Example.tryDeserialize(nativeBytes);
        nativeBytes.delete();
        return nativeExample ? new CarlExample(nativeExample) : null;
    }

    tryDeserializeDefinition(jsBytes) {
        const nativeBytes = new this._carl.SerializedBytes();
        nativeBytes.resize(jsBytes.length);
        for (let idx = 0; idx < jsBytes.length; ++idx) {
            nativeBytes.set(idx, jsBytes[idx]);
        }
        const nativeDefinition = this._carl.Definition.tryDeserialize(nativeBytes);
        nativeBytes.delete();
        return nativeDefinition ? new CarlDefinition(nativeDefinition) : null;
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
            request.onerror = () => reject("Error deleting database: we are in an unknown state");
            request.onsuccess = () => resolve();
        });
    }

    async addAsync(obj) {
        const transaction = this._db.transaction(SerializationsDB._OBJECT_STORE_NAME, "readwrite")
        const os = transaction.objectStore(SerializationsDB._OBJECT_STORE_NAME);

        return await new Promise((resolve, reject) => {
            const request = os.add(obj);
            request.onerror = () => reject("Error storing CARL object");
            request.onsuccess = () => resolve(request.result);
        });
    }

    async fetchAllAsync() {
        const transaction = this._db.transaction(SerializationsDB._OBJECT_STORE_NAME, "readonly")
        const os = transaction.objectStore(SerializationsDB._OBJECT_STORE_NAME);

        return await new Promise((resolve, reject) => {
            const request = os.getAll();
            request.onerror = () => reject("Error loading CARL objects");
            request.onsuccess = (evt) => {
                resolve(evt.target.result);
            };
        });
    }

    async updateAsync(id, updates) {
        const transaction = this._db.transaction(SerializationsDB._OBJECT_STORE_NAME, "readwrite")
        const os = transaction.objectStore(SerializationsDB._OBJECT_STORE_NAME);

        let fetched = await new Promise((resolve, reject) => {
            const request = os.get(id);
            request.onerror = () => reject("Error in get() request");
            request.onsuccess = (evt) => {
                resolve(evt.target.result);
            };
        });
        fetched = { ...fetched, ...updates };

        return await new Promise((resolve, reject) => {
            const request = os.put(fetched);
            request.onerror = () => reject("Error in put() request");
            request.onsuccess = () => resolve();
        });
    }

    async deleteAsync(id) {
        const transaction = this._db.transaction(SerializationsDB._OBJECT_STORE_NAME, "readwrite")
        const os = transaction.objectStore(SerializationsDB._OBJECT_STORE_NAME);

        return await new Promise((resolve, reject) => {
            const request = os.delete(id);
            request.onerror = () => reject("Error deleting CARL objects");
            request.onsuccess = () => resolve();
        });
    }
}

function App() {
    const canvasRef = React.useRef(null);
    const immersiveExperienceRef = React.useRef(null);
    const dbRef = React.useRef(null);
    const carlRef = React.useRef(null);
    const [examples, setExamples] = useState([]);
    const [definitions, setDefinitions] = useState([]);

    async function reloadFromDatabase() {
        let xmpls = [];
        let defs = [];

        const records = await dbRef.current.fetchAllAsync();
        records.forEach(record => {
            switch (record.type) {
                case "example":
                    xmpls.push(record);
                    break;
                case "definition":
                    defs.push(record);
                    break;
                default:
                    console.error("Unknown record:");
                    console.error(record);
                    break;
            }
        });

        setExamples(xmpls);
        setDefinitions(defs);
    }

    async function bootstrapAsync() {
        dbRef.current = await SerializationsDB.loadAsync();
        carlRef.current = await CarlIntegration.CreateAsync();
        reloadFromDatabase();

        carlRef.current.onExampleCreated = (example, metadata) => {
            const inspector = example.getRecordingInspector();
            const recStart = inspector.getStartTimestamp();
            const recEnd = inspector.getEndTimestamp();
            inspector.dispose();
            const start = example.getStartTimestamp();
            const end = example.getEndTimestamp();
            console.log(metadata);
            dbRef.current.addAsync({
                type: "example",
                name: `Example ${Date.now()}`,
                recordingStart: recStart,
                recordingEnd: recEnd,
                exampleStart: start,
                exampleEnd: end,
                startTime: start - recStart,
                endTime: end - recStart,
                duration: recEnd - recStart,
                color: '#FF0000',
                showInXR: true,
                bytes: example.serialize(),
                ...metadata,
            }).then(id => {
                reloadFromDatabase();
                example.onDisposed = () => {
                    dbRef.current.deleteAsync(id);
                    reloadFromDatabase();
                };
            });
        };
        carlRef.current.onDefinitionCreated = (definition, metadata) => {
            dbRef.current.addAsync({
                type: "definition",
                name: `Definition ${Date.now()}`,
                actionType: carlRef.current.getActionTypesMap()[definition.getActionType()].name,
                examplesCount: definition.getExamplesCount(),
                counterexamplesCount: definition.getCounterexamplesCount(),
                defaultSensitivity: definition.getDefaultSensitivity(),
                color: '#FF0000',
                showInXR: true,
                bytes: definition.serialize(),
                ...metadata,
            }).then(id => {
                reloadFromDatabase();
                definition.onDisposed = () => {
                    dbRef.current.deleteAsync(id);
                    reloadFromDatabase();
                };
            });
        };
    }
    React.useEffect(() => {
        bootstrapAsync();
    }, []);

    const handleRecordNewActions = async () => {
        if (!carlRef.current || !canvasRef.current) return;
        if (immersiveExperienceRef.current) return;

        const initialExamples = [];
        for (const record of examples) {
            if (!record.showInXR) continue;
            const deserialized = carlRef.current.tryDeserializeExample(record.bytes);
            if (deserialized) {
                initialExamples.push({ value: deserialized, color: record.color });
            }
        }

        const initialDefinitions = [];
        for (const record of definitions) {
            if (!record.showInXR) continue;
            const deserialized = carlRef.current.tryDeserializeDefinition(record.bytes);
            if (deserialized) {
                initialDefinitions.push({ value: deserialized, color: record.color });
            }
        }

        const onSessionEnded = () => {
            immersiveExperienceRef.current = null;
            reloadFromDatabase();
        };

        immersiveExperienceRef.current = await initializeImmersiveExperienceAsync(
            canvasRef.current,
            carlRef.current,
            { initialExamples, initialDefinitions, onSessionEnded },
        );
    };

    // Handler for updating an example
    const updateExample = (id, updates) => {
        if (dbRef.current) {
            examples.map((example) => {
                if (example.id === id) {
                    let shouldReserialize = false;
                    if (updates.startTime) {
                        updates.exampleStart = updates.startTime + example.recordingStart;
                        shouldReserialize = true;
                    }

                    if (updates.endTime) {
                        updates.exampleEnd = updates.endTime + example.recordingStart;
                        shouldReserialize = true;
                    }

                    if (shouldReserialize) {
                        let ex = carlRef.current.tryDeserializeExample(example.bytes);
                        ex.setStartTimestamp(updates.exampleStart ?? example.exampleStart);
                        ex.setEndTimestamp(updates.exampleEnd ?? example.exampleEnd);
                        updates.bytes = ex.serialize();
                        ex.dispose();
                    }

                    dbRef.current.updateAsync(id, updates).then(() => {
                        reloadFromDatabase();
                    });
                }
            });
        }
    };

    // Handler for downloading an example
    const downloadExample = (name, bytes) => {
        if (dbRef.current) {
            downloadSerializedBytes(bytes, `example_${name}.carl`);
        }
    }

    // Handler for deleting an example
    const deleteExample = (id) => {
        if (dbRef.current) {
            dbRef.current.deleteAsync(id).then(() => {
                reloadFromDatabase();
            });
        }
    };

    // Handler for updating a definition
    const updateDefinition = (id, updates) => {
        if (dbRef.current) {
            dbRef.current.updateAsync(id, updates).then(() => {
                reloadFromDatabase();
            });
        }
    };

    // Handler for downloading a definition
    const downloadDefinition = (name, bytes) => {
        if (dbRef.current) {
            downloadSerializedBytes(bytes, `definition_${name}.carl`);
        }
    }

    // Handler for deleting a definition
    const deleteDefinition = (id) => {
        if (dbRef.current) {
            dbRef.current.deleteAsync(id).then(() => {
                reloadFromDatabase();
            });
        }
    };

    // Handler for unpacking a definition
    const unpackDefinition = (id) => {
        // const definition = definitions.find(def => def.id === id);
        // if (definition) {
        //     console.log(`Unpacking definition: ${definition.name}`);
        //     // Examples are already in the examples array, just remove the definition
        //     deleteDefinition(id);
        // }
        console.error("TODO: Implement");
    };

    const handleFilesDropped = async (files) => {
        if (!dbRef.current || !carlRef.current) {
            // nothing to do if integration or DB isn't ready yet
            return;
        }

        // helper that reads a file into a Uint8Array and tries to parse it
        const processFile = (file) => {
            return new Promise(resolve => {
                const reader = new FileReader();
                reader.onload = () => {
                    const arrayBuffer = reader.result;
                    const array = new Uint8Array(arrayBuffer);
                    let parsed = null;
                    let type = null;

                    const example = carlRef.current.tryDeserializeExample(array);
                    if (example) {
                        parsed = example;
                        type = "example";
                    } else {
                        const definition = carlRef.current.tryDeserializeDefinition(array);
                        if (definition) {
                            parsed = definition;
                            type = "definition";
                        }
                    }

                    resolve({ file, type, parsed });
                };
                reader.onerror = () => {
                    console.error(`Failed reading dropped file: ${file.name}`);
                    resolve({ file, type: null, parsed: null });
                };
                reader.readAsArrayBuffer(file);
            });
        };

        const results = await Promise.all(Array.from(files).map(processFile));

        for (const { file, type, parsed } of results) {
            if (!parsed) {
                continue; // skip unknown blobs
            }

            if (type === "example") {
                const example = parsed; // CarlExample
                const inspector = example.getRecordingInspector();
                const recStart = inspector.getStartTimestamp();
                const recEnd = inspector.getEndTimestamp();
                inspector.dispose();
                const start = example.getStartTimestamp();
                const end = example.getEndTimestamp();

                await dbRef.current.addAsync({
                    type: "example",
                    name: file.name,
                    recordingStart: recStart,
                    recordingEnd: recEnd,
                    exampleStart: start,
                    exampleEnd: end,
                    startTime: start - recStart,
                    endTime: end - recStart,
                    duration: recEnd - recStart,
                    color: '#FF0000',
                    showInXR: true,
                    bytes: example.serialize(),
                    fileName: file.name,
                }).then(id => {
                    example.onDisposed = () => {
                        dbRef.current.deleteAsync(id);
                        reloadFromDatabase();
                    };
                });
            } else if (type === "definition") {
                const definition = parsed; // CarlDefinition
                await dbRef.current.addAsync({
                    type: "definition",
                    name: file.name,
                    actionType: carlRef.current.getActionTypesMap()[definition.getActionType()].name,
                    examplesCount: definition.getExamplesCount(),
                    counterexamplesCount: definition.getCounterexamplesCount(),
                    defaultSensitivity: definition.getDefaultSensitivity(),
                    color: '#FF0000',
                    showInXR: true,
                    bytes: definition.serialize(),
                    fileName: file.name,
                }).then(id => {
                    definition.onDisposed = () => {
                        dbRef.current.deleteAsync(id);
                        reloadFromDatabase();
                    };
                });
            }
        }

        // once all files have been persisted, refresh view
        reloadFromDatabase();
    };

    // Handler for creating a new definition
    const createDefinition = (newDefinition) => {
        // const id = `${Date.now()}`;
        // setDefinitions([...definitions, { ...newDefinition, id }]);
        console.error("TODO: Implement");
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
                                onDownloadExample={downloadExample}
                                onDeleteExample={deleteExample}
                                onUpdateDefinition={updateDefinition}
                                onDownloadDefinition={downloadDefinition}
                                onDeleteDefinition={deleteDefinition}
                                onUnpackDefinition={unpackDefinition}
                                onFilesDropped={handleFilesDropped}
                            />
                        }
                    />
                    <Route
                        path="/preview/:exampleId?"
                        element={ examples.length > 0 ?
                            <PreviewMode
                                examples={examples}
                                onUpdateExample={updateExample}
                                onDeleteExample={deleteExample}
                            /> :
                            <Navigate to="/library" replace />
                        }
                    />
                    <Route
                        path="/definition-builder"
                        element={ examples.length > 0 ?
                            <DefinitionBuilder
                                examples={examples}
                                definitions={definitions}
                                onCreateDefinition={createDefinition}
                            /> :
                            <Navigate to="/library" replace />
                        }
                    />
                    <Route path="/" element={<Navigate to="/library" replace />} />
                    <Route path="*" element={<Navigate to="/library" replace />} />
                </Routes>
            </div>
            <canvas ref={canvasRef} width={16} height={16} style={{ display: 'none' }}></canvas>
        </Router>
    );
}

export default App;

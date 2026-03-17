/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * Application root.  Owns routing, top-level React state, IndexedDB
 * bootstrapping, and all handler functions passed down to child routes.
 *
 * @module App
 */

import React, { useState } from 'react';
import { BrowserRouter as Router, Routes, Route, Navigate } from 'react-router-dom';
import Navigation from './components/Navigation';
import LibraryView from './components/LibraryView';
import PreviewMode from './components/PreviewMode';
import DefinitionBuilder from './components/DefinitionBuilder';
import { initializeImmersiveExperienceAsync } from 'carl-actionstudio-immersiveexperience';
import { downloadSerializedBytes } from './lib/utils.js';
import { SerializationsDB } from './lib/serializationsDB.js';
import { CarlIntegration } from './lib/carlIntegration.js';
import './styles/App.css';

function App() {
    // --- Refs ---
    const canvasRef = React.useRef(null);
    const immersiveExperienceRef = React.useRef(null);
    const dbRef = React.useRef(null);
    const carlRef = React.useRef(null);

    // --- State ---
    const [examples, setExamples] = useState([]);
    const [definitions, setDefinitions] = useState([]);
    const [isLoaded, setIsLoaded] = useState(false);

    // --- Database bootstrap ---
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
        setIsLoaded(true);
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
                immersiveExperienceRef.current?.notifyExampleSaved(example, id);
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

    // --- Example handlers ---
    const handleRecordNewActions = async () => {
        if (!carlRef.current || !canvasRef.current) return;
        if (immersiveExperienceRef.current) return;

        const initialExamples = [];
        for (const record of examples) {
            if (!record.showInXR) continue;
            const deserialized = carlRef.current.tryDeserializeExample(record.bytes);
            if (deserialized) {
                initialExamples.push({ value: deserialized, color: record.color, id: record.id });
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

        const onExampleUpdated = (id, example) => {
            if (!dbRef.current) return;
            const inspector = example.getRecordingInspector();
            const recordingStart = inspector.getStartTimestamp();
            inspector.dispose();
            const newStart = example.getStartTimestamp();
            const newEnd = example.getEndTimestamp();
            dbRef.current.updateAsync(id, {
                exampleStart: newStart,
                exampleEnd: newEnd,
                startTime: newStart - recordingStart,
                endTime: newEnd - recordingStart,
                bytes: example.serialize(),
            }).then(() => reloadFromDatabase());
        };

        immersiveExperienceRef.current = await initializeImmersiveExperienceAsync(
            canvasRef.current,
            carlRef.current,
            { initialExamples, initialDefinitions, onSessionEnded, onExampleUpdated },
        );
    };

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

    const downloadExample = (name, bytes) => {
        if (dbRef.current) {
            downloadSerializedBytes(bytes, `example_${name}.carl`);
        }
    }

    const deleteExample = (id) => {
        if (dbRef.current) {
            dbRef.current.deleteAsync(id).then(() => {
                reloadFromDatabase();
            });
        }
    };

    // --- Definition handlers ---
    const updateDefinition = (id, updates) => {
        if (dbRef.current) {
            dbRef.current.updateAsync(id, updates).then(() => {
                reloadFromDatabase();
            });
        }
    };

    const downloadDefinition = (name, bytes) => {
        if (dbRef.current) {
            downloadSerializedBytes(bytes, `definition_${name}.carl`);
        }
    }

    const deleteDefinition = (id) => {
        if (dbRef.current) {
            dbRef.current.deleteAsync(id).then(() => {
                reloadFromDatabase();
            });
        }
    };

    const unpackDefinition = async (id) => {
        const record = definitions.find(d => d.id === id);
        if (!record || !carlRef.current) return;

        const definition = carlRef.current.tryDeserializeDefinition(record.bytes);
        if (!definition) return;

        const examplesCount = definition.getExamplesCount();
        const counterexamplesCount = definition.getCounterexamplesCount();

        for (let i = 0; i < examplesCount; i++) {
            const example = definition.getExample(i);
            const inspector = example.getRecordingInspector();
            const recStart = inspector.getStartTimestamp();
            const recEnd = inspector.getEndTimestamp();
            inspector.dispose();
            const start = example.getStartTimestamp();
            const end = example.getEndTimestamp();
            await dbRef.current.addAsync({
                type: "example",
                name: `${record.name}_example_${i}`,
                recordingStart: recStart,
                recordingEnd: recEnd,
                exampleStart: start,
                exampleEnd: end,
                startTime: start - recStart,
                endTime: end - recStart,
                duration: recEnd - recStart,
                color: record.color,
                showInXR: record.showInXR,
                bytes: example.serialize(),
            });
            example.dispose();
        }

        for (let i = 0; i < counterexamplesCount; i++) {
            const example = definition.getCounterexample(i);
            const inspector = example.getRecordingInspector();
            const recStart = inspector.getStartTimestamp();
            const recEnd = inspector.getEndTimestamp();
            inspector.dispose();
            const start = example.getStartTimestamp();
            const end = example.getEndTimestamp();
            await dbRef.current.addAsync({
                type: "example",
                name: `${record.name}_counterexample_${i}`,
                recordingStart: recStart,
                recordingEnd: recEnd,
                exampleStart: start,
                exampleEnd: end,
                startTime: start - recStart,
                endTime: end - recStart,
                duration: recEnd - recStart,
                color: record.color,
                showInXR: record.showInXR,
                bytes: example.serialize(),
            });
            example.dispose();
        }

        definition.dispose();
        reloadFromDatabase();
    };

    const createDefinition = (newDefinition) => {
        if (!carlRef.current) return;

        const actionTypeId = carlRef.current.getActionTypesMap().find(t => t.name === newDefinition.actionType)?.typeId;
        if (actionTypeId === undefined) {
            console.error(`Unknown action type: ${newDefinition.actionType}`);
            return;
        }

        const exampleObjects = newDefinition.examples
            .map(id => examples.find(ex => ex.id === id))
            .filter(Boolean)
            .map(r => carlRef.current.tryDeserializeExample(r.bytes))
            .filter(Boolean);

        const counterexampleObjects = newDefinition.counterexamples
            .map(id => examples.find(ex => ex.id === id))
            .filter(Boolean)
            .map(r => carlRef.current.tryDeserializeExample(r.bytes))
            .filter(Boolean);

        const definition = carlRef.current.draftDefinition(actionTypeId, exampleObjects, counterexampleObjects);
        definition.setDefaultSensitivity(newDefinition.defaultSensitivity);
        carlRef.current.finalizeDefinition(definition, {
            name: newDefinition.name,
            color: newDefinition.color,
            showInXR: newDefinition.showInXR,
        });

        [...exampleObjects, ...counterexampleObjects].forEach(e => e.dispose());
    };

    // --- File drop handler ---
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

    // --- Routing ---
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
                        element={ !isLoaded ? null : examples.length > 0 ?
                            <PreviewMode
                                examples={examples}
                                onUpdateExample={updateExample}
                                onDeleteExample={deleteExample}
                                carl={carlRef.current}
                            /> :
                            <Navigate to="/library" replace />
                        }
                    />
                    <Route
                        path="/definition-builder"
                        element={ !isLoaded ? null : examples.length > 0 ?
                            <DefinitionBuilder
                                examples={examples}
                                definitions={definitions}
                                onCreateDefinition={createDefinition}
                                carl={carlRef.current}
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

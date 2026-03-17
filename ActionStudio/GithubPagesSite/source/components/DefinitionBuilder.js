/**
 * Full-page builder for creating a new CARL Definition from a set of
 * examples and counterexamples.  Also provides an inline test runner that
 * plots recognizer response curves using real CARL scoring.
 *
 * @module components/DefinitionBuilder
 */

import React, { useState, useEffect } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import { actionTypes } from '../data/mockData';
import { formatDuration } from '../lib/utils.js';
import '../styles/DefinitionBuilder.css';

function DefinitionBuilder({ examples, definitions, onCreateDefinition, carl }) {
  const navigate = useNavigate();

  // --- State ---
  const [name, setName] = useState('New Definition');
  const [actionType, setActionType] = useState(actionTypes[0]);
  const [sensitivity, setSensitivity] = useState(5);
  const [color, setColor] = useState('#5B9FFF');
  const [selectedExamples, setSelectedExamples] = useState([]);
  const [selectedCounterexamples, setSelectedCounterexamples] = useState([]);
  const [testExamples, setTestExamples] = useState([]);
  const [testResults, setTestResults] = useState([]);

  // --- Handlers ---
  const handleSave = () => {
    const definitionData = {
      name: name.trim(),
      actionType,
      examples: selectedExamples,
      counterexamples: selectedCounterexamples,
      defaultSensitivity: sensitivity,
      color,
      showInXR: true,
    };

    onCreateDefinition(definitionData);
    
    navigate('/library');
  };

  const handleCancel = () => {
    navigate('/library');
  };

  const toggleExample = (exampleId) => {
    setSelectedExamples(prev => 
      prev.includes(exampleId) 
        ? prev.filter(id => id !== exampleId)
        : [...prev, exampleId]
    );
  };

  const toggleCounterexample = (exampleId) => {
    setSelectedCounterexamples(prev => 
      prev.includes(exampleId) 
        ? prev.filter(id => id !== exampleId)
        : [...prev, exampleId]
    );
  };

  const handleRunTest = () => {
    if (testExamples.length === 0) {
      alert('Please select examples to test');
      return;
    }
    if (selectedExamples.length === 0) {
      alert('Please select at least one example to build the definition before testing');
      return;
    }
    if (!carl) return;

    const actionTypeId = carl.getActionTypesMap().find(t => t.name === actionType)?.typeId;
    if (actionTypeId === undefined) return;

    const exampleObjects = selectedExamples
      .map(id => examples.find(ex => ex.id === id))
      .filter(Boolean)
      .map(r => carl.tryDeserializeExample(r.bytes))
      .filter(Boolean);

    const counterexampleObjects = selectedCounterexamples
      .map(id => examples.find(ex => ex.id === id))
      .filter(Boolean)
      .map(r => carl.tryDeserializeExample(r.bytes))
      .filter(Boolean);
      
    const tempDefinition = carl.draftDefinition(actionTypeId, exampleObjects, counterexampleObjects);
    tempDefinition.setDefaultSensitivity(sensitivity);

    const testCarlEntries = testExamples
      .map(id => ({ id, record: examples.find(ex => ex.id === id) }))
      .filter(e => e.record)
      .map(e => ({ ...e, carlExample: carl.tryDeserializeExample(e.record.bytes) }))
      .filter(e => e.carlExample);

    let rawResults;
    try {
      rawResults = carl.testDefinition(tempDefinition, testCarlEntries.map(e => e.carlExample), sensitivity);
    } catch (err) {
      console.error('testDefinition failed:', err);
      [...exampleObjects, ...counterexampleObjects, ...testCarlEntries.map(e => e.carlExample)].forEach(e => e.dispose());
      tempDefinition.dispose?.();
      return;
    }

    const results = testCarlEntries.map((entry, idx) => ({
      exampleId:    entry.id,
      exampleName:  entry.record.name,
      exampleColor: entry.record.color,
      data: rawResults[idx] ?? [],
    }));

    setTestResults(results);

    [...exampleObjects, ...counterexampleObjects, ...testCarlEntries.map(e => e.carlExample)].forEach(e => e.dispose());
    tempDefinition.dispose?.();
  };

  const handleTestExampleToggle = (exampleId) => {
    setTestExamples(prev =>
      prev.includes(exampleId)
        ? prev.filter(id => id !== exampleId)
        : [...prev, exampleId]
    );
  };

  // --- Render ---
  return (
    <div className="definition-builder">
      <div className="builder-header">
        <input
          type="text"
          className="definition-name-input"
          value={name}
          onChange={(e) => setName(e.target.value)}
          placeholder="Definition Name"
        />
        <div className="header-actions">
          <button className="cancel-btn" onClick={handleCancel}>Cancel</button>
          <button className="save-btn" onClick={handleSave}>
            ✓ Save Definition
          </button>
        </div>
      </div>

      <div className="builder-content">
        <aside className="config-panel">
          <h3>Configuration</h3>
          
          <div className="config-section">
            <label>Action Type</label>
            <select 
              className="action-type-select"
              value={actionType}
              onChange={(e) => setActionType(e.target.value)}
            >
              {actionTypes.map(type => (
                <option key={type} value={type}>{type}</option>
              ))}
            </select>
          </div>

          <div className="config-section">
            <label>Default Sensitivity</label>
            <div className="sensitivity-control">
              <input
                type="range"
                className="sensitivity-slider"
                min="0"
                max="10"
                step="0.1"
                value={sensitivity}
                style={{ background: `linear-gradient(to right, #5B9FFF 0%, #5B9FFF ${sensitivity * 10}%, #2a2a3e ${sensitivity * 10}%, #2a2a3e 100%)` }}
                onChange={(e) => setSensitivity(parseFloat(e.target.value))}
              />
              <input
                type="number"
                className="sensitivity-input"
                min="0"
                max="10"
                step="0.1"
                value={sensitivity}
                onChange={(e) => setSensitivity(parseFloat(e.target.value))}
              />
            </div>
            <div className="sensitivity-labels">
              <span>0</span>
              <span>5</span>
              <span>10</span>
            </div>
          </div>

          <div className="config-section">
            <label>Color</label>
            <input
              type="color"
              className="color-picker"
              value={color}
              onChange={(e) => setColor(e.target.value)}
            />
          </div>
        </aside>

        <main className="selection-panel">
          <section className="selection-section">
            <div className="section-header">
              <h3>Examples</h3>
              <span className="count-badge">{selectedExamples.length}</span>
            </div>
            <div className="example-grid">
              {examples.map(example => (
                <div
                  key={example.id}
                  className={`selectable-card ${selectedExamples.includes(example.id) ? 'selected' : ''}`}
                  onClick={() => toggleExample(example.id)}
                  style={{ borderLeftColor: example.color }}
                >
                  {selectedExamples.includes(example.id) && (
                    <div className="selection-check">✓</div>
                  )}
                  <div className="card-preview-small" style={{ background: `linear-gradient(135deg, ${example.color}22, ${example.color}44)` }}>
                    <span className="preview-icon-small">🎬</span>
                  </div>
                  <div className="card-info">
                    <span className="card-name-small">{example.name}</span>
                    <span className="card-duration-small">{formatDuration(example.duration)}</span>
                  </div>
                </div>
              ))}
            </div>
          </section>

          <section className="selection-section">
            <div className="section-header">
              <h3>Counterexamples</h3>
              <span className="count-badge counterexample-badge">{selectedCounterexamples.length}</span>
            </div>
            <div className="example-grid">
              {examples.map(example => (
                <div
                  key={example.id}
                  className={`selectable-card ${selectedCounterexamples.includes(example.id) ? 'selected counterexample' : ''}`}
                  onClick={() => toggleCounterexample(example.id)}
                  style={{ borderLeftColor: example.color }}
                >
                  {selectedCounterexamples.includes(example.id) && (
                    <div className="selection-check counterexample-check">✓</div>
                  )}
                  <div className="card-preview-small" style={{ background: `linear-gradient(135deg, ${example.color}22, ${example.color}44)` }}>
                    <span className="preview-icon-small">🎬</span>
                  </div>
                  <div className="card-info">
                    <span className="card-name-small">{example.name}</span>
                    <span className="card-duration-small">{formatDuration(example.duration)}</span>
                  </div>
                </div>
              ))}
            </div>
          </section>
        </main>

        <aside className="testing-panel">
          <h3>Testing</h3>
          
          <div className="test-controls">
            <label>Select Examples to Test</label>
            <div className="test-example-list">
              {examples.map(example => (
                <label key={example.id} className="test-example-item">
                  <input
                    type="checkbox"
                    checked={testExamples.includes(example.id)}
                    onChange={() => handleTestExampleToggle(example.id)}
                  />
                  <span style={{ color: example.color }}>●</span>
                  <span>{example.name}</span>
                </label>
              ))}
            </div>
            <button className="run-test-btn" onClick={handleRunTest}>
              Run Test
            </button>
          </div>

          <div className="test-results">
            {testResults.length > 0 ? (
              <div className="graph-container">
                <svg className="test-graph" viewBox="0 0 400 200">
                  {/* Grid lines */}
                  <line x1="40" y1="20" x2="40" y2="180" stroke="#444" strokeWidth="1" />
                  <line x1="40" y1="180" x2="380" y2="180" stroke="#444" strokeWidth="1" />
                  
                  {/* Y-axis labels */}
                  <text x="30" y="25" fill="#888" fontSize="10">1</text>
                  <text x="30" y="105" fill="#888" fontSize="10">.5</text>
                  <text x="30" y="185" fill="#888" fontSize="10">0</text>
                  
                  {/* X-axis label */}
                  <text x="200" y="195" fill="#888" fontSize="10" textAnchor="middle">Time (s)</text>
                  <text x="20" y="100" fill="#888" fontSize="10" transform="rotate(-90 20 100)">Score</text>
                  
                  {/* Half score line */}
                  <line 
                    x1="40" 
                    y1="100"
                    x2="380" 
                    y2="100" 
                    stroke="#FBBF24" 
                    strokeWidth="1" 
                    strokeDasharray="4"
                  />
                  
                  {/* Plot data */}
                  {testResults.map((result, idx) => {
                    const points = result.data.map((point, i) => {
                      const x = 40 + (point.time / result.data[result.data.length - 1].time) * 340;
                      const y = 180 - point.score * 160;
                      return `${x},${y}`;
                    }).join(' ');
                    
                    return (
                      <polyline
                        key={result.exampleId}
                        points={points}
                        fill="none"
                        stroke={result.exampleColor}
                        strokeWidth="2"
                      />
                    );
                  })}
                </svg>
                
                <div className="graph-legend">
                  {testResults.map(result => (
                    <div key={result.exampleId} className="legend-item">
                      <span className="legend-color" style={{ backgroundColor: result.exampleColor }}></span>
                      <span className="legend-label">{result.exampleName}</span>
                    </div>
                  ))}
                </div>
              </div>
            ) : (
              <div className="no-results">
                <p>Select examples and run a test to see results</p>
              </div>
            )}
          </div>
        </aside>
      </div>
    </div>
  );
}

export default DefinitionBuilder;

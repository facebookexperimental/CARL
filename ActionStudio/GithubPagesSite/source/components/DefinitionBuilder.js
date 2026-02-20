import React, { useState, useEffect } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import { actionTypes } from '../data/mockData';
import '../styles/DefinitionBuilder.css';

function DefinitionBuilder({ examples, definitions, onCreateDefinition, onUpdateDefinition }) {
  const { definitionId } = useParams();
  const navigate = useNavigate();
  
  const [name, setName] = useState('New Definition');
  const [actionType, setActionType] = useState('Strike');
  const [sensitivity, setSensitivity] = useState(7.5);
  const [color, setColor] = useState('#5B9FFF');
  const [selectedExamples, setSelectedExamples] = useState([]);
  const [selectedCounterexamples, setSelectedCounterexamples] = useState([]);
  const [testExamples, setTestExamples] = useState([]);
  const [testResults, setTestResults] = useState([]);

  useEffect(() => {
    if (definitionId) {
      const definition = definitions.find(def => def.id === definitionId);
      if (definition) {
        setName(definition.name);
        setActionType(definition.actionType);
        setSensitivity(definition.defaultSensitivity);
        setColor(definition.color);
        setSelectedExamples(definition.examples);
        setSelectedCounterexamples(definition.counterexamples);
      }
    }
  }, [definitionId, definitions]);

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

    if (definitionId) {
      onUpdateDefinition(definitionId, definitionData);
    } else {
      onCreateDefinition(definitionData);
    }
    
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

    // Generate mock test results
    const results = testExamples.map(exId => {
      const example = examples.find(ex => ex.id === exId);
      // Generate mock sensitivity curve data
      const dataPoints = [];
      const duration = example.duration;
      const steps = 50;
      
      for (let i = 0; i <= steps; i++) {
        const time = (duration * i) / steps;
        // Mock sensitivity response with some randomness
        const baseResponse = Math.sin((i / steps) * Math.PI * 2) * 3 + sensitivity;
        const noise = (Math.random() - 0.5) * 1.5;
        const response = Math.max(0, Math.min(10, baseResponse + noise));
        dataPoints.push({ time, response });
      }
      
      return {
        exampleId: exId,
        exampleName: example.name,
        exampleColor: example.color,
        data: dataPoints,
      };
    });

    setTestResults(results);
  };

  const handleTestExampleToggle = (exampleId) => {
    setTestExamples(prev =>
      prev.includes(exampleId)
        ? prev.filter(id => id !== exampleId)
        : [...prev, exampleId]
    );
  };

  const formatDuration = (seconds) => {
    const mins = Math.floor(seconds / 60);
    const secs = Math.floor(seconds % 60);
    return `${mins.toString().padStart(2, '0')}:${secs.toString().padStart(2, '0')}`;
  };

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
                  <text x="25" y="25" fill="#888" fontSize="10">10</text>
                  <text x="30" y="105" fill="#888" fontSize="10">5</text>
                  <text x="30" y="185" fill="#888" fontSize="10">0</text>
                  
                  {/* X-axis label */}
                  <text x="200" y="195" fill="#888" fontSize="10" textAnchor="middle">Time (s)</text>
                  <text x="20" y="100" fill="#888" fontSize="10" transform="rotate(-90 20 100)">Sensitivity</text>
                  
                  {/* Sensitivity threshold line */}
                  <line 
                    x1="40" 
                    y1={180 - (sensitivity / 10) * 160} 
                    x2="380" 
                    y2={180 - (sensitivity / 10) * 160} 
                    stroke="#FBBF24" 
                    strokeWidth="1" 
                    strokeDasharray="4"
                  />
                  
                  {/* Plot data */}
                  {testResults.map((result, idx) => {
                    const points = result.data.map((point, i) => {
                      const x = 40 + (point.time / result.data[result.data.length - 1].time) * 340;
                      const y = 180 - (point.response / 10) * 160;
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

          <details className="test-history">
            <summary>Test History</summary>
            <p className="history-placeholder">Previous test results will appear here</p>
          </details>
        </aside>
      </div>
    </div>
  );
}

export default DefinitionBuilder;

import React, { useState } from 'react';
import '../styles/Card.css';

function DefinitionCard({ definition, examples, onUpdate, onDelete, onUnpack }) {
  const [isEditing, setIsEditing] = useState(false);
  const [name, setName] = useState(definition.name);

  const handleNameSubmit = () => {
    if (name.trim() && name !== definition.name) {
      onUpdate({ name: name.trim() });
    }
    setIsEditing(false);
  };

  const handleNameKeyDown = (e) => {
    if (e.key === 'Enter') {
      handleNameSubmit();
    } else if (e.key === 'Escape') {
      setName(definition.name);
      setIsEditing(false);
    }
  };

  const handleDownload = () => {
    console.log(`Downloading definition: ${definition.name}`);
    // TODO: Implement actual download functionality
    alert(`Download functionality will be implemented with your data format.\nDefinition: ${definition.name}`);
  };

  const handleUnpack = () => {
    if (window.confirm(`Unpack "${definition.name}"? This will extract all examples and delete the definition.`)) {
      onUnpack();
    }
  };

  return (
    <div className="card definition-card" style={{ borderLeftColor: definition.color }}>
      <div className="card-content">
        {isEditing ? (
          <input
            type="text"
            className="card-name-input"
            value={name}
            onChange={(e) => setName(e.target.value)}
            onBlur={handleNameSubmit}
            onKeyDown={handleNameKeyDown}
            autoFocus
          />
        ) : (
          <h3 className="card-name" onClick={() => setIsEditing(true)}>
            {definition.name}
          </h3>
        )}
        
        <div className="definition-meta">
          <div className="meta-row">
            <span className="action-type-badge">{definition.actionType}</span>
            <span className="sensitivity">Sensitivity: {definition.defaultSensitivity.toFixed(1)}</span>
          </div>
          <div className="meta-row">
            <span className="example-count">
              {definition.examplesCount} examples, {definition.counterexamplesCount} counterexamples
            </span>
          </div>
        </div>
      </div>

      <div className="card-actions">
        <button 
          className="action-btn" 
          onClick={() => onUpdate({ showInXR: !definition.showInXR })}
          title={definition.showInXR ? "Hide in XR" : "Show in XR"}
        >
          {definition.showInXR ? '👁' : '🚫'}
        </button>
        <button 
          className="action-btn" 
          onClick={handleDownload}
          title="Download"
        >
          ⬇
        </button>
        <button 
          className="action-btn" 
          onClick={handleUnpack}
          title="Unpack"
        >
          📦
        </button>
        {/* <button // TODO: Replace this with a "compress" feature?
          className="action-btn" 
          onClick={onEdit}
          title="Edit"
        >
          ✏️
        </button> */}
        <button 
          className="action-btn danger" 
          onClick={onDelete}
          title="Delete"
        >
          🗑
        </button>
      </div>
    </div>
  );
}

export default DefinitionCard;

/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * Card component for a single CARL Definition in the library grid.
 * Displays action type, sensitivity, and example/counterexample counts.
 * Supports inline rename, XR visibility toggle, download, unpack, and delete.
 *
 * @module components/DefinitionCard
 */

import React, { useState } from 'react';
import '../styles/Card.css';

function DefinitionCard({ definition, examples, onUpdate, onDownload, onDelete, onUnpack }) {
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

  const handleUnpack = () => {
    if (window.confirm(`Unpack "${definition.name}"? This will extract all examples and counterexamples as individual records.`)) {
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
          onClick={onDownload}
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

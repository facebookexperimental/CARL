/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * Main library view showing all stored examples and definitions.
 * Provides search/filter controls, drag-and-drop import of `.carl` files,
 * and navigation to preview and definition-builder routes.
 *
 * @module components/LibraryView
 */

import React, { useState } from 'react';
import { useNavigate } from 'react-router-dom';
import ExampleCard from './ExampleCard';
import DefinitionCard from './DefinitionCard';
import '../styles/LibraryView.css';

function LibraryView({ 
  examples, 
  definitions, 
  onRecordNewActions,
  onUpdateExample,
  onDownloadExample,
  onDeleteExample,
  onUpdateDefinition,
  onDownloadDefinition,
  onDeleteDefinition,
  onUnpackDefinition,
  onFilesDropped,
}) {
  const navigate = useNavigate();

  // --- Filtering ---
  const [searchQuery, setSearchQuery] = useState('');

  // --- Drag and drop ---
  const [isDragging, setIsDragging] = useState(false);

  const handleDragOver = (e) => {
    e.preventDefault();
    setIsDragging(true);
  };

  const handleDragLeave = (e) => {
    e.preventDefault();
    setIsDragging(false);
  };

  const handleDrop = (e) => {
    e.preventDefault();
    setIsDragging(false);
    const files = Array.from(e.dataTransfer.files);
    console.log('Files dropped:', files);
    if (onFilesDropped) {
      onFilesDropped(files);
    }
  };

  const handlePreviewExample = (exampleId) => {
    navigate(`/preview/${exampleId}`);
  };

  const handleCreateDefinition = () => {
    navigate('/definition-builder');
  };

  // --- Render ---
  const filteredExamples = examples.filter(ex =>
    ex.name.toLowerCase().includes(searchQuery.toLowerCase())
  );

  const filteredDefinitions = definitions.filter(def => {
    const matchesSearch = def.name.toLowerCase().includes(searchQuery.toLowerCase());
    return matchesSearch;
  });

  return (
    <div 
      className="library-view"
      onDragOver={handleDragOver}
      onDragLeave={handleDragLeave}
      onDrop={handleDrop}
    >
      {isDragging && (
        <div className="drag-overlay">
          <div className="drag-message">
            Drop Actions or Definitions here
          </div>
        </div>
      )}

      <aside className="library-sidebar">
        <button className="record-button" onClick={onRecordNewActions}>
          <span className="record-icon">🥽</span>
          Record New Actions
        </button>

        <div className="filter-controls">
          <h3>Filter Controls</h3>
          
          <input
            type="text"
            className="search-input"
            placeholder="Search..."
            value={searchQuery}
            onChange={(e) => setSearchQuery(e.target.value)}
          />
        </div>
      </aside>

      <main className="library-main">
        <section className="library-section">
          <div className="section-header">
            <h2>Examples</h2>
            <span className="count-badge">{filteredExamples.length}</span>
          </div>
          <div className="card-grid">
            {filteredExamples.map(example => (
              <ExampleCard
                key={example.id}
                example={example}
                onPreview={() => handlePreviewExample(example.id)}
                onUpdate={(updates) => onUpdateExample(example.id, updates)}
                onDownload={() => onDownloadExample(example.name, example.bytes)}
                onDelete={() => onDeleteExample(example.id)}
              />
            ))}
          </div>
        </section>

        <section className="library-section">
          <div className="section-header">
            <h2>Definitions</h2>
            <span className="count-badge">{filteredDefinitions.length}</span>
            <button className="create-definition-btn" onClick={handleCreateDefinition}>
              + New Definition
            </button>
          </div>
          <div className="card-grid definitions-grid">
            {filteredDefinitions.map(definition => (
              <DefinitionCard
                key={definition.id}
                definition={definition}
                examples={examples}
                onUpdate={(updates) => onUpdateDefinition(definition.id, updates)}
                onDownload={() => onDownloadDefinition(definition.name, definition.bytes)}
                onDelete={() => onDeleteDefinition(definition.id)}
                onUnpack={() => onUnpackDefinition(definition.id)}
              />
            ))}
          </div>
        </section>
      </main>
    </div>
  );
}

export default LibraryView;

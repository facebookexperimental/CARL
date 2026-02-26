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
  onUnpackDefinition
}) {
  const navigate = useNavigate();
  const [searchQuery, setSearchQuery] = useState('');
  const [showArchived, setShowArchived] = useState(true);
  const [filterType, setFilterType] = useState('All Types');
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
    // TODO: Implement file parsing and import
    alert(`Dropped ${files.length} file(s). File import will be implemented with your data format.`);
  };

  const handlePreviewExample = (exampleId) => {
    navigate(`/preview/${exampleId}`);
  };

  const handleCreateDefinition = () => {
    navigate('/definition-builder');
  };

  const filteredExamples = examples.filter(ex => 
    ex.name.toLowerCase().includes(searchQuery.toLowerCase())
  );

  const filteredDefinitions = definitions.filter(def => {
    const matchesSearch = def.name.toLowerCase().includes(searchQuery.toLowerCase());
    const matchesType = filterType === 'All Types' || def.actionType === filterType;
    return matchesSearch && matchesType;
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

          <label className="toggle-label">
            <input
              type="checkbox"
              checked={showArchived}
              onChange={(e) => setShowArchived(e.target.checked)}
            />
            <span>Show Archived</span>
          </label>

          <label className="filter-label">
            Filter by Type
            <select 
              className="filter-select"
              value={filterType}
              onChange={(e) => setFilterType(e.target.value)}
            >
              <option>All Types</option>
              <option>Strike</option>
              <option>Grapple</option>
              <option>Defense</option>
              <option>Movement</option>
              <option>Gesture</option>
            </select>
          </label>
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

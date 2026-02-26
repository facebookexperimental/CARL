import React, { useState } from 'react';
import '../styles/Card.css';

function ExampleCard({ example, onPreview, onUpdate, onDownload, onDelete }) {
  const [isEditing, setIsEditing] = useState(false);
  const [name, setName] = useState(example.name);

  const handleNameSubmit = () => {
    if (name.trim() && name !== example.name) {
      onUpdate({ name: name.trim() });
    }
    setIsEditing(false);
  };

  const handleNameKeyDown = (e) => {
    if (e.key === 'Enter') {
      handleNameSubmit();
    } else if (e.key === 'Escape') {
      setName(example.name);
      setIsEditing(false);
    }
  };

  const formatDuration = (seconds) => {
    const mins = Math.floor(seconds / 60);
    const secs = Math.floor(seconds % 60);
    return `${mins.toString().padStart(2, '0')}:${secs.toString().padStart(2, '0')}`;
  };

  return (
    <div className="card example-card" style={{ borderLeftColor: example.color }}>
      <div className="card-preview">
        <div className="preview-placeholder" style={{ background: `linear-gradient(135deg, ${example.color}22, ${example.color}44)` }}>
          <div className="preview-icon">🎬</div>
        </div>
      </div>
      
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
            {example.name}
          </h3>
        )}
        
        <div className="card-meta">
          <span className="duration">{formatDuration(example.duration)}</span>
        </div>
      </div>

      <div className="card-actions">
        <button 
          className="action-btn" 
          onClick={onPreview}
          title="Preview"
        >
          📹
        </button>
        <button 
          className="action-btn" 
          onClick={() => onUpdate({ showInXR: !example.showInXR })}
          title={example.showInXR ? "Hide in XR" : "Show in XR"}
        >
          {example.showInXR ? '👁' : '🚫'}
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
          onClick={() => setIsEditing(true)}
          title="Rename"
        >
          ✏️
        </button>
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

export default ExampleCard;

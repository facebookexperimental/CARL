/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * Full-page preview mode for inspecting and trimming a single CARL example.
 * Renders a Babylon.js 3D scene that scrubs through the recorded hand-joint
 * data.  The sidebar panel lets users rename, recolor, set XR visibility,
 * and adjust start/end trim bounds.
 *
 * Note: `formatTime` is intentionally local — it uses a higher-precision
 * MM:SS.ss format distinct from the MM:SS `formatDuration` used on cards.
 *
 * @module components/PreviewMode
 */

import React, { useState, useEffect, useRef } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import { initializePreviewExperienceAsync } from 'carl-actionstudio-immersiveexperience';
import '../styles/PreviewMode.css';

function PreviewMode({ examples, onUpdateExample, onDeleteExample, carl }) {
  const { exampleId } = useParams();
  const navigate = useNavigate();
  const canvasRef = useRef(null);

  // --- State ---
  const [currentExample, setCurrentExample] = useState(null);
  const [isPlaying, setIsPlaying] = useState(false);
  const [currentTime, setCurrentTime] = useState(0);
  const [startTime, setStartTime] = useState(0);
  const [endTime, setEndTime] = useState(0);
  const [isEditingName, setIsEditingName] = useState(false);
  const [name, setName] = useState('');
  const [color, setColor] = useState('');
  const [searchQuery, setSearchQuery] = useState('');
  const [isDraggingStart, setIsDraggingStart] = useState(false);
  const [isDraggingEnd, setIsDraggingEnd] = useState(false);
  const timelineRef = useRef(null);
  const experienceHandleRef = useRef(null);

  // --- Effects ---
  useEffect(() => {
    if (exampleId) {
      const example = examples.find(ex => ex.id === Number.parseInt(exampleId));
      if (example) {
        setCurrentExample(example);
        setName(example.name);
        setColor(example.color);
        setStartTime(example.startTime);
        setEndTime(example.endTime);
        setCurrentTime(example.startTime);
      } else if (examples.length > 0) {
        navigate('/library', { replace: true });
      }
    } else if (examples.length > 0) {
      // Load first example if no ID specified
      navigate(`/preview/${examples[0].id}`, { replace: true });
    }
  }, [exampleId, examples, navigate]);

  // Effect A — initialize/teardown Babylon.js scene when example changes
  useEffect(() => {
    if (!canvasRef.current || !currentExample || !carl) return;
    let handle;
    const carlExample = carl.tryDeserializeExample(currentExample.bytes);
    if (!carlExample) return;
    initializePreviewExperienceAsync(canvasRef.current, carlExample)
      .then(h => {
        handle = h;
        experienceHandleRef.current = h;
        h.onTimeUpdate = (t) => setCurrentTime(t);
        h.onPlaybackEnd = () => setIsPlaying(false);
        h.setTime(currentExample.startTime);
      });
    return () => {
      experienceHandleRef.current = null;
      handle?.dispose();
      carlExample.dispose();
    };
  }, [currentExample, carl]);

  // Effect B — scrub to current time when not playing
  useEffect(() => {
    if (!isPlaying) {
      experienceHandleRef.current?.setTime(currentTime);
    }
  }, [currentTime, isPlaying]);

  // --- Handlers ---
  const handleExampleSelect = (example) => {
    navigate(`/preview/${example.id}`);
  };

  const handlePlayPause = () => {
    const newPlaying = !isPlaying;
    setIsPlaying(newPlaying);
    experienceHandleRef.current?.setPlaying(newPlaying);
  };

  const handleTimelineChange = (e) => {
    const newTime = parseFloat(e.target.value);
    setCurrentTime(newTime);
  };

  const handleStartTimeChange = (value) => {
    const newStart = Math.max(0, Math.min(value, endTime - 0.1));
    setStartTime(newStart);
    setCurrentTime(newStart);
    if (isPlaying) {
      setIsPlaying(false);
      experienceHandleRef.current?.setPlaying(false);
    }
    experienceHandleRef.current?.setTrimBounds(newStart, endTime);
  };

  const handleEndTimeChange = (value) => {
    const newEnd = Math.max(startTime + 0.1, Math.min(value, currentExample.duration));
    setEndTime(newEnd);
    setCurrentTime(newEnd);
    if (isPlaying) {
      setIsPlaying(false);
      experienceHandleRef.current?.setPlaying(false);
    }
    experienceHandleRef.current?.setTrimBounds(startTime, newEnd);
  };

  const handlePointerMove = (e) => {
    if (!timelineRef.current || !currentExample) return;
    
    const rect = timelineRef.current.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const percentage = Math.max(0.0, Math.min(1.0, x / rect.width));
    const time = percentage * currentExample.duration;
    
    if (isDraggingStart) {
      handleStartTimeChange(time);
    } else if (isDraggingEnd) {
      handleEndTimeChange(time);
    }
  };

  const handlePointerUp = () => {
    setIsDraggingStart(false);
    setIsDraggingEnd(false);
  };

  useEffect(() => {
    if (isDraggingStart || isDraggingEnd) {
      window.addEventListener('pointermove', handlePointerMove);
      window.addEventListener('pointerup', handlePointerUp);
      return () => {
        window.removeEventListener('pointermove', handlePointerMove);
        window.removeEventListener('pointerup', handlePointerUp);
      };
    }
  }, [isDraggingStart, isDraggingEnd, currentExample]);

  const handleSave = () => {
    if (currentExample) {
      onUpdateExample(currentExample.id, {
        name: name.trim(),
        color,
        startTime,
        endTime,
      });
    }
  };

  const handleDelete = () => {
    if (currentExample && window.confirm(`Delete "${currentExample.name}"?`)) {
      onDeleteExample(currentExample.id);
      navigate('/library');
    }
  };

  const handleNameSubmit = () => {
    if (name.trim() && name !== currentExample.name) {
      onUpdateExample(currentExample.id, { name: name.trim() });
    }
    setIsEditingName(false);
  };

  const handleToggleXR = () => {
    if (currentExample) {
      onUpdateExample(currentExample.id, { showInXR: !currentExample.showInXR });
    }
  };

  const formatTime = (seconds) => {
    const mins = Math.floor(seconds / 60);
    const secs = (seconds % 60).toFixed(2);
    return `${mins.toString().padStart(2, '0')}:${parseFloat(secs).toFixed(2).padStart(5, '0')}`;
  };

  // --- Render ---
  const filteredExamples = examples.filter(ex =>
    ex.name.toLowerCase().includes(searchQuery.toLowerCase())
  );

  if (!currentExample) {
    return (
      <div className="preview-mode">
        <div className="no-examples">
          <p>No examples available. Record some actions first!</p>
          <button onClick={() => navigate('/library')}>Go to Library</button>
        </div>
      </div>
    );
  }

  return (
    <div className="preview-mode">
      <aside className="examples-browser">
        <div className="browser-header">
          <h3>Examples</h3>
          <input
            type="text"
            className="browser-search"
            placeholder="Search..."
            value={searchQuery}
            onChange={(e) => setSearchQuery(e.target.value)}
          />
        </div>
        <div className="examples-list">
          {filteredExamples.map(example => (
            <div
              key={example.id}
              className={`example-item ${example.id === currentExample.id ? 'active' : ''}`}
              onClick={() => handleExampleSelect(example)}
            >
              <div className="item-color" style={{ backgroundColor: example.color }}></div>
              <div className="item-info">
                <span className="item-name">{example.name}</span>
                <span className="item-duration">{formatTime(example.duration)}</span>
              </div>
            </div>
          ))}
        </div>
      </aside>

      <main className="preview-main">
        <div className="preview-header">
          {isEditingName ? (
            <input
              type="text"
              className="preview-name-input"
              value={name}
              onChange={(e) => setName(e.target.value)}
              onBlur={handleNameSubmit}
              onKeyDown={(e) => {
                if (e.key === 'Enter') handleNameSubmit();
                if (e.key === 'Escape') {
                  setName(currentExample.name);
                  setIsEditingName(false);
                }
              }}
              autoFocus
            />
          ) : (
            <h2 onClick={() => setIsEditingName(true)}>{currentExample.name}</h2>
          )}
          <button className="close-btn" onClick={() => navigate('/library')}>✕</button>
        </div>

        <div className="preview-canvas-container">
          <canvas
            ref={canvasRef}
            className="preview-canvas"
            width={1200}
            height={600}
          />
        </div>

        <div className="playback-controls">
          <button className="play-btn" onClick={handlePlayPause}>
            {isPlaying ? '⏸' : '▶'}
          </button>
          <span className="time-display">{formatTime(currentTime)}</span>
          <span className="time-separator">/</span>
          <span className="time-display">{formatTime(currentExample.duration)}</span>
        </div>

        <div className="timeline-container">
          <div className="timeline-labels">
            <span>Start: {formatTime(startTime)}</span>
            <span className="current-time-label">{formatTime(currentTime)}</span>
            <span>End: {formatTime(endTime)}</span>
          </div>
          
          <div className="timeline" ref={timelineRef}>
            <input
              type="range"
              className="timeline-scrubber"
              min="0"
              max={currentExample.duration}
              step="0.01"
              value={currentTime}
              style={{ background: `linear-gradient(to right, #5B9FFF 0%, #5B9FFF ${(currentTime / currentExample.duration) * 100}%, #2a2a3e ${(currentTime / currentExample.duration) * 100}%, #2a2a3e 100%)` }}
              onChange={handleTimelineChange}
            />
            <div className="trim-handles">
              <div 
                className="trim-handle start-handle"
                style={{ left: `${(startTime / currentExample.duration) * 100}%` }}
                onPointerDown={() => setIsDraggingStart(true)}
              >
                <span className="handle-label">START</span>
                <div className="handle-grabber"></div>
              </div>
              <div 
                className="trim-handle end-handle"
                style={{ left: `${(endTime / currentExample.duration) * 100}%` }}
                onPointerDown={() => setIsDraggingEnd(true)}
              >
                <span className="handle-label">END</span>
                <div className="handle-grabber"></div>
              </div>
            </div>
          </div>
          
          <div className="timeline-timestamps">
            <input
              type="number"
              className="timestamp-input"
              value={startTime.toFixed(2)}
              onChange={(e) => handleStartTimeChange(parseFloat(e.target.value))}
              step="0.01"
            />
            <input
              type="number"
              className="timestamp-input"
              value={endTime.toFixed(2)}
              onChange={(e) => handleEndTimeChange(parseFloat(e.target.value))}
              step="0.01"
            />
          </div>
        </div>
      </main>

      <aside className="preview-sidebar">
        <h3>Metadata</h3>
        
        <div className="metadata-section">
          <label>Color</label>
          <input
            type="color"
            className="color-picker"
            value={color}
            onChange={(e) => setColor(e.target.value)}
          />
        </div>

        <div className="metadata-section">
          <label className="toggle-label">
            <input
              type="checkbox"
              checked={currentExample.showInXR}
              onChange={handleToggleXR}
            />
            <span>Show in XR Experience</span>
          </label>
        </div>

        <div className="metadata-section">
          <label>Start Time</label>
          <input
            type="number"
            className="time-input"
            value={startTime.toFixed(2)}
            onChange={(e) => handleStartTimeChange(parseFloat(e.target.value))}
            step="0.01"
          />
        </div>

        <div className="metadata-section">
          <label>End Time</label>
          <input
            type="number"
            className="time-input"
            value={endTime.toFixed(2)}
            onChange={(e) => handleEndTimeChange(parseFloat(e.target.value))}
            step="0.01"
          />
        </div>

        <button className="save-btn" onClick={handleSave}>Save Changes</button>
        <button className="delete-btn" onClick={handleDelete}>Delete Example</button>
      </aside>
    </div>
  );
}

export default PreviewMode;

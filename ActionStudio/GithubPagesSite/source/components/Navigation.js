/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * Top navigation bar with links to Library, Preview, and Definition Builder.
 *
 * @module components/Navigation
 */

import React from 'react';
import { NavLink, useLocation } from 'react-router-dom';
import '../styles/Navigation.css';

function Navigation() {
  const location = useLocation();

  const isActive = (path) => {
    return location.pathname.startsWith(path);
  };

  return (
    <nav className="navigation">
      <div className="nav-logo">
        <div className="logo-icon">CA</div>
        <div className="logo-text">
          <span className="logo-title">CARL</span>
          <span className="logo-subtitle">Action Studio</span>
        </div>
      </div>
      
      <div className="nav-tabs">
        <NavLink 
          to="/library" 
          className={`nav-tab ${isActive('/library') ? 'active' : ''}`}
        >
          <span className="nav-icon">📚</span>
          Library
        </NavLink>
        <NavLink 
          to="/preview" 
          className={`nav-tab ${isActive('/preview') ? 'active' : ''}`}
        >
          <span className="nav-icon">👁</span>
          Preview
        </NavLink>
        <NavLink 
          to="/definition-builder" 
          className={`nav-tab ${isActive('/definition-builder') ? 'active' : ''}`}
        >
          <span className="nav-icon">🔧</span>
          Definition Builder
        </NavLink>
      </div>

      <div className="nav-user">
        <div className="user-avatar">U</div>
      </div>
    </nav>
  );
}

export default Navigation;

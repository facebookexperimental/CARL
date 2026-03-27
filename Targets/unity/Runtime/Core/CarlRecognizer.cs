/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

using System;
using Carl.Native;

namespace Carl
{
    /// <summary>
    /// A recognizer that evaluates live input against a Definition and produces
    /// a confidence score (0.0–1.0+). Created asynchronously via
    /// <see cref="CarlSession.CreateRecognizerAsync"/>.
    /// </summary>
    public sealed class CarlRecognizer : IDisposable
    {
        ulong _handle;
        readonly CarlSession _session;

        internal CarlRecognizer(ulong handle, CarlSession session)
        {
            _handle = handle;
            _session = session;
        }

        internal ulong Handle
        {
            get
            {
                ThrowIfDisposed();
                return _handle;
            }
        }

        /// <summary>
        /// The current recognition score. Returns 0.0 when no match,
        /// approaching 1.0+ for strong matches.
        /// Thread-safe — can be polled from any thread.
        /// </summary>
        public double CurrentScore
        {
            get
            {
                ThrowIfDisposed();
                return CarlNative.carl_getCurrentScore(_handle);
            }
        }

        public void SetSensitivity(double sensitivity)
        {
            ThrowIfDisposed();
            CarlNative.carl_setSensitivity(_handle, sensitivity);
        }

        /// <summary>
        /// Returns a RecordingInspector for the canonical (idealized) recording
        /// that this recognizer was built from.
        /// The caller is responsible for disposing the returned inspector.
        /// </summary>
        public CarlRecordingInspector GetCanonicalRecordingInspector()
        {
            ThrowIfDisposed();
            ulong inspectorPtr = CarlNative.carl_getCanonicalRecordingInspector(_handle);
            return new CarlRecordingInspector(inspectorPtr, null);
        }

        void ThrowIfDisposed()
        {
            if (_handle == 0)
                throw new ObjectDisposedException(nameof(CarlRecognizer));
        }

        /// <summary>
        /// Disposes this recognizer. Deletion is scheduled asynchronously on
        /// the session's processing thread.
        /// </summary>
        public void Dispose()
        {
            if (_handle != 0)
            {
                _session.DisposeRecognizer(this);
                _handle = 0;
            }
        }

        internal void DisposeNative(ulong sessionHandle)
        {
            if (_handle != 0)
            {
                CarlNative.carl_disposeRecognizer(sessionHandle, _handle);
                _handle = 0;
            }
        }

        // Note: No finalizer. Recognizer disposal requires the session handle
        // (carl_disposeRecognizer takes both sessionPtr and recognizerPtr), so it
        // cannot be safely cleaned up from a finalizer where the session may already
        // be disposed. Recognizers are tracked in CarlSession._recognizers and
        // cleaned up when the session is disposed.
    }
}

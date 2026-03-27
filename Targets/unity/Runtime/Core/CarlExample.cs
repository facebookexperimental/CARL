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
    /// An example of an action occurrence within a Recording,
    /// defined by start and end timestamps.
    /// </summary>
    public sealed class CarlExample : IDisposable
    {
        ulong _handle;

        internal CarlExample(ulong handle)
        {
            _handle = handle;
        }

        internal ulong Handle
        {
            get
            {
                ThrowIfDisposed();
                return _handle;
            }
        }

        public static CarlExample Create(CarlRecording recording, double startTimestamp, double endTimestamp)
        {
            if (recording == null) throw new ArgumentNullException(nameof(recording));
            ulong handle = CarlNative.carl_createExample(recording.Handle, startTimestamp, endTimestamp);
            return new CarlExample(handle);
        }

        public static CarlExample CreateAutoTrimmed(CarlRecognizer recognizer, CarlRecording recording)
        {
            if (recognizer == null) throw new ArgumentNullException(nameof(recognizer));
            if (recording == null) throw new ArgumentNullException(nameof(recording));
            ulong handle = CarlNative.carl_createAutoTrimmedExample(recognizer.Handle, recording.Handle);
            return new CarlExample(handle);
        }

        /// <summary>
        /// Returns a new copy of the Recording associated with this Example.
        /// The caller is responsible for disposing the returned Recording.
        /// </summary>
        public CarlRecording GetRecording()
        {
            ThrowIfDisposed();
            ulong recordingPtr = CarlNative.carl_getRecording(_handle);
            return new CarlRecording(recordingPtr);
        }

        public double StartTimestamp
        {
            get
            {
                ThrowIfDisposed();
                return CarlNative.carl_getExampleStartTimestamp(_handle);
            }
        }

        public double EndTimestamp
        {
            get
            {
                ThrowIfDisposed();
                return CarlNative.carl_getExampleEndTimestamp(_handle);
            }
        }

        /// <summary>
        /// Loads an Example from a file. Returns null if the file could not be loaded.
        /// </summary>
        public static CarlExample LoadFromFile(string path)
        {
            ulong handle = CarlNative.carl_loadExampleFromFile(path);
            return handle != 0 ? new CarlExample(handle) : null;
        }

        public void SaveToFile(string path)
        {
            ThrowIfDisposed();
            CarlNative.carl_saveExampleToFile(_handle, path);
        }

        void ThrowIfDisposed()
        {
            if (_handle == 0)
                throw new ObjectDisposedException(nameof(CarlExample));
        }

        public void Dispose()
        {
            if (_handle != 0)
            {
                CarlNative.carl_disposeExample(_handle);
                _handle = 0;
            }
        }

        ~CarlExample()
        {
            if (_handle != 0)
                CarlNative.carl_disposeExample(_handle);
        }
    }
}

/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

using System;
using System.Runtime.InteropServices;
using Carl.Native;

namespace Carl
{
    /// <summary>
    /// A definition of an action, consisting of examples (and optional counterexamples)
    /// of the action being performed. Used to create a <see cref="CarlRecognizer"/>.
    /// </summary>
    public sealed class CarlDefinition : IDisposable
    {
        ulong _handle;

        public CarlDefinition(ActionType actionType)
        {
            _handle = CarlNative.carl_createDefinition((ulong)actionType);
        }

        internal CarlDefinition(ulong handle)
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

        public void AddExample(CarlRecording recording, double startTimestamp, double endTimestamp)
        {
            ThrowIfDisposed();
            if (recording == null) throw new ArgumentNullException(nameof(recording));
            CarlNative.carl_addExample(_handle, recording.Handle, startTimestamp, endTimestamp);
        }

        public void AddCounterexample(CarlRecording recording, double startTimestamp, double endTimestamp)
        {
            ThrowIfDisposed();
            if (recording == null) throw new ArgumentNullException(nameof(recording));
            CarlNative.carl_addCounterexample(_handle, recording.Handle, startTimestamp, endTimestamp);
        }

        public double DefaultSensitivity
        {
            get
            {
                ThrowIfDisposed();
                return CarlNative.carl_getDefaultSensitivity(_handle);
            }
            set
            {
                ThrowIfDisposed();
                CarlNative.carl_setDefaultSensitivity(_handle, value);
            }
        }

        public ulong ExamplesCount
        {
            get
            {
                ThrowIfDisposed();
                return CarlNative.carl_getExamplesCount(_handle);
            }
        }

        public ulong CounterexamplesCount
        {
            get
            {
                ThrowIfDisposed();
                return CarlNative.carl_getCounterexamplesCount(_handle);
            }
        }

        /// <summary>
        /// Returns a new copy of the Example at the given index.
        /// The caller is responsible for disposing the returned Example.
        /// </summary>
        public CarlExample GetExampleAt(int index)
        {
            ThrowIfDisposed();
            ulong handle = CarlNative.carl_getExampleAtIdx(_handle, (ulong)index);
            return new CarlExample(handle);
        }

        /// <summary>
        /// Returns a new copy of the Counterexample at the given index.
        /// The caller is responsible for disposing the returned Example.
        /// </summary>
        public CarlExample GetCounterexampleAt(int index)
        {
            ThrowIfDisposed();
            ulong handle = CarlNative.carl_getCounterexampleAtIdx(_handle, (ulong)index);
            return new CarlExample(handle);
        }

        public byte[] Serialize()
        {
            ThrowIfDisposed();
            ulong bytesPtr = CarlNative.carl_serializeDefinition(_handle);
            return CarlBytes.CopyAndFree(bytesPtr);
        }

        public static CarlDefinition Deserialize(byte[] data)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            GCHandle pinned = GCHandle.Alloc(data, GCHandleType.Pinned);
            try
            {
                ulong handle = CarlNative.carl_deserializeDefinition(pinned.AddrOfPinnedObject(), (ulong)data.Length);
                return new CarlDefinition(handle);
            }
            finally
            {
                pinned.Free();
            }
        }

        /// <summary>
        /// Loads a Definition from a file. Returns null if the file could not be loaded.
        /// </summary>
        public static CarlDefinition LoadFromFile(string path)
        {
            ulong handle = CarlNative.carl_loadDefinitionFromFile(path);
            return handle != 0 ? new CarlDefinition(handle) : null;
        }

        public void SaveToFile(string path)
        {
            ThrowIfDisposed();
            CarlNative.carl_saveDefinitionToFile(_handle, path);
        }

        void ThrowIfDisposed()
        {
            if (_handle == 0)
                throw new ObjectDisposedException(nameof(CarlDefinition));
        }

        public void Dispose()
        {
            if (_handle != 0)
            {
                CarlNative.carl_disposeDefinition(_handle);
                _handle = 0;
            }
        }

        ~CarlDefinition()
        {
            if (_handle != 0)
                CarlNative.carl_disposeDefinition(_handle);
        }
    }
}

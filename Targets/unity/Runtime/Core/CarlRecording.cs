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
    /// A completed recording of InputSamples.
    /// Can be serialized/deserialized and used to create Examples and Definitions.
    /// </summary>
    public sealed class CarlRecording : IDisposable
    {
        ulong _handle;

        internal CarlRecording(ulong handle)
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

        public byte[] Serialize()
        {
            ThrowIfDisposed();
            ulong bytesPtr = CarlNative.carl_serializeRecording(_handle);
            return CarlBytes.CopyAndFree(bytesPtr);
        }

        public static CarlRecording Deserialize(byte[] data)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            GCHandle pinned = GCHandle.Alloc(data, GCHandleType.Pinned);
            try
            {
                ulong handle = CarlNative.carl_deserializeRecording(pinned.AddrOfPinnedObject(), (ulong)data.Length);
                return new CarlRecording(handle);
            }
            finally
            {
                pinned.Free();
            }
        }

        /// <summary>
        /// Creates a RecordingInspector for scrubbing through this recording.
        /// The inspector is a non-owning view — this Recording must remain alive
        /// for the lifetime of the returned inspector.
        /// </summary>
        public CarlRecordingInspector GetInspector()
        {
            ThrowIfDisposed();
            ulong inspectorPtr = CarlNative.carl_getRecordingInspector(_handle);
            return new CarlRecordingInspector(inspectorPtr, this);
        }

        void ThrowIfDisposed()
        {
            if (_handle == 0)
                throw new ObjectDisposedException(nameof(CarlRecording));
        }

        public void Dispose()
        {
            if (_handle != 0)
            {
                CarlNative.carl_disposeRecording(_handle);
                _handle = 0;
            }
        }

        ~CarlRecording()
        {
            if (_handle != 0)
                CarlNative.carl_disposeRecording(_handle);
        }
    }
}

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
    /// Allows scrubbing through a <see cref="CarlRecording"/> by timestamp.
    /// This is a non-owning view: the parent Recording must remain alive.
    /// </summary>
    public sealed class CarlRecordingInspector : IDisposable
    {
        ulong _handle;
        // Hold a reference to the parent recording to prevent GC collection,
        // since the native inspector is a non-owning view.
        readonly CarlRecording _parentRecording;

        internal CarlRecordingInspector(ulong handle, CarlRecording parentRecording)
        {
            _handle = handle;
            _parentRecording = parentRecording;
        }

        public double StartTimestamp
        {
            get
            {
                ThrowIfDisposed();
                return CarlNative.carl_getStartTimestamp(_handle);
            }
        }

        public double EndTimestamp
        {
            get
            {
                ThrowIfDisposed();
                return CarlNative.carl_getEndTimestamp(_handle);
            }
        }

        /// <summary>
        /// Inspects the recording at the given timestamp.
        /// Returns the serialized InputSample bytes at that point in time.
        /// </summary>
        public byte[] Inspect(double timestamp)
        {
            ThrowIfDisposed();
            ulong bytesPtr = CarlNative.carl_inspect(_handle, timestamp);
            return CarlBytes.CopyAndFree(bytesPtr);
        }

        void ThrowIfDisposed()
        {
            if (_handle == 0)
                throw new ObjectDisposedException(nameof(CarlRecordingInspector));
        }

        public void Dispose()
        {
            if (_handle != 0)
            {
                CarlNative.carl_disposeRecordingInspector(_handle);
                _handle = 0;
            }
        }

        ~CarlRecordingInspector()
        {
            if (_handle != 0)
                CarlNative.carl_disposeRecordingInspector(_handle);
        }
    }
}

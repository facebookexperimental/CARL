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
    /// An in-progress recording that accumulates InputSamples over time.
    /// Call <see cref="Finish"/> to produce a completed <see cref="CarlRecording"/>.
    /// </summary>
    public sealed class CarlInProgressRecording : IDisposable
    {
        ulong _handle;

        public CarlInProgressRecording(ulong maxSeconds = ulong.MaxValue)
        {
            _handle = CarlNative.carl_startRecording(maxSeconds);
        }

        public void AddSample(ref CarlInputSampleInterop sample)
        {
            ThrowIfDisposed();
            CarlNative.carl_recordObjectInputSample(_handle, ref sample);
        }

        /// <summary>
        /// Finishes recording and returns the completed Recording.
        /// This instance becomes disposed after this call (ownership is transferred).
        /// </summary>
        public CarlRecording Finish()
        {
            ThrowIfDisposed();
            ulong recordingPtr = CarlNative.carl_finishRecording(_handle);
            _handle = 0;
            return new CarlRecording(recordingPtr);
        }

        void ThrowIfDisposed()
        {
            if (_handle == 0)
                throw new ObjectDisposedException(nameof(CarlInProgressRecording));
        }

        public void Dispose()
        {
            if (_handle != 0)
            {
#if UNITY_5_3_OR_NEWER
                UnityEngine.Debug.LogWarning(
                    "[CARL] CarlInProgressRecording disposed without calling Finish(). " +
                    "The native recording buffer has been leaked. Always call Finish() to " +
                    "produce a CarlRecording before disposing.");
#endif
                _handle = 0;
            }
        }
    }
}

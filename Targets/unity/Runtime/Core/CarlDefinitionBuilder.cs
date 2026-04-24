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
    /// Holds the results of <see cref="CarlDefinitionBuilder.CreateExamplesFromRecordings"/>,
    /// providing access to auto-trimmed <see cref="CarlExample"/>s discovered from
    /// multiple raw recordings.
    /// </summary>
    public sealed class CarlBuiltExamples : IDisposable
    {
        ulong _handle;

        CarlBuiltExamples(ulong handle)
        {
            _handle = handle;
        }

        /// <summary>
        /// The number of auto-trimmed examples in this result.
        /// </summary>
        public int Count
        {
            get
            {
                ThrowIfDisposed();
                return (int)CarlNative.carl_getBuiltExamplesCount(_handle);
            }
        }

        /// <summary>
        /// Returns a new copy of the Example at the given index.
        /// The caller is responsible for disposing the returned Example.
        /// </summary>
        public CarlExample GetExampleAt(int index)
        {
            ThrowIfDisposed();
            ulong handle = CarlNative.carl_getBuiltExampleAtIdx(_handle, (ulong)index);
            return new CarlExample(handle);
        }

        void ThrowIfDisposed()
        {
            if (_handle == 0)
                throw new ObjectDisposedException(nameof(CarlBuiltExamples));
        }

        public void Dispose()
        {
            if (_handle != 0)
            {
                CarlNative.carl_disposeBuiltExamples(_handle);
                _handle = 0;
            }
        }

        ~CarlBuiltExamples()
        {
            if (_handle != 0)
                CarlNative.carl_disposeBuiltExamples(_handle);
        }

        /// <summary>
        /// Provides the internal handle for use by <see cref="CarlDefinitionBuilder"/>.
        /// </summary>
        internal static CarlBuiltExamples FromHandle(ulong handle)
        {
            return new CarlBuiltExamples(handle);
        }
    }

    /// <summary>
    /// Builds auto-trimmed Examples from multiple raw Recordings by
    /// discovering the common action across all recordings. No
    /// pre-existing Definition or Recognizer is required.
    /// </summary>
    public static class CarlDefinitionBuilder
    {
        /// <summary>
        /// Analyzes multiple recordings to find and extract aligned Examples
        /// of a common action. The recordings should each contain the action
        /// somewhere within them, potentially with leading and trailing excess.
        /// </summary>
        /// <param name="actionType">The type of action being demonstrated.</param>
        /// <param name="recordings">Two or more recordings, each containing the action.</param>
        /// <param name="expectedDuration">
        /// Optional hint for the expected duration of the action in seconds.
        /// Pass 0 for unconstrained.
        /// </param>
        /// <returns>
        /// A <see cref="CarlBuiltExamples"/> containing one trimmed Example per
        /// input Recording. The caller is responsible for disposing the result.
        /// </returns>
        public static CarlBuiltExamples CreateExamplesFromRecordings(
            ActionType actionType,
            CarlRecording[] recordings,
            double expectedDuration = 0.0)
        {
            if (recordings == null) throw new ArgumentNullException(nameof(recordings));

            ulong[] handles = new ulong[recordings.Length];
            for (int i = 0; i < recordings.Length; i++)
            {
                if (recordings[i] == null) throw new ArgumentNullException($"recordings[{i}]");
                handles[i] = recordings[i].Handle;
            }

            GCHandle pinned = GCHandle.Alloc(handles, GCHandleType.Pinned);
            try
            {
                ulong resultPtr = CarlNative.carl_createExamplesFromRecordings(
                    (ulong)actionType,
                    pinned.AddrOfPinnedObject(),
                    (ulong)recordings.Length,
                    expectedDuration);
                return CarlBuiltExamples.FromHandle(resultPtr);
            }
            finally
            {
                pinned.Free();
            }
        }
    }
}

/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

using System;
using System.Runtime.InteropServices;

namespace Carl.Native
{
    internal static class CarlBytes
    {
        /// <summary>
        /// Copies bytes from a native CARL byte buffer and frees the native allocation.
        /// Uses the two-call pattern: first call with size=0 to query length,
        /// second call to copy data and free the native buffer.
        /// </summary>
        public static byte[] CopyAndFree(ulong bytesPtr)
        {
            if (bytesPtr == 0)
                throw new ArgumentException("Invalid bytes pointer.", nameof(bytesPtr));

            ulong size = CarlNative.carl_getBytes(bytesPtr, IntPtr.Zero, 0);
            if (size == 0)
                return Array.Empty<byte>();

            byte[] result = new byte[size];
            GCHandle handle = GCHandle.Alloc(result, GCHandleType.Pinned);
            try
            {
                CarlNative.carl_getBytes(bytesPtr, handle.AddrOfPinnedObject(), size);
            }
            finally
            {
                handle.Free();
            }
            return result;
        }
    }
}

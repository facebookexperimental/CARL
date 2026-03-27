/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using Carl.Native;
using AOT;

namespace Carl
{
    /// <summary>
    /// The central CARL session that manages input processing and recognizers.
    /// Feed InputSamples each frame and call <see cref="TickCallbacks"/> from the main thread.
    /// </summary>
    public sealed class CarlSession : IDisposable
    {
        ulong _handle;
        readonly List<CarlRecognizer> _recognizers = new List<CarlRecognizer>();

        // Static callback routing for IL2CPP/AOT compatibility.
        // ConcurrentDictionary used for logger callbacks since the native logger
        // may fire from a background thread while the main thread mutates the dictionary.
        static readonly ConcurrentDictionary<ulong, Action<string>> s_loggerCallbacks = new ConcurrentDictionary<ulong, Action<string>>();
        static readonly Dictionary<ulong, (CarlSession session, CarlDefinition definition, Action<CarlRecognizer> callback)> s_recognizerCallbacks =
            new Dictionary<ulong, (CarlSession, CarlDefinition, Action<CarlRecognizer>)>();
        static ulong s_nextRequestId = 1;

        // Must hold references to prevent GC collection of delegates passed to native code.
        static readonly CarlLoggerCallback s_loggerTrampoline = LoggerTrampoline;
        static readonly CarlRecognizerCreatedCallback s_recognizerTrampoline = RecognizerCreatedTrampoline;

        public CarlSession(bool singleThreaded = false)
        {
            _handle = singleThreaded
                ? CarlNative.carl_createSingleThreadedSession()
                : CarlNative.carl_createSession();
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
        /// Sets a logger callback for diagnostic messages from the CARL engine.
        /// </summary>
        public void SetLogger(Action<string> logger)
        {
            ThrowIfDisposed();
            if (logger != null)
            {
                s_loggerCallbacks[_handle] = logger;
                CarlNative.carl_setSessionLogger(_handle, s_loggerTrampoline);
            }
            else
            {
                s_loggerCallbacks.TryRemove(_handle, out _);
                CarlNative.carl_setSessionLogger(_handle, null);
            }
        }

        /// <summary>
        /// Processes pending callbacks on the main thread.
        /// Must be called every frame (typically from MonoBehaviour.Update).
        /// </summary>
        public void TickCallbacks()
        {
            ThrowIfDisposed();
            CarlNative.carl_tickCallbacks(_handle);
        }

        /// <summary>
        /// Feeds an input sample to the session for processing.
        /// </summary>
        public void AddInputSample(ref CarlInputSampleInterop sample)
        {
            ThrowIfDisposed();
            CarlNative.carl_addInputSample(_handle, ref sample);
        }

        /// <summary>
        /// Creates a recognizer asynchronously. The callback fires on the main thread
        /// after <see cref="TickCallbacks"/> is called.
        /// The definition is kept alive (preventing GC/finalization) until the
        /// native async operation completes.
        /// </summary>
        public void CreateRecognizerAsync(CarlDefinition definition, Action<CarlRecognizer> callback)
        {
            ThrowIfDisposed();
            if (definition == null) throw new ArgumentNullException(nameof(definition));
            if (callback == null) throw new ArgumentNullException(nameof(callback));

            ulong requestId = s_nextRequestId++;
            // Store the definition reference to prevent GC from finalizing it
            // while the native async operation is still using the native handle.
            s_recognizerCallbacks[requestId] = (this, definition, callback);
            CarlNative.carl_createRecognizerAsync(_handle, definition.Handle, requestId, s_recognizerTrampoline);
        }

        internal void DisposeRecognizer(CarlRecognizer recognizer)
        {
            if (_handle != 0)
            {
                _recognizers.Remove(recognizer);
                recognizer.DisposeNative(_handle);
            }
        }

        void ThrowIfDisposed()
        {
            if (_handle == 0)
                throw new ObjectDisposedException(nameof(CarlSession));
        }

        public void Dispose()
        {
            if (_handle != 0)
            {
                // Dispose all owned recognizers first.
                for (int i = _recognizers.Count - 1; i >= 0; i--)
                {
                    _recognizers[i].DisposeNative(_handle);
                }
                _recognizers.Clear();

                // Remove pending async callbacks that reference this session
                // to prevent stale references and use-after-dispose.
                var staleKeys = new List<ulong>();
                foreach (var kvp in s_recognizerCallbacks)
                {
                    if (kvp.Value.session == this)
                        staleKeys.Add(kvp.Key);
                }
                foreach (var key in staleKeys)
                    s_recognizerCallbacks.Remove(key);

                s_loggerCallbacks.TryRemove(_handle, out _);

                CarlNative.carl_disposeSession(_handle);
                _handle = 0;
            }
        }

        ~CarlSession()
        {
            if (_handle != 0)
            {
                s_loggerCallbacks.TryRemove(_handle, out _);
                CarlNative.carl_disposeSession(_handle);
            }
        }

        // --- Static trampolines for native callbacks (IL2CPP/AOT safe) ---

        [MonoPInvokeCallback(typeof(CarlLoggerCallback))]
        static void LoggerTrampoline(IntPtr messagePtr)
        {
            // We cannot identify which session this came from via the callback signature,
            // so we iterate. In practice there is typically one session.
            string message = Marshal.PtrToStringAnsi(messagePtr);
            foreach (var kvp in s_loggerCallbacks)
            {
                try
                {
                    kvp.Value?.Invoke(message);
                }
                catch (Exception)
                {
                    // Swallow exceptions in user callbacks to avoid crashing native code.
                }
            }
        }

        [MonoPInvokeCallback(typeof(CarlRecognizerCreatedCallback))]
        static void RecognizerCreatedTrampoline(ulong requestId, ulong recognizerPtr)
        {
            if (s_recognizerCallbacks.TryGetValue(requestId, out var entry))
            {
                s_recognizerCallbacks.Remove(requestId);
                // entry.definition reference is released here, allowing GC if no other refs exist.
                var recognizer = new CarlRecognizer(recognizerPtr, entry.session);
                entry.session._recognizers.Add(recognizer);
                try
                {
                    entry.callback?.Invoke(recognizer);
                }
                catch (Exception e)
                {
                    UnityEngine.Debug.LogError($"[CARL] Exception in recognizer creation callback: {e.GetType().Name}: {e.Message}");
                }
            }
            else
            {
                UnityEngine.Debug.LogWarning($"[CARL] RecognizerCreatedTrampoline: unknown requestId={requestId}. Callback may have been cleaned up.");
            }
        }
    }
}

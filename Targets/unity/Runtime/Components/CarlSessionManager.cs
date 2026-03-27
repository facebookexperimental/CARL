/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

using UnityEngine;

namespace Carl
{
    /// <summary>
    /// Singleton MonoBehaviour that manages the CARL session lifecycle.
    /// Place on a GameObject in your scene to initialize CARL.
    /// Calls TickCallbacks() every frame to dispatch async results.
    /// </summary>
    [DefaultExecutionOrder(-100)]
    [AddComponentMenu("CARL/Session Manager")]
    public class CarlSessionManager : MonoBehaviour
    {
        [Tooltip("Enable diagnostic logging from the CARL engine.")]
        [SerializeField] bool enableLogging = true;

        [Tooltip("Use single-threaded mode. Disables background processing — " +
                 "all work happens inline. Useful for debugging.")]
        [SerializeField] bool singleThreaded;

        public static CarlSessionManager Instance { get; private set; }
        public CarlSession Session { get; private set; }

        void Awake()
        {
            if (Instance != null && Instance != this)
            {
                Debug.LogWarning("[CARL] Duplicate CarlSessionManager detected. Destroying this instance.");
                Destroy(gameObject);
                return;
            }

            Instance = this;
            Session = new CarlSession(singleThreaded);

            if (enableLogging)
            {
                Session.SetLogger(message => Debug.Log($"[CARL] {message}"));
            }
        }

        void Update()
        {
            Session?.TickCallbacks();
        }

        void OnDestroy()
        {
            if (Instance == this)
            {
                Session?.Dispose();
                Session = null;
                Instance = null;
            }
        }
    }
}

/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

using Carl.Native;
using UnityEngine;
using UnityEngine.Events;

namespace Carl
{
    /// <summary>
    /// Records InputSamples at runtime to create Examples for gesture definitions.
    /// Requires a <see cref="CarlXRInputProvider"/> in the scene to capture input.
    /// </summary>
    [AddComponentMenu("CARL/Gesture Recorder")]
    public class CarlGestureRecorder : MonoBehaviour
    {
        [Header("Events")]
        public UnityEvent OnRecordingStarted;
        public UnityEvent<CarlRecording> OnRecordingFinished;

        CarlInProgressRecording _recording;
        CarlXRInputProvider _inputProvider;
        double _recordingStartTime;

        /// <summary>
        /// Whether a recording is currently in progress.
        /// </summary>
        public bool IsRecording => _recording != null;

        /// <summary>
        /// The time elapsed since recording started, in seconds.
        /// </summary>
        public double ElapsedTime => IsRecording ? Time.realtimeSinceStartupAsDouble - _recordingStartTime : 0;

        void Awake()
        {
            _inputProvider = FindObjectOfType<CarlXRInputProvider>();
        }

        /// <summary>
        /// Begins recording InputSamples.
        /// </summary>
        /// <param name="maxSeconds">Maximum recording duration in seconds. Default is unlimited.</param>
        public void StartRecording(ulong maxSeconds = ulong.MaxValue)
        {
            if (_recording != null)
            {
                Debug.LogWarning("[CARL] Recording already in progress.");
                return;
            }

            _recording = new CarlInProgressRecording(maxSeconds);
            _recordingStartTime = Time.realtimeSinceStartupAsDouble;
            OnRecordingStarted?.Invoke();
        }

        /// <summary>
        /// Stops recording and returns the completed Recording.
        /// </summary>
        public CarlRecording StopRecording()
        {
            if (_recording == null)
            {
                Debug.LogWarning("[CARL] No recording in progress.");
                return null;
            }

            CarlRecording recording = _recording.Finish();
            _recording = null;
            OnRecordingFinished?.Invoke(recording);
            return recording;
        }

        /// <summary>
        /// Creates an Example from a recording with the given start and end timestamps.
        /// </summary>
        public CarlExample CreateExample(CarlRecording recording, double startTimestamp, double endTimestamp)
        {
            return CarlExample.Create(recording, startTimestamp, endTimestamp);
        }

        void LateUpdate()
        {
            if (_recording == null || _inputProvider == null)
                return;

            ref var sample = ref _inputProvider.CurrentSample;
            _recording.AddSample(ref sample);
        }

        void OnDestroy()
        {
            _recording?.Dispose();
            _recording = null;
        }
    }
}

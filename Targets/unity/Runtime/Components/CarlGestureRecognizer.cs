/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

using UnityEngine;
using UnityEngine.Events;

namespace Carl
{
    /// <summary>
    /// Wraps a CARL recognizer with threshold-based gesture detection.
    /// Fires UnityEvents on gesture detection and loss with hysteresis
    /// to prevent rapid toggling near the threshold boundary.
    /// </summary>
    [AddComponentMenu("CARL/Gesture Recognizer")]
    public class CarlGestureRecognizer : MonoBehaviour
    {
        [Header("Definition")]
        [Tooltip("A CarlDefinitionAsset containing the serialized gesture definition.")]
        [SerializeField] CarlDefinitionAsset definitionAsset;

        [Tooltip("Path to a .carl definition file. Used if definitionAsset is not set.")]
        [SerializeField] string definitionFilePath;

        [Header("Recognition")]
        [Tooltip("Score above which the gesture is considered detected.")]
        [SerializeField] float activationThreshold = 0.8f;

        [Tooltip("Score below which a detected gesture is considered lost. " +
                 "Should be lower than activationThreshold for hysteresis.")]
        [SerializeField] float deactivationThreshold = 0.6f;

        [Tooltip("Override the definition's default sensitivity. " +
                 "Set to -1 to use the definition's built-in sensitivity.")]
        [SerializeField] double sensitivity = -1;

        [Header("Events")]
        public UnityEvent<float> OnScoreChanged;
        public UnityEvent OnGestureDetected;
        public UnityEvent OnGestureLost;

        CarlRecognizer _recognizer;
        CarlDefinition _definition;
        bool _isRecognized;
        float _lastScore;

        /// <summary>
        /// Whether the gesture is currently detected.
        /// </summary>
        public bool IsRecognized => _isRecognized;

        /// <summary>
        /// The most recent recognition score (0.0–1.0+).
        /// </summary>
        public float LastScore => _lastScore;

        void Start()
        {
            if (CarlSessionManager.Instance == null)
            {
                Debug.LogError("[CARL] CarlGestureRecognizer requires a CarlSessionManager in the scene.");
                enabled = false;
                return;
            }

            CarlDefinition definition = null;

            if (definitionAsset != null)
            {
                definition = definitionAsset.Load();
            }
            else if (!string.IsNullOrEmpty(definitionFilePath))
            {
                definition = CarlDefinition.LoadFromFile(definitionFilePath);
            }

            if (definition == null)
            {
                Debug.LogError("[CARL] CarlGestureRecognizer: No valid definition provided.");
                enabled = false;
                return;
            }

            CarlSessionManager.Instance.Session.CreateRecognizerAsync(definition, recognizer =>
            {
                _recognizer = recognizer;
                if (sensitivity >= 0)
                {
                    _recognizer.SetSensitivity(sensitivity);
                }
            });

            _definition = definition;
        }

        void Update()
        {
            if (_recognizer == null)
                return;

            float score = (float)_recognizer.CurrentScore;

            if (score != _lastScore)
            {
                _lastScore = score;
                OnScoreChanged?.Invoke(score);
            }

            if (!_isRecognized && score >= activationThreshold)
            {
                _isRecognized = true;
                OnGestureDetected?.Invoke();
            }
            else if (_isRecognized && score < deactivationThreshold)
            {
                _isRecognized = false;
                OnGestureLost?.Invoke();
            }
        }

        void OnDestroy()
        {
            _recognizer?.Dispose();
            _recognizer = null;
            _definition?.Dispose();
            _definition = null;
        }
    }
}

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

            // No definition source configured via Inspector — the caller will
            // use one of the Initialize methods to set up the recognizer at runtime.
            if (definition == null)
                return;

            InitializeInternal(definition);
        }

        /// <summary>
        /// Initializes (or re-initializes) the recognizer at runtime with the given
        /// <see cref="CarlDefinition"/>. This component takes ownership of the
        /// definition and will dispose it in <see cref="OnDestroy"/>.
        /// Any previously active recognizer is disposed first.
        /// </summary>
        public void Initialize(CarlDefinition definition)
        {
            if (definition == null)
            {
                Debug.LogError("[CARL] CarlGestureRecognizer.Initialize: definition is null.");
                return;
            }

            if (CarlSessionManager.Instance == null)
            {
                Debug.LogError("[CARL] CarlGestureRecognizer requires a CarlSessionManager in the scene.");
                return;
            }

            DisposeRecognizer();
            InitializeInternal(definition);
        }

        /// <summary>
        /// Initializes (or re-initializes) the recognizer at runtime from a
        /// <see cref="CarlDefinitionAsset"/>. A new <see cref="CarlDefinition"/>
        /// is deserialized from the asset; the asset itself is not retained.
        /// Any previously active recognizer is disposed first.
        /// </summary>
        public void Initialize(CarlDefinitionAsset asset)
        {
            if (asset == null || !asset.HasData)
            {
                Debug.LogError("[CARL] CarlGestureRecognizer.Initialize: asset is null or has no data.");
                return;
            }

            Initialize(asset.Load());
        }

        /// <summary>
        /// Creates a new definition from raw recordings using
        /// <see cref="CarlDefinitionBuilder"/>, then initializes the recognizer
        /// with it. This is the all-in-one path for runtime gesture training:
        /// record → build definition → create recognizer.
        /// Any previously active recognizer is disposed first.
        /// </summary>
        /// <param name="actionType">The type of action demonstrated in the recordings.</param>
        /// <param name="recordings">Two or more recordings, each containing the action.</param>
        /// <param name="expectedDuration">
        /// Optional hint for the expected action duration in seconds. Pass 0 for unconstrained.
        /// </param>
        /// <returns>
        /// The created <see cref="CarlDefinition"/>, which is also owned by this component.
        /// Useful for persisting the definition via
        /// <see cref="CarlDefinitionAsset.SetFromDefinition"/> or
        /// <see cref="CarlDefinition.SaveToFile"/>.
        /// </returns>
        public CarlDefinition InitializeFromRecordings(
            ActionType actionType,
            CarlRecording[] recordings,
            double expectedDuration = 0.0,
            double defaultSensitivity = 5.0)
        {
            if (CarlSessionManager.Instance == null)
            {
                Debug.LogError("[CARL] CarlGestureRecognizer requires a CarlSessionManager in the scene.");
                return null;
            }

            var definition = CarlDefinition.CreateFromRecordings(actionType, recordings, expectedDuration);
            definition.DefaultSensitivity = defaultSensitivity;
            DisposeRecognizer();
            InitializeInternal(definition);
            return definition;
        }

        void InitializeInternal(CarlDefinition definition)
        {
            _definition = definition;

            CarlSessionManager.Instance.Session.CreateRecognizerAsync(definition, recognizer =>
            {
                _recognizer = recognizer;
                if (sensitivity >= 0)
                {
                    _recognizer.SetSensitivity(sensitivity);
                }
            });
        }

        void DisposeRecognizer()
        {
            _recognizer?.Dispose();
            _recognizer = null;
            _definition?.Dispose();
            _definition = null;
            _isRecognized = false;
            _lastScore = 0f;
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
            DisposeRecognizer();
        }
    }
}

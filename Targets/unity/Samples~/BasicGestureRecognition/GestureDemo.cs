/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

using Carl;
using Carl.Native;
using UnityEngine;

/// <summary>
/// Sample script demonstrating basic CARL gesture recognition setup.
/// Attach to any GameObject in a scene that has CarlSessionManager and CarlXRInputProvider.
/// </summary>
public class GestureDemo : MonoBehaviour
{
    [Header("Definition")]
    [Tooltip("Assign a CarlDefinitionAsset with a gesture definition.")]
    [SerializeField] CarlDefinitionAsset definitionAsset;

    [Tooltip("Or specify a file path to a .carl definition file.")]
    [SerializeField] string definitionFilePath;

    [Header("Recognition")]
    [SerializeField] float activationThreshold = 0.8f;

    CarlRecognizer _recognizer;
    CarlDefinition _definition;
    bool _isRecognized;

    void Start()
    {
        if (CarlSessionManager.Instance == null)
        {
            Debug.LogError("[GestureDemo] No CarlSessionManager found in scene. " +
                           "Add one to a GameObject before using CARL.");
            enabled = false;
            return;
        }

        // Load the definition from asset or file.
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
            Debug.LogError("[GestureDemo] No valid definition provided.");
            enabled = false;
            return;
        }

        Debug.Log($"[GestureDemo] Definition loaded with {definition.ExamplesCount} examples. " +
                  $"Default sensitivity: {definition.DefaultSensitivity:F2}");

        // Create the recognizer asynchronously.
        CarlSessionManager.Instance.Session.CreateRecognizerAsync(definition, recognizer =>
        {
            _recognizer = recognizer;
            Debug.Log("[GestureDemo] Recognizer created and ready.");
        });

        // Keep definition alive until async creation completes.
        _definition = definition;
    }

    void Update()
    {
        if (_recognizer == null)
            return;

        float score = (float)_recognizer.CurrentScore;

        if (!_isRecognized && score >= activationThreshold)
        {
            _isRecognized = true;
            Debug.Log($"[GestureDemo] Gesture DETECTED! Score: {score:F3}");
        }
        else if (_isRecognized && score < activationThreshold * 0.75f)
        {
            _isRecognized = false;
            Debug.Log($"[GestureDemo] Gesture lost. Score: {score:F3}");
        }
    }

    void OnDestroy()
    {
        _recognizer?.Dispose();
        _recognizer = null;
        _definition?.Dispose();
        _definition = null;
    }

    void OnGUI()
    {
        if (_recognizer == null)
        {
            GUI.Label(new Rect(10, 10, 400, 30), "[CARL] Initializing recognizer...");
            return;
        }

        float score = (float)_recognizer.CurrentScore;
        string status = _isRecognized ? "DETECTED" : "---";
        GUI.Label(new Rect(10, 10, 400, 30), $"[CARL] Score: {score:F3} | Status: {status}");
    }
}

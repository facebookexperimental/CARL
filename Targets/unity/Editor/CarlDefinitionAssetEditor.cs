/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

using UnityEditor;
using UnityEngine;

namespace Carl.Editor
{
    [CustomEditor(typeof(CarlDefinitionAsset))]
    public class CarlDefinitionAssetEditor : UnityEditor.Editor
    {
        public override void OnInspectorGUI()
        {
            var asset = (CarlDefinitionAsset)target;

            EditorGUILayout.LabelField("Action Type", asset.ActionType.ToString());
            EditorGUILayout.LabelField("Description", asset.Description);
            EditorGUILayout.Space();

            if (asset.HasData)
            {
                EditorGUILayout.HelpBox(
                    $"Definition data loaded ({asset.SerializedDefinition.Length:N0} bytes).",
                    MessageType.Info);
            }
            else
            {
                EditorGUILayout.HelpBox(
                    "No definition data. Use 'Import from File' to load a .carl definition.",
                    MessageType.Warning);
            }

            EditorGUILayout.Space();

            if (GUILayout.Button("Import from File"))
            {
                string path = EditorUtility.OpenFilePanel("Import CARL Definition", "", "carl");
                if (!string.IsNullOrEmpty(path))
                {
                    if (asset.ImportFromFile(path))
                    {
                        EditorUtility.SetDirty(asset);
                        AssetDatabase.SaveAssets();
                        Debug.Log($"[CARL] Imported definition from {path}");
                    }
                    else
                    {
                        EditorUtility.DisplayDialog("Import Failed",
                            $"Could not load CARL definition from:\n{path}", "OK");
                    }
                }
            }

            if (asset.HasData && GUILayout.Button("Export to File"))
            {
                string path = EditorUtility.SaveFilePanel("Export CARL Definition", "", "definition", "carl");
                if (!string.IsNullOrEmpty(path))
                {
                    using (var def = asset.Load())
                    {
                        def?.SaveToFile(path);
                        Debug.Log($"[CARL] Exported definition to {path}");
                    }
                }
            }

            EditorGUILayout.Space();
            DrawDefaultInspector();
        }
    }
}

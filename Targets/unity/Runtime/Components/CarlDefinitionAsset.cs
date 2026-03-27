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
    /// ScriptableObject that wraps a serialized CARL gesture definition,
    /// allowing it to be stored as a Unity asset and assigned in the Inspector.
    /// </summary>
    [CreateAssetMenu(fileName = "NewCarlDefinition", menuName = "CARL/Definition Asset")]
    public class CarlDefinitionAsset : ScriptableObject
    {
        [Tooltip("The serialized CARL definition bytes.")]
        [HideInInspector]
        [SerializeField] byte[] serializedDefinition;

        [Tooltip("The action type of this definition.")]
        [SerializeField] ActionType actionType;

        [Tooltip("A human-readable description of this gesture.")]
        [TextArea]
        [SerializeField] string description;

        /// <summary>The action type of the gesture definition.</summary>
        public ActionType ActionType => actionType;

        /// <summary>A human-readable description of this gesture.</summary>
        public string Description => description;

        /// <summary>Whether this asset contains valid serialized definition data.</summary>
        public bool HasData => serializedDefinition != null && serializedDefinition.Length > 0;

        /// <summary>The raw serialized definition bytes.</summary>
        public byte[] SerializedDefinition
        {
            get => serializedDefinition;
            set => serializedDefinition = value;
        }

        /// <summary>
        /// Deserializes and returns a new CarlDefinition instance.
        /// The caller is responsible for disposing the returned Definition.
        /// Returns null if no serialized data is present.
        /// </summary>
        public CarlDefinition Load()
        {
            if (!HasData)
                return null;
            return CarlDefinition.Deserialize(serializedDefinition);
        }

        /// <summary>
        /// Imports a definition from a file, storing its serialized bytes in this asset.
        /// </summary>
        public bool ImportFromFile(string path)
        {
            using (var def = CarlDefinition.LoadFromFile(path))
            {
                if (def == null) return false;
                serializedDefinition = def.Serialize();
                return true;
            }
        }

        /// <summary>
        /// Sets the definition data from a CarlDefinition instance.
        /// </summary>
        public void SetFromDefinition(CarlDefinition definition, ActionType type, string desc = null)
        {
            serializedDefinition = definition.Serialize();
            actionType = type;
            if (desc != null) description = desc;
        }
    }
}

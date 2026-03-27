/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

using UnityEditor;
using UnityEditor.Build;
using UnityEditor.Build.Reporting;
using UnityEngine;

namespace Carl.Editor
{
    /// <summary>
    /// Validates that the native CARL plugin binary exists for the target platform before building.
    /// </summary>
    public class CarlBuildPostprocessor : IPreprocessBuildWithReport
    {
        public int callbackOrder => 0;

        public void OnPreprocessBuild(BuildReport report)
        {
            string pluginPath = null;

            switch (report.summary.platform)
            {
                case BuildTarget.StandaloneWindows64:
                    pluginPath = "Packages/com.meta.carl/Runtime/Plugins/x86_64/carl.dll";
                    break;
                case BuildTarget.Android:
                    pluginPath = "Packages/com.meta.carl/Runtime/Plugins/Android/arm64-v8a/libcarl.so";
                    break;
            }

            if (pluginPath != null)
            {
                var importer = AssetImporter.GetAtPath(pluginPath);
                if (importer == null)
                {
                    Debug.LogWarning(
                        $"[CARL] Native plugin not found at '{pluginPath}'. " +
                        $"Build may fail at runtime. See the CARL README for build instructions.");
                }
            }
        }
    }
}

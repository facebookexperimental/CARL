/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * Public surface of the ImmersiveExperience package.
 *
 * Re-exports the CARL interface types, the main XR experience entry point, and the
 * standalone 2D preview experience so that GithubPagesSite can import them as a single module.
 */
export * from "./carlInterfaces";
export * from "./main";
export * from "./previewExperience";

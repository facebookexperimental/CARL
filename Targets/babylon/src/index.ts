/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * Public surface of `@meta-experimental/carl-babylon`.
 *
 * Bundles everything a Babylon.js project needs to integrate CARL:
 *  - WASM runtime loader (`initializeCarl`) and the shipped `carl.js` / `carl.wasm` assets.
 *  - Concrete `ICarl` implementation (`CarlIntegration` and friends).
 *  - The `ICarl*` interface contract.
 *  - Reusable Babylon.js glue: input-sample helpers, recognition graph, input puppet,
 *    and a standalone 2D preview experience.
 */
export * from "./carlInterfaces";
export * from "./wasmLoader";
export * from "./carlIntegration";
export * from "./utils";
export * from "./recognitionGraph";
export * from "./inputPuppet";
export * from "./previewExperience";

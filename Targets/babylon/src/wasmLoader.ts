/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

/**
 * Loader for the CARL Emscripten WASM runtime.
 *
 * `carl.js` is built with Emscripten's MODULARIZE option: loading it defines a global
 * `CARL` factory that, when invoked, resolves to the instantiated WASM module.  The factory
 * auto-resolves `carl.wasm` relative to the location of `carl.js` (via Emscripten's
 * `scriptDirectory`/`locateFile`), so we inject `carl.js` as a `<script>` tag rather than
 * `require()`-ing it into the consumer bundle.  Bundling `carl.js` would both bloat the
 * consumer bundle and break the automatic `scriptDirectory` detection of `carl.wasm`.
 */

/**
 * The instantiated native CARL module.  This is the untyped Emscripten module surface;
 * the strongly-typed wrappers live in `carlIntegration.ts`.
 */
// eslint-disable-next-line @typescript-eslint/no-explicit-any
export type NativeCarlModule = any;

export interface CarlInitOptions {
    /**
     * Base URL (relative or absolute) of the directory containing `carl.js` and `carl.wasm`.
     * Defaults to `"assets"`, matching the conventional `docs/assets/` deployment layout.
     */
    assetsUrl?: string;
}

type CarlFactory = () => Promise<NativeCarlModule>;

/**
 * Loads and instantiates the CARL WASM runtime.
 *
 * @param options.assetsUrl Directory containing `carl.js` + `carl.wasm` (default `"assets"`).
 * @returns A promise that resolves to the instantiated native CARL module.
 */
export function initializeCarl(options?: CarlInitOptions): Promise<NativeCarlModule> {
    const assetsUrl = options?.assetsUrl ?? "assets";
    return new Promise<NativeCarlModule>((resolve, reject) => {
        const script = document.createElement("script");
        script.type = "text/javascript";
        script.src = `${assetsUrl}/carl.js`;
        script.onload = () => {
            const factory = (window as unknown as { CARL?: CarlFactory }).CARL;
            if (!factory) {
                reject(new Error("CARL factory was not defined after loading carl.js."));
                return;
            }
            factory().then(resolve).catch(reject);
        };
        script.onerror = () => reject(new Error(`Failed to load carl.js from "${script.src}".`));
        document.body.appendChild(script);
    });
}

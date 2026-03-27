# Changelog

## [0.1.0] - 2026-03-26

### Added
- Initial release of CARL Unity integration package.
- P/Invoke bindings for the complete CARL C API.
- Managed C# wrapper classes with IDisposable lifecycle management.
- Unity components: CarlSessionManager, CarlXRInputProvider, CarlGestureRecognizer, CarlGestureRecorder.
- CarlDefinitionAsset ScriptableObject for embedding gesture definitions as Unity assets.
- Editor tooling: custom inspector for CarlDefinitionAsset, build preprocessor for native plugin validation.
- Sample scene demonstrating basic gesture recognition setup.

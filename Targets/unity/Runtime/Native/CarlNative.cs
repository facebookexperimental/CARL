/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

using System;
using System.Runtime.InteropServices;

namespace Carl.Native
{
    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    public delegate void CarlLoggerCallback(IntPtr messagePtr);

    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    public delegate void CarlRecognizerCreatedCallback(ulong requestId, ulong recognizerPtr);

    internal static class CarlNative
    {
        const string LibName = "carl";

        // --- Byte buffer ---
        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong carl_getBytes(ulong bytesPtr, IntPtr destination, ulong size);

        // --- Recording ---
        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong carl_startRecording(ulong maxSeconds);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void carl_recordObjectInputSample(ulong inProgressRecordingPtr, ref CarlInputSampleInterop sample);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void carl_recordInputSample(ulong inProgressRecordingPtr, IntPtr bytes, ulong size);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong carl_finishRecording(ulong inProgressRecordingPtr);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong carl_serializeRecording(ulong recordingPtr);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong carl_deserializeRecording(IntPtr bytes, ulong size);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void carl_disposeRecording(ulong recordingPtr);

        // --- RecordingInspector ---
        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong carl_getRecordingInspector(ulong recordingPtr);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong carl_inspect(ulong recordingInspectorPtr, double timestamp);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern double carl_getStartTimestamp(ulong recordingInspectorPtr);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern double carl_getEndTimestamp(ulong recordingInspectorPtr);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void carl_disposeRecordingInspector(ulong recordingInspectorPtr);

        // --- Example ---
        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong carl_createExample(ulong recordingPtr, double startTimestamp, double endTimestamp);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong carl_createAutoTrimmedExample(ulong recognizerPtr, ulong recordingPtr);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong carl_getRecording(ulong examplePtr);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern double carl_getExampleStartTimestamp(ulong examplePtr);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern double carl_getExampleEndTimestamp(ulong examplePtr);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void carl_disposeExample(ulong examplePtr);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong carl_loadExampleFromFile([MarshalAs(UnmanagedType.LPUTF8Str)] string path);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void carl_saveExampleToFile(ulong examplePtr, [MarshalAs(UnmanagedType.LPUTF8Str)] string path);

        // --- Definition ---
        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong carl_createDefinition(ulong descriptorType);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void carl_addExample(ulong definitionPtr, ulong recordingPtr, double startTimestamp, double endTimestamp);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void carl_addCounterexample(ulong definitionPtr, ulong recordingPtr, double startTimestamp, double endTimestamp);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern double carl_getDefaultSensitivity(ulong definitionPtr);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void carl_setDefaultSensitivity(ulong definitionPtr, double defaultSensitivity);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong carl_serializeDefinition(ulong definitionPtr);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong carl_deserializeDefinition(IntPtr bytes, ulong size);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong carl_loadDefinitionFromFile([MarshalAs(UnmanagedType.LPUTF8Str)] string path);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void carl_saveDefinitionToFile(ulong definitionPtr, [MarshalAs(UnmanagedType.LPUTF8Str)] string path);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong carl_getExamplesCount(ulong definitionPtr);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong carl_getExampleAtIdx(ulong definitionPtr, ulong idx);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong carl_getCounterexamplesCount(ulong definitionPtr);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong carl_getCounterexampleAtIdx(ulong definitionPtr, ulong idx);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void carl_disposeDefinition(ulong definitionPtr);

        // --- Session ---
        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong carl_createSession();

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong carl_createSingleThreadedSession();

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void carl_setSessionLogger(ulong sessionPtr, CarlLoggerCallback callback);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void carl_tickCallbacks(ulong sessionPtr);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void carl_addSerializedInputSample(ulong sessionPtr, IntPtr bytes, ulong size);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void carl_addInputSample(ulong sessionPtr, ref CarlInputSampleInterop sample);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void carl_disposeSession(ulong sessionPtr);

        // --- Recognizer ---
        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void carl_createRecognizerAsync(ulong sessionPtr, ulong definitionPtr, ulong requestId, CarlRecognizerCreatedCallback callback);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern double carl_getCurrentScore(ulong recognizerPtr);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void carl_setSensitivity(ulong recognizerPtr, double sensitivity);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong carl_getCanonicalRecordingInspector(ulong recognizerPtr);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void carl_disposeRecognizer(ulong sessionPtr, ulong recognizerPtr);

        // --- DefinitionBuilder ---
        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong carl_createExamplesFromRecordings(ulong actionType, IntPtr recordingPtrs, ulong count, double expectedDuration);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong carl_getBuiltExamplesCount(ulong resultPtr);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong carl_getBuiltExampleAtIdx(ulong resultPtr, ulong idx);

        [DllImport(LibName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void carl_disposeBuiltExamples(ulong resultPtr);
    }
}

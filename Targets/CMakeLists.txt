# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

add_library(carl_static_library "c_api/main.cpp")
target_link_libraries(carl_static_library PUBLIC carl_core)
if(EMSCRIPTEN)
    target_link_options(carl_static_library PUBLIC -s EXPORTED_FUNCTIONS=["_getBytes","_startRecording","_recordInputSample","_finishRecording","_serializeRecording","_deserializeRecording","_disposeRecording","_getRecordingInspector","_inspect","_getStartTimestamp","_getEndTimestamp","_disposeRecordingInspector","_createdDefinition","_addExample","_getDefaultSensitivity","_setDefaultSensitivity","_serializeDefinition","_deserializeDefinition","_getExamplesCount","_getRecordingFromExampleAtIdx","_getStartTimestampFromExampleAtIdx","_getEndTimestampFromExampleAtIdx","_disposeDefinition","_createSession","_addInputSample","_disposeSession","_createRecognizer","_getCurrentScore","_setSensitivity","_getCanonicalRecordingInspector","_disposeRecognizer","_malloc","_free"])
endif()
set_property(TARGET carl_static_library PROPERTY FOLDER "Targets")

if(NOT EMSCRIPTEN)
    add_library(carl_shared_library SHARED "c_api/main.cpp")
    target_link_libraries(carl_shared_library PUBLIC carl_core)
    set_property(TARGET carl_shared_library PROPERTY FOLDER "Targets")
    if(WIN32 OR WINDOWS_STORE)
        target_compile_definitions(carl_shared_library PRIVATE CARL_PLATFORM_WINDOWS)
    endif()
    if(ANDROID)
        set_property(TARGET carl_shared_library PROPERTY POSITION_INDEPENDENT_CODE ON)
    endif()
endif()

if(WIN32)
    add_executable(carl_test_console "console/main.cpp")
    target_link_libraries(carl_test_console PRIVATE carl_core)
    set_property(TARGET carl_test_console PROPERTY FOLDER "Targets")
    set_property(DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT carl_test_console)
endif()

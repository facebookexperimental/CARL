# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

if(ANDROID)
    find_library(android_logging log)
endif()

if(NOT EMSCRIPTEN)
    add_library(carl_static 
        "c_api/main.cpp"
        "c_api/include/carl.h")
    target_link_libraries(carl_static PUBLIC carl_core carl_utilities)
    if(WIN32 OR WINDOWS_STORE)
        target_compile_definitions(carl_static PRIVATE CARL_PLATFORM_WINDOWS)
    elseif(ANDROID)
        target_compile_definitions(carl_static PRIVATE CARL_PLATFORM_ANDROID)
        target_link_libraries(carl_static PRIVATE ${android_logging})
    endif()
    set_property(TARGET carl_static PROPERTY FOLDER "Targets")
endif()

if(EMSCRIPTEN)
    add_executable(carl
        "c_api/main.cpp"
        "c_api/include/carl.h")
    target_compile_definitions(carl PRIVATE CARL_PLATFORM_EMSCRIPTEN)
    target_link_libraries(carl PRIVATE PUBLIC carl_core carl_utilities)
    target_link_options(carl PRIVATE -sMODULARIZE=1)
    target_link_options(carl PRIVATE -sEXPORT_NAME="CARL")
    target_link_options(carl PRIVATE -sALLOW_MEMORY_GROWTH=1)
    target_link_options(carl PRIVATE --bind)
    
    # Compile optimizations
    target_compile_options(carl_input_sample PRIVATE -Oz)
    target_compile_options(carl_core PRIVATE -Oz)
    target_compile_options(carl_utilities PRIVATE -Oz)
    target_compile_options(carl PRIVATE -Oz)

    # Link optimizatins
    target_link_options(carl_input_sample PRIVATE -Oz)
    target_link_options(carl_core PRIVATE -Oz)
    target_link_options(carl_utilities PRIVATE -Oz)
    target_link_options(carl PRIVATE -Oz)
    target_link_options(carl_input_sample PRIVATE -flto)
    target_link_options(carl_core PRIVATE -flto)
    target_link_options(carl_utilities PRIVATE -flto)
    target_link_options(carl PRIVATE -flto)
else()
    add_library(carl SHARED "c_api/main.cpp")
    target_link_libraries(carl PUBLIC carl_core carl_utilities)
    set_property(TARGET carl PROPERTY FOLDER "Targets")
    if(WIN32 OR WINDOWS_STORE)
        target_compile_definitions(carl PRIVATE CARL_PLATFORM_WINDOWS)
    elseif(ANDROID)
        target_compile_definitions(carl PRIVATE CARL_PLATFORM_ANDROID)
        target_link_libraries(carl PRIVATE ${android_logging})
        set_property(TARGET carl PROPERTY POSITION_INDEPENDENT_CODE ON)
    endif()
endif()

if(WIN32)
    add_executable(carl_test_console "console/main.cpp")
    target_link_libraries(carl_test_console PRIVATE carl_core carl_utilities)
    set_property(TARGET carl_test_console PROPERTY FOLDER "Targets")
    set_property(DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT carl_test_console)
endif()

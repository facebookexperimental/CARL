# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

set(SOURCES
    "include/carl/InputSample.h"
    "source/InputSample.cpp")

add_library(carl_input_sample ${SOURCES})
target_include_directories(carl_input_sample PUBLIC "include")
target_link_libraries(carl_input_sample PUBLIC carl_core_utilities)
set_property(TARGET carl_input_sample PROPERTY FOLDER "CARL")
if (MSVC)
    target_compile_options(carl_input_sample PRIVATE /bigobj)
endif()
if(ANDROID)
    set_property(TARGET carl_input_sample PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

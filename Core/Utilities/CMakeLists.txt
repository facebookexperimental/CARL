# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

set(SOURCES
    "include/carl/utilities/FileSerialization.h"
    "source/FileSerialization.cpp")

add_library(carl_utilities ${SOURCES})
target_include_directories(carl_utilities PUBLIC "include")
target_link_libraries(carl_utilities PUBLIC carl_core)
set_property(TARGET carl_utilities PROPERTY FOLDER "CARL")
if (MSVC)
    target_compile_options(carl_utilities PRIVATE /bigobj)
endif()
if(ANDROID)
    set_property(TARGET carl_utilities PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

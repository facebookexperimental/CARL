# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

add_library(carl_core_utilities INTERFACE)
target_include_directories(carl_core_utilities INTERFACE "include")
target_link_libraries(carl_core_utilities INTERFACE arcana carl_core_foundation carl_input_sample)
set_property(TARGET carl_core_utilities PROPERTY FOLDER "CARL")
if (MSVC)
    target_compile_options(carl_core_utilities INTERFACE /bigobj)
endif()
if(ANDROID)
    set_property(TARGET carl_core_utilities PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

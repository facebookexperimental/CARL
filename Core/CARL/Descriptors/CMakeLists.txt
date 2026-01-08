# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

add_library(carl_descriptors INTERFACE)
target_include_directories(carl_descriptors INTERFACE "include")
target_link_libraries(carl_descriptors INTERFACE carl_core_utilities carl_input_sample)
set_property(TARGET carl_descriptors PROPERTY FOLDER "CARL")
if (MSVC)
    target_compile_options(carl_descriptors INTERFACE /bigobj)
endif()
if(ANDROID)
    set_property(TARGET carl_descriptors PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

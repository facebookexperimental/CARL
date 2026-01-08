# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

add_library(carl_action_type_definitions INTERFACE)
target_include_directories(carl_action_type_definitions INTERFACE "include")
target_link_libraries(carl_action_type_definitions INTERFACE carl_action_type carl_descriptors)
set_property(TARGET carl_action_type_definitions PROPERTY FOLDER "CARL")
if (MSVC)
    target_compile_options(carl_action_type_definitions INTERFACE /bigobj)
endif()
if(ANDROID)
    set_property(TARGET carl_action_type_definitions PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

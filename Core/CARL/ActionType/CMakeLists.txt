# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

add_library(carl_action_type INTERFACE)
target_include_directories(carl_action_type INTERFACE "include")
set_property(TARGET carl_action_type PROPERTY FOLDER "CARL")
if (MSVC)
    target_compile_options(carl_action_type INTERFACE /bigobj)
endif()
if(ANDROID)
    set_property(TARGET carl_action_type PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

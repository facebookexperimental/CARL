# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

cmake_minimum_required(VERSION 3.18)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

project(carl)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(Dependencies)
add_subdirectory(Core)
add_subdirectory(Targets)

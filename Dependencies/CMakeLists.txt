# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "Build pkg-config .pc file for Eigen")
set(EIGEN_SPLIT_LARGE_TESTS OFF CACHE BOOL "Split large tests into smaller executables")
set(EIGEN_BUILD_DOC OFF CACHE BOOL "Enable creation of Eigen documentation")
set(BUILD_TESTING OFF CACHE BOOL "Enable creation of Eigen tests.")
add_subdirectory(eigen)
# Prevents the inclusion of any non-permissively licensed code from Eigen.
target_compile_definitions(eigen INTERFACE EIGEN_MPL2_ONLY)
set_property(TARGET eigen PROPERTY FOLDER Dependencies)
set_property(TARGET uninstall PROPERTY FOLDER Dependencies)

# Arcana.cpp
add_subdirectory(arcana.cpp)
set_property(TARGET arcana PROPERTY FOLDER Dependencies)

# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

set(SOURCES
    "include/carl/Carl.h"
    "include/carl/Definition.h"
    "include/carl/Example.h"
    "include/carl/InputSample.h"
    "include/carl/Recognizer.h"
    "include/carl/Recording.h"
    "include/carl/Serialization.h"
    "include/carl/Session.h"
    "include/carl/Signaling.h"
    "include/carl/Types.h"
    "source/Definition.cpp"
    "source/Descriptor.h"
    "source/DynamicTimeWarping.h"
    "source/EgocentricTemporalSpace.h"
    "source/Example.cpp"
    "source/InputSample.cpp"
    "source/Recognizer.cpp"
    "source/Recording.cpp"
    "source/Session.cpp"
    "source/SessionImpl.h"
    "source/TransitoryTemporalSpace.h")

add_library(carl_core ${SOURCES})
target_include_directories(carl_core PUBLIC "include")
target_link_libraries(carl_core PUBLIC arcana eigen)
set_property(TARGET carl_core PROPERTY FOLDER "CARL")
if (MSVC)
    target_compile_options(carl_core PRIVATE /bigobj)
endif()
if(ANDROID)
    set_property(TARGET carl_core PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

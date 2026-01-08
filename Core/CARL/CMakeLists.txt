# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

set(SOURCES
    "include/carl/ActionType.h"
    "include/carl/Carl.h"
    "include/carl/ContractId.h"
    "include/carl/Definition.h"
    "include/carl/Example.h"
    "include/carl/InputSample.h"
    "include/carl/Recognizer.h"
    "include/carl/Recording.h"
    "include/carl/Serialization.h"
    "include/carl/Session.h"
    "include/carl/Signaling.h"
    "include/carl/Types.h"
    "source/descriptor/ControllerState.h"
    "source/descriptor/CombinedDescriptor.h"
    "source/descriptor/Custom.h"
    "source/descriptor/DescriptorUtils.h"
    "source/descriptor/EgocentricRelativeWristPosition.h"
    "source/descriptor/EgocentricTemporalSpace.h"
    "source/descriptor/EgocentricWristDisplacement.h"
    "source/descriptor/EgocentricWristOrientation.h"
    "source/descriptor/EgocentricWristTranslation.h"
    "source/descriptor/HandShape.h"
    "source/descriptor/SequenceOperations.h"
    "source/descriptor/TimePoint.h"
    "source/descriptor/TimestampedDescriptor.h"
    "source/descriptor/Trivial.h"
    "source/descriptor/Tuning.h"
    "source/descriptor/WristRotation.h"
    "source/ActionTypeDefinitions.h"
    "source/Definition.cpp"
    "source/Descriptor.h"
    "source/DynamicTimeWarping.h"
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

/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include "include/carl.h"

#include <carl/InputSample.h>
#include <carl/Types.h>

namespace carl::capi
{
    carl::InputSample convert(const carl_InputSample& input);
    carl_InputSample convert(const carl::InputSample& input);
}

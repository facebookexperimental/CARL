/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include <filesystem>
#include <optional>
#include <vector>

#include <gsl/span>

namespace carl::action
{
    class Example;
    class Definition;
}

namespace carl::utilities
{
    template<typename T>
    std::vector<uint8_t> Serialize(const T&);

    template<typename T>
    void SerializeToFile(const T&, std::filesystem::path);

    template<typename T>
    std::optional<T> TryDeserialize(gsl::span<const uint8_t>);

    template<typename T>
    std::optional<T> TryDeserializeFromFile(std::filesystem::path);

    template<typename T>
    std::optional<T> TryDeserializeLegacy(gsl::span<const uint8_t>);

    template<typename T>
    std::optional<T> TryDeserializeLegacyFile(std::filesystem::path);

    template<>
    std::vector<uint8_t> Serialize(const action::Example&);

    template<>
    void SerializeToFile(const action::Example&, std::filesystem::path);

    template<>
    std::optional<action::Example> TryDeserialize(gsl::span<const uint8_t>);

    template<>
    std::optional<action::Example> TryDeserializeFromFile(std::filesystem::path);

    template<>
    std::optional<action::Example> TryDeserializeLegacy(gsl::span<const uint8_t>);

    template<>
    std::optional<action::Example> TryDeserializeLegacyFile(std::filesystem::path);

    template<>
    std::vector<uint8_t> Serialize(const action::Definition&);

    template<>
    void SerializeToFile(const action::Definition&, std::filesystem::path);

    template<>
    std::optional<action::Definition> TryDeserialize(gsl::span<const uint8_t>);

    template<>
    std::optional<action::Definition> TryDeserializeFromFile(std::filesystem::path);

    template<>
    std::optional<action::Definition> TryDeserializeLegacy(gsl::span<const uint8_t>);

    template<>
    std::optional<action::Definition> TryDeserializeLegacyFile(std::filesystem::path);
}

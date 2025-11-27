#pragma once

#include <filesystem>
#include <optional>

namespace carl::action
{
    class Example;
    class Definition;
}

namespace carl::utilities
{
    template<typename T>
    void SerializeToFile(const T&, std::filesystem::path);

    template<typename T>
    std::optional<T> TryDeserializeFromFile(std::filesystem::path);

    template<typename T>
    std::optional<T> TryDeserializeLegacyFile(std::filesystem::path);

    template<>
    void SerializeToFile(const action::Example&, std::filesystem::path);

    template<>
    std::optional<action::Example> TryDeserializeFromFile(std::filesystem::path);

    template<>
    std::optional<action::Example> TryDeserializeLegacyFile(std::filesystem::path);

    template<>
    void SerializeToFile(const action::Definition&, std::filesystem::path);

    template<>
    std::optional<action::Definition> TryDeserializeFromFile(std::filesystem::path);

    template<>
    std::optional<action::Definition> TryDeserializeLegacyFile(std::filesystem::path);
}
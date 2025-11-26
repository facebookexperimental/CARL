#pragma once

#include <filesystem>
#include <optional>

namespace carl::utilities
{
    template<typename T>
    void SerializeToFile(const T&, std::filesystem::path);

    template<typename T>
    std::optional<T> TryDeserializeFromFile(std::filesystem::path);

    template<typename T>
    std::optional<T> TryDeserializeLegacyFile(std::filesystem::path);
}
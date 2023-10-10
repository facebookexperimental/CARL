/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <carl/carl.h>

#include <filesystem>
#include <fstream>
#include <iostream>

namespace
{
    carl::action::Definition loadDefinition(const char* rawPath)
    {
        std::filesystem::path path{ rawPath };
        auto length = std::filesystem::file_size(path);

        std::ifstream fileStream{ path.c_str(), std::ios::in | std::ios::binary };
        std::vector<uint8_t> bytes{};
        bytes.resize(static_cast<size_t>(length));
        fileStream.seekg(std::ios::beg);
        fileStream.read(reinterpret_cast<char*>(bytes.data()), bytes.size());

        carl::Deserialization deserialization{ bytes.data() };
        deserialization >> std::array<uint8_t, 4>{};
        return { carl::Deserialization{deserialization} };
    }

    carl::action::Example loadExample(const char* rawPath)
    {
        std::filesystem::path path{ rawPath };
        auto length = std::filesystem::file_size(path);

        std::ifstream fileStream{ path.c_str(), std::ios::in | std::ios::binary };
        std::vector<uint8_t> bytes{};
        bytes.resize(static_cast<size_t>(length));
        fileStream.seekg(std::ios::beg);
        fileStream.read(reinterpret_cast<char*>(bytes.data()), bytes.size());

        carl::Deserialization deserialization{ bytes.data() };
        float startTimestamp;
        deserialization >> startTimestamp;
        float endTimestamp;
        deserialization >> endTimestamp;
        deserialization >> std::array<uint8_t, 4>{};
        carl::action::Recording recording{ deserialization };
        return { recording, static_cast<double>(startTimestamp), static_cast<double>(endTimestamp) };
    }
}

void test()
{
    auto example0{ loadExample("C:\\scratch\\CARLFiles\\recording_0.bin") };
    auto example1{ loadExample("C:\\scratch\\CARLFiles\\recording_1.bin") };

    carl::action::Definition definition{
        carl::action::Definition::ActionType::RightHandGesture };
    definition.addExample(example0);

    carl::Session session{};
    carl::action::Recognizer recognizer{ session, definition };

    auto autoTrimmedExample = recognizer.createAutoTrimmedExample(example1.getRecording());

    for (const auto& sample : example1.getRecording().getSamples())
    {
        session.addInput(sample);
        std::cout << recognizer.currentScore() << std::endl;
    }
}

void main()
{
    test();
}

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
    std::vector<carl::action::Example> pullRecordings{};
    for (size_t idx = 0; idx < 10; ++idx)
    {
        std::stringstream ss{};
        ss << "C:\\scratch\\CARLFiles\\pull_recordings\\recording_" << idx << ".bin";
        pullRecordings.emplace_back(loadExample(ss.str().c_str()));
    }

    carl::action::Definition pullDefinition{
        carl::action::Definition::ActionType::RightHandGesture };
    for (const auto& example : pullRecordings)
    {
        pullDefinition.addExample(example);
    }

    std::vector<carl::action::Example> nodRecordings{};
    for (size_t idx = 0; idx < 10; ++idx)
    {
        std::stringstream ss{};
        ss << "C:\\scratch\\CARLFiles\\nod_recordings\\recording_" << idx << ".bin";
        nodRecordings.emplace_back(loadExample(ss.str().c_str()));
    }

    carl::action::Definition nodDefinition{
        carl::action::Definition::ActionType::RightHandGesture };
    for (const auto& example : nodRecordings)
    {
        nodDefinition.addExample(example);
    }

    carl::action::Definition sumDefinition{
        carl::action::Definition::ActionType::RightHandGesture };
    for (const auto& example : pullRecordings)
    {
        sumDefinition.addExample(example);
    }
    for (const auto& example : nodRecordings)
    {
        sumDefinition.addExample(example);
    }

    carl::Session session{};
    //carl::action::Recognizer pullRecognizer{ session, pullDefinition };
    //carl::action::Recognizer nodRecognizer{ session, nodDefinition };
    carl::action::Recognizer summRecognizer{ session, sumDefinition };

    //auto autoTrimmedExample = recognizer.createAutoTrimmedExample(example1.getRecording());

    /*for (const auto& sample : pullRecordings[0].getRecording().getSamples())
    {
        session.addInput(sample);
        std::cout << pullRecognizer.currentScore() << std::endl;
    }*/
}

void test2()
{
    auto definition = loadDefinition("C:\\scratch\\CARLFiles\\definition_2.bin");
    auto recording = definition.getExamples().front().getRecording();

    int idx = 0;
    for (; recording.getSamples()[idx + 1].Timestamp < recording.getInspector().endTimestamp(); ++idx);
    auto sample = recording.getSamples()[idx];

    carl::Session session{};
    carl::action::Recognizer recognizer{ session, definition };
    for (idx = 0; idx < 1000; ++idx)
    {
        if (idx == 900)
        {
            std::cout << "Ready to test" << std::endl;
        }

        auto newSample{ sample };
        newSample.Timestamp = idx * 0.05;
        session.addInput(newSample);
    }
}

void test3()
{
    auto example0 = loadExample("C:\\scratch\\CARLFiles\\push_recordings\\recording_0.bin");
    auto example1 = loadExample("C:\\scratch\\CARLFiles\\push_recordings\\recording_1.bin");

    carl::action::Definition definition{ carl::action::Definition::ActionType::RightHandGesture };
    definition.addExample(example0);

    carl::Session session{};
    
    carl::action::Recognizer recognizer{ session, definition };

    recognizer.analyzeRecording(example0.getRecording(), std::cout);
}

void main()
{
    //test();
    //test2();
    test3();
}

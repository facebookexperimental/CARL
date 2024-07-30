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
    carl::action::Definition loadDefinition(std::filesystem::path path)
    {
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

    carl::action::Example loadExample(std::filesystem::path path)
    {
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

    void ensureDirectory(const std::filesystem::path& path)
    {
        if (!std::filesystem::is_directory(path))
        {
            auto parent = path.parent_path();
            ensureDirectory(parent);
            std::filesystem::create_directory(path);
        }
    }

    void runAnalysis(int argc, const char** argv, std::filesystem::path executionDirectory)
    {
        std::filesystem::path outputDirectory{ argv[0] };
        if (!outputDirectory.is_absolute())
        {
            outputDirectory = executionDirectory / outputDirectory;
        }

        carl::action::Definition::ActionType definitionType{};
        std::string definitionTypeArg{ argv[1] };
        if (definitionTypeArg == "lhp")
        {
            definitionType = carl::action::Definition::ActionType::LeftHandPose;
        }
        else if (definitionTypeArg == "rhp")
        {
            definitionType = carl::action::Definition::ActionType::RightHandPose;
        }
        else if (definitionTypeArg == "lhg")
        {
            definitionType = carl::action::Definition::ActionType::LeftHandGesture;
        }
        else if (definitionTypeArg == "rhg")
        {
            definitionType = carl::action::Definition::ActionType::RightHandGesture;
        }
        else if (definitionTypeArg == "thg")
        {
            definitionType = carl::action::Definition::ActionType::TwoHandGesture;
        }
        else if (definitionTypeArg == "lcg")
        {
            definitionType = carl::action::Definition::ActionType::LeftControllerGesture;
        }
        else if (definitionTypeArg == "rcg")
        {
            definitionType = carl::action::Definition::ActionType::RightControllerGesture;
        }
        else if (definitionTypeArg == "tcg")
        {
            definitionType = carl::action::Definition::ActionType::TwoControllerGesture;
        }

        enum class State
        {
            Examples,
            Sessions,
            Ambient,
        };
        auto currentState = State::Examples;

        std::vector<carl::action::Example> exampleExamples{};
        std::vector<carl::action::Example> sessionExamples{};
        std::vector<carl::action::Example> ambientExamples{};
        for (size_t idx = 2; idx < static_cast<size_t>(argc); ++idx)
        {
            std::string arg{ argv[idx] };
            if (arg == "examples")
            {
                currentState = State::Examples;
            }
            else if (arg == "sessions")
            {
                currentState = State::Sessions;
            }
            else if (arg == "ambient")
            {
                currentState = State::Ambient;
            }
            else
            {
                std::filesystem::path path{ arg };
                if (path.is_relative())
                {
                    path = executionDirectory / path;
                }

                auto storeExample = [&](std::filesystem::path examplePath) {
                    auto example = loadExample(examplePath);
                    switch (currentState)
                    {
                    case State::Examples:
                        exampleExamples.push_back(std::move(example));
                        break;
                    case State::Sessions:
                        sessionExamples.push_back(std::move(example));
                        break;
                    case State::Ambient:
                        ambientExamples.push_back(std::move(example));
                        break;
                    }
                };

                if (arg[arg.size() - 1] == '*')
                {
                    // Expand the wildcard
                    auto prefix = path.filename().string();
                    prefix = prefix.substr(0, prefix.size() - 1);
                    for (const auto& entry : std::filesystem::directory_iterator{ path.parent_path() })
                    {
                        if (entry.path().has_filename() && entry.path().filename().string().find(prefix) == 0)
                        {
                            storeExample(entry.path());
                        }
                    }
                }
                else
                {
                    storeExample(path);
                }
            }
        }

        constexpr auto analyze = [](std::filesystem::path outputDirectory, carl::action::Recognizer& recognizer, std::vector<carl::action::Example>& examples) {
            ensureDirectory(outputDirectory);
            size_t preExistingFilesCount = 0;
            for (const auto& _ : std::filesystem::directory_iterator(outputDirectory))
            {
                ++preExistingFilesCount;
            }

            for (size_t idx = 0; idx < examples.size(); ++idx)
            {
                const auto& recording = examples[idx].getRecording();
                auto outputFilePath = outputDirectory / (std::to_string(preExistingFilesCount + idx) + ".csv");
                std::ofstream outStream{ outputFilePath };
                recognizer.analyzeRecording(recording, outStream);
            }
        };

        carl::Session session{ true };

        // Test each example against all the others
        for (size_t i = 0; i < exampleExamples.size() - 2; ++i)
        {
            carl::action::Definition def{ definitionType };
            def.addExample(exampleExamples[i]);
            carl::action::Recognizer rec{ session, def };
            analyze(outputDirectory / "examples", rec, exampleExamples);
        }

        // Build a definition from the examples
        carl::action::Definition definition{ definitionType };
        for (const auto& example : exampleExamples)
        {
            definition.addExample(example);
        }
        carl::action::Recognizer recognizer{ session, definition };

        analyze(outputDirectory / "sessions", recognizer, sessionExamples);
        analyze(outputDirectory / "ambient", recognizer, ambientExamples);
    }

    void runAnalysisFile(std::string analysisFilePath)
    {
        std::filesystem::path path{ analysisFilePath };
        std::ifstream fileStream{ path.c_str(), std::ios::in | std::ios::binary };

        for (std::string line{}; std::getline(fileStream, line);)
        {
            static const std::string disallowedChars{ "\r\n" };
            static const auto isCharDisallowed = [](const char c) {
                return disallowedChars.find(c) < disallowedChars.size();
            };

            std::vector<std::string> tokens{};
            bool shouldAddToken = true;
            for (const char c : line)
            {
                if (std::isspace(c))
                {
                    shouldAddToken = true;
                }
                else if (!isCharDisallowed(c))
                {
                    if (shouldAddToken)
                    {
                        tokens.emplace_back();
                        shouldAddToken = false;
                    }

                    tokens.back().push_back(c);
                }
            }

            std::vector<const char*> args{};
            for (const auto& token : tokens)
            {
                args.push_back(token.data());
            }
            runAnalysis(static_cast<int>(args.size()), args.data(), path.parent_path());
        }
    }
}

/*
Tuning profile: Let the user provide a set of definitions (all of which contain examples) as
well as independent examples. The definitions (obviously) represent their own action, and the
independent examples represent no action. The system then analyzes every example against
every other example, determining the expected distances across descriptor dimensions, trying
to arrive at an optimal tuning.

Instead of DTW, do a naive sequence match? Fail out if any single connection is too large?
*/

void main(int argc, const char** argv)
{
    try
    {
        switch (argc)
        {
        case 2:
        {
            runAnalysisFile(argv[1]);
            break;
        }
        default:
        {
            std::filesystem::path exePath{ argv[0] };
            runAnalysis(argc - 1, argv + 1, exePath.parent_path());
            break;
        }
        }
    }
    catch (...)
    {
        std::cout << "TODO: Print usage message" << std::endl;
    }
}

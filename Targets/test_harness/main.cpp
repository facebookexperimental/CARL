/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <carl/Definition.h>
#include <carl/DefinitionBuilder.h>
#include <carl/Example.h>
#include <carl/Recognizer.h>
#include <carl/Recording.h>
#include <carl/Session.h>
#include <carl/utilities/FileSerialization.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

struct ExampleInfo
{
    carl::action::Example example;
    std::string filePath;
    std::string folderName;
};

struct DefinitionInfo
{
    carl::action::Definition definition;
    std::string filePath;
    std::string folderName;
};

struct PerformanceMetrics
{
    size_t frameCount{};
    double averageFrameTimeUs{};
    double maxFrameTimeUs{};
    double p95FrameTimeUs{};
    double p99FrameTimeUs{};
    double minFrameTimeUs{};
    double totalTimeMs{};
};

struct AssociatedRecognitionMetrics
{
    double truePositiveWindowStart{};
    double truePositiveWindowEnd{};
    double maxConfidenceInWindow{};
    double timestampOfMaxInWindow{};
    double maxConfidenceOutsideWindow{};
};

struct NonAssociatedRecognitionMetrics
{
    double maxConfidence{};
    double timestampOfMaxConfidence{};
};

struct TestResult
{
    std::string definitionPath;
    std::string examplePath;
    bool associated{};
    PerformanceMetrics performance;
    std::optional<AssociatedRecognitionMetrics> associatedRecognition;
    std::optional<NonAssociatedRecognitionMetrics> nonAssociatedRecognition;
};

double computePercentile(std::vector<double>& sortedValues, double percentile)
{
    if (sortedValues.empty()) return 0.0;
    double idx = (percentile / 100.0) * static_cast<double>(sortedValues.size() - 1);
    size_t lower = static_cast<size_t>(std::floor(idx));
    size_t upper = static_cast<size_t>(std::ceil(idx));
    if (lower == upper) return sortedValues[lower];
    double frac = idx - static_cast<double>(lower);
    return sortedValues[lower] * (1.0 - frac) + sortedValues[upper] * frac;
}

PerformanceMetrics computePerformanceMetrics(std::vector<double>& frameTimes)
{
    PerformanceMetrics metrics;
    metrics.frameCount = frameTimes.size();
    if (frameTimes.empty()) return metrics;

    double total = std::accumulate(frameTimes.begin(), frameTimes.end(), 0.0);
    metrics.totalTimeMs = total / 1000.0;
    metrics.averageFrameTimeUs = total / static_cast<double>(frameTimes.size());

    std::sort(frameTimes.begin(), frameTimes.end());
    metrics.minFrameTimeUs = frameTimes.front();
    metrics.maxFrameTimeUs = frameTimes.back();
    metrics.p95FrameTimeUs = computePercentile(frameTimes, 95.0);
    metrics.p99FrameTimeUs = computePercentile(frameTimes, 99.0);

    return metrics;
}

std::string makeRelativePath(const fs::path& filePath, const fs::path& baseDir)
{
    auto rel = filePath.lexically_relative(baseDir);
    return rel.generic_string();
}

std::string escapeJsonString(const std::string& s)
{
    std::string result;
    result.reserve(s.size());
    for (char c : s)
    {
        switch (c)
        {
        case '"': result += "\\\""; break;
        case '\\': result += "\\\\"; break;
        case '\n': result += "\\n"; break;
        case '\r': result += "\\r"; break;
        case '\t': result += "\\t"; break;
        default: result += c;
        }
    }
    return result;
}

std::string getCurrentTimestamp()
{
    auto now = std::chrono::system_clock::now();
    auto time = std::chrono::system_clock::to_time_t(now);
    std::tm tm{};
#ifdef _WIN32
    gmtime_s(&tm, &time);
#else
    gmtime_r(&time, &tm);
#endif
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%dT%H:%M:%SZ");
    return oss.str();
}

void writeJsonNumber(std::ostream& out, const char* key, double value, bool comma = true)
{
    out << "        \"" << key << "\": ";
    if (std::isnan(value) || std::isinf(value))
    {
        out << "0";
    }
    else
    {
        out << std::fixed << std::setprecision(4) << value;
    }
    if (comma) out << ",";
    out << "\n";
}

TestResult runPair(
    const DefinitionInfo& defInfo,
    const ExampleInfo& exInfo,
    double tolerance,
    const fs::path& testDataDir)
{
    TestResult result;
    result.definitionPath = makeRelativePath(defInfo.filePath, testDataDir);
    result.examplePath = makeRelativePath(exInfo.filePath, testDataDir);
    result.associated = (defInfo.folderName == exInfo.folderName);

    carl::Session session(true);
    carl::action::Recognizer recognizer(session, defInfo.definition);

    const auto& samples = exInfo.example.getRecording().getSamples();
    std::vector<double> frameTimes;
    frameTimes.reserve(samples.size());

    struct ScoreAtTime
    {
        double timestamp;
        double score;
    };
    std::vector<ScoreAtTime> scores;
    scores.reserve(samples.size());

    arcana::cancellation_source cancellationSource{};

    for (const auto& sample : samples)
    {
        auto start = std::chrono::steady_clock::now();
        session.addInput(sample);
        session.tickCallbacks(cancellationSource);
        auto end = std::chrono::steady_clock::now();

        double elapsedUs = std::chrono::duration<double, std::micro>(end - start).count();
        frameTimes.push_back(elapsedUs);
        scores.push_back({ sample.Timestamp, recognizer.currentScore() });
    }

    result.performance = computePerformanceMetrics(frameTimes);

    if (result.associated)
    {
        AssociatedRecognitionMetrics rec;
        rec.truePositiveWindowStart = exInfo.example.getEndTimestamp();
        rec.truePositiveWindowEnd = exInfo.example.getEndTimestamp() + tolerance;
        rec.maxConfidenceInWindow = 0.0;
        rec.timestampOfMaxInWindow = rec.truePositiveWindowStart;
        rec.maxConfidenceOutsideWindow = 0.0;

        for (const auto& s : scores)
        {
            if (s.timestamp >= rec.truePositiveWindowStart && s.timestamp <= rec.truePositiveWindowEnd)
            {
                if (s.score > rec.maxConfidenceInWindow)
                {
                    rec.maxConfidenceInWindow = s.score;
                    rec.timestampOfMaxInWindow = s.timestamp;
                }
            }
            else
            {
                rec.maxConfidenceOutsideWindow = std::max(rec.maxConfidenceOutsideWindow, s.score);
            }
        }

        result.associatedRecognition = rec;
    }
    else
    {
        NonAssociatedRecognitionMetrics rec;
        rec.maxConfidence = 0.0;
        rec.timestampOfMaxConfidence = 0.0;

        for (const auto& s : scores)
        {
            if (s.score > rec.maxConfidence)
            {
                rec.maxConfidence = s.score;
                rec.timestampOfMaxConfidence = s.timestamp;
            }
        }

        result.nonAssociatedRecognition = rec;
    }

    return result;
}

bool writeResultsJson(
    const std::vector<TestResult>& results,
    const std::string& outputPath,
    const std::string& commitHash,
    double tolerance)
{
    std::ofstream out(outputPath);
    if (!out.is_open())
    {
        std::cerr << "Error: Could not open output file: " << outputPath << "\n";
        return false;
    }

    out << "{\n";
    out << "  \"version\": \"1.0\",\n";
    out << "  \"timestamp\": \"" << getCurrentTimestamp() << "\",\n";
    out << "  \"commitHash\": \"" << escapeJsonString(commitHash) << "\",\n";
    out << "  \"truePositiveWindowToleranceSeconds\": " << std::fixed << std::setprecision(1) << tolerance << ",\n";
    out << "  \"results\": [\n";

    for (size_t i = 0; i < results.size(); ++i)
    {
        const auto& r = results[i];
        out << "    {\n";
        out << "      \"definition\": \"" << escapeJsonString(r.definitionPath) << "\",\n";
        out << "      \"example\": \"" << escapeJsonString(r.examplePath) << "\",\n";
        out << "      \"associated\": " << (r.associated ? "true" : "false") << ",\n";

        out << "      \"performance\": {\n";
        out << "        \"frameCount\": " << r.performance.frameCount << ",\n";
        writeJsonNumber(out, "averageFrameTimeUs", r.performance.averageFrameTimeUs);
        writeJsonNumber(out, "maxFrameTimeUs", r.performance.maxFrameTimeUs);
        writeJsonNumber(out, "p95FrameTimeUs", r.performance.p95FrameTimeUs);
        writeJsonNumber(out, "p99FrameTimeUs", r.performance.p99FrameTimeUs);
        writeJsonNumber(out, "minFrameTimeUs", r.performance.minFrameTimeUs);
        writeJsonNumber(out, "totalTimeMs", r.performance.totalTimeMs, false);
        out << "      },\n";

        out << "      \"recognition\": {\n";
        if (r.associated && r.associatedRecognition)
        {
            const auto& rec = *r.associatedRecognition;
            writeJsonNumber(out, "truePositiveWindowStart", rec.truePositiveWindowStart);
            writeJsonNumber(out, "truePositiveWindowEnd", rec.truePositiveWindowEnd);
            writeJsonNumber(out, "maxConfidenceInWindow", rec.maxConfidenceInWindow);
            writeJsonNumber(out, "timestampOfMaxInWindow", rec.timestampOfMaxInWindow);
            writeJsonNumber(out, "maxConfidenceOutsideWindow", rec.maxConfidenceOutsideWindow, false);
        }
        else if (r.nonAssociatedRecognition)
        {
            const auto& rec = *r.nonAssociatedRecognition;
            writeJsonNumber(out, "maxConfidence", rec.maxConfidence);
            writeJsonNumber(out, "timestampOfMaxConfidence", rec.timestampOfMaxConfidence, false);
        }
        out << "      }\n";

        out << "    }";
        if (i + 1 < results.size()) out << ",";
        out << "\n";
    }

    out << "  ]\n";
    out << "}\n";

    std::cout << "Results written to: " << outputPath << "\n";
    return true;
}

struct ParsedArgs
{
    std::string testDataDir;
    std::string outputPath = "test_results.json";
    double tolerance = 0.2;
    std::string commitHash = "unknown";
    uint64_t buildActionType = 3; // RightHandGesture
};

std::optional<ParsedArgs> parseArgs(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr << "Usage: carl_test_harness <test_data_dir> [options]\n"
                  << "Options:\n"
                  << "  --output <path>      Output JSON path (default: test_results.json)\n"
                  << "  --tolerance <sec>    True positive window tolerance in seconds (default: 1.0)\n"
                  << "  --commit <hash>      Commit hash to include in output (default: unknown)\n"
                  << "  --build-action-type <id>  ActionType for build_ tests (default: 3 = RightHandGesture)\n";
        return std::nullopt;
    }

    ParsedArgs args;
    args.testDataDir = argv[1];

    for (int i = 2; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg == "--output" && i + 1 < argc)
        {
            args.outputPath = argv[++i];
        }
        else if (arg == "--tolerance" && i + 1 < argc)
        {
            try
            {
                args.tolerance = std::stod(argv[++i]);
            }
            catch (const std::exception&)
            {
                std::cerr << "Invalid tolerance value: " << argv[i] << "\n";
                return std::nullopt;
            }
        }
        else if (arg == "--commit" && i + 1 < argc)
        {
            args.commitHash = argv[++i];
        }
        else if (arg == "--build-action-type" && i + 1 < argc)
        {
            try
            {
                args.buildActionType = std::stoull(argv[++i]);
            }
            catch (const std::exception&)
            {
                std::cerr << "Invalid build-action-type value: " << argv[i] << "\n";
                return std::nullopt;
            }
        }
        else
        {
            std::cerr << "Unknown argument: " << arg << "\n";
            return std::nullopt;
        }
    }

    return args;
}

int main(int argc, char* argv[])
{
    auto parsedArgs = parseArgs(argc, argv);
    if (!parsedArgs) return 1;
    const auto& args = *parsedArgs;

    fs::path testDataDir(args.testDataDir);
    if (!fs::exists(testDataDir) || !fs::is_directory(testDataDir))
    {
        std::cerr << "Error: Test data directory does not exist: " << testDataDir << "\n";
        return 1;
    }

    std::vector<DefinitionInfo> allDefinitions;
    std::vector<ExampleInfo> allExamples;

    std::cout << "Scanning test data directory: " << testDataDir << "\n";

    for (const auto& subdir : fs::directory_iterator(testDataDir))
    {
        if (!subdir.is_directory()) continue;

        std::string folderName = subdir.path().filename().string();
        std::cout << "  Scanning folder: " << folderName << "\n";

        for (const auto& entry : fs::directory_iterator(subdir.path()))
        {
            if (!entry.is_regular_file()) continue;
            if (entry.path().extension() != ".carl") continue;

            auto defOpt = carl::utilities::TryDeserializeFromFile<carl::action::Definition>(entry.path());
            if (defOpt.has_value())
            {
                std::cout << "    Loaded Definition: " << entry.path().filename() << "\n";
                allDefinitions.push_back({ std::move(*defOpt), entry.path().string(), folderName });
                continue;
            }

            auto exOpt = carl::utilities::TryDeserializeFromFile<carl::action::Example>(entry.path());
            if (exOpt.has_value())
            {
                std::cout << "    Loaded Example: " << entry.path().filename()
                          << " (samples: " << exOpt->getRecording().getSamples().size()
                          << ", action: " << exOpt->getStartTimestamp() << "s - " << exOpt->getEndTimestamp() << "s)\n";
                allExamples.push_back({ std::move(*exOpt), entry.path().string(), folderName });
                continue;
            }

            std::cerr << "    Warning: Could not load " << entry.path().filename() << " as Definition or Example\n";
        }
    }

    std::cout << "\nLoaded " << allDefinitions.size() << " definitions and "
              << allExamples.size() << " examples.\n";

    // Process build_ directories: auto-build definitions using DefinitionBuilder
    constexpr auto BUILD_PREFIX = "build_";
    constexpr size_t BUILD_PREFIX_LEN = 6;
    constexpr size_t BUILD_GROUP_SIZE = 3;

    auto buildActionType = static_cast<carl::action::ActionType>(args.buildActionType);

    for (const auto& subdir : fs::directory_iterator(testDataDir))
    {
        if (!subdir.is_directory()) continue;

        std::string folderName = subdir.path().filename().string();
        if (folderName.rfind(BUILD_PREFIX, 0) != 0) continue;

        std::string identity = folderName.substr(BUILD_PREFIX_LEN);
        std::cout << "\n  DefinitionBuilder test: " << folderName << " (identity: " << identity << ")\n";

        // Collect filenames from this build_ directory, sorted
        std::vector<fs::path> buildFiles;

        for (const auto& entry : fs::directory_iterator(subdir.path()))
        {
            if (!entry.is_regular_file()) continue;
            if (entry.path().extension() != ".carl") continue;
            buildFiles.push_back(entry.path());
        }

        std::sort(buildFiles.begin(), buildFiles.end());

        // Load examples in sorted order
        std::vector<carl::action::Example> buildExamples;
        for (const auto& filePath : buildFiles)
        {
            auto exOpt = carl::utilities::TryDeserializeFromFile<carl::action::Example>(filePath);
            if (exOpt.has_value())
            {
                std::cout << "    Loaded: " << filePath.filename() << "\n";
                buildExamples.push_back(std::move(*exOpt));

                // Also add to allExamples for the test matrix, using the post-prefix identity
                auto exOpt2 = carl::utilities::TryDeserializeFromFile<carl::action::Example>(filePath);
                if (exOpt2.has_value())
                {
                    allExamples.push_back({ std::move(*exOpt2), filePath.string(), identity });
                }
            }
        }

        if (buildExamples.size() < BUILD_GROUP_SIZE)
        {
            std::cerr << "    Warning: " << folderName << " has " << buildExamples.size()
                      << " examples, need at least " << BUILD_GROUP_SIZE << ". Skipping.\n";
            continue;
        }

        // Extract raw recordings (full recording, ignoring example boundaries)
        std::vector<carl::action::Recording> allRecordings;
        allRecordings.reserve(buildExamples.size());
        for (const auto& example : buildExamples)
        {
            allRecordings.push_back(example.getRecording());
        }

        // Group into non-overlapping triples and build definitions
        size_t groupIdx = 0;
        for (size_t i = 0; i + BUILD_GROUP_SIZE <= allRecordings.size(); i += BUILD_GROUP_SIZE)
        {
            gsl::span<const carl::action::Recording> group{ &allRecordings[i], BUILD_GROUP_SIZE };

            std::cout << "    Building definition from recordings " << i << "-" << (i + BUILD_GROUP_SIZE - 1) << "... ";
            std::cout.flush();

            auto start = std::chrono::steady_clock::now();
            auto builtExamples = carl::action::DefinitionBuilder::createExamplesFromRecordings(
                buildActionType, group, 0.);
            auto end = std::chrono::steady_clock::now();
            double elapsedMs = std::chrono::duration<double, std::milli>(end - start).count();

            std::cout << "done (" << builtExamples.size() << " examples, "
                      << std::fixed << std::setprecision(1) << elapsedMs << "ms)\n";

            if (builtExamples.empty())
            {
                std::cerr << "    Warning: DefinitionBuilder returned no examples for group " << groupIdx << "\n";
                ++groupIdx;
                continue;
            }

            carl::action::Definition def{ buildActionType };
            for (auto& ex : builtExamples)
            {
                def.addExample(std::move(ex));
            }

            std::string syntheticPath = (subdir.path() / ("auto_built_" + std::to_string(groupIdx))).string();
            allDefinitions.push_back({ std::move(def), syntheticPath, identity });
            ++groupIdx;
        }
    }

    std::cout << "\nTotal after builder: " << allDefinitions.size() << " definitions and "
              << allExamples.size() << " examples.\n";

    if (allDefinitions.empty() || allExamples.empty())
    {
        std::cout << "No pairs to test. Writing empty results.\n";
        return writeResultsJson({}, args.outputPath, args.commitHash, args.tolerance) ? 0 : 1;
    }

    size_t totalPairs = allDefinitions.size() * allExamples.size();
    std::cout << "Running " << totalPairs << " test pairs ("
              << allDefinitions.size() << " definitions x "
              << allExamples.size() << " examples)...\n\n";

    std::vector<TestResult> results;
    results.reserve(totalPairs);

    size_t pairIdx = 0;
    for (const auto& defInfo : allDefinitions)
    {
        for (const auto& exInfo : allExamples)
        {
            ++pairIdx;
            bool associated = (defInfo.folderName == exInfo.folderName);
            std::cout << "[" << pairIdx << "/" << totalPairs << "] "
                      << defInfo.folderName << " vs " << exInfo.folderName
                      << " (" << (associated ? "associated" : "non-associated") << ")... ";
            std::cout.flush();

            auto result = runPair(defInfo, exInfo, args.tolerance, testDataDir);

            std::cout << "done. Frames: " << result.performance.frameCount
                      << ", Avg: " << std::fixed << std::setprecision(1)
                      << result.performance.averageFrameTimeUs << "us";
            if (result.associated && result.associatedRecognition)
            {
                std::cout << ", MaxInWindow: " << std::setprecision(3)
                          << result.associatedRecognition->maxConfidenceInWindow;
            }
            else if (result.nonAssociatedRecognition)
            {
                std::cout << ", MaxConfidence: " << std::setprecision(3)
                          << result.nonAssociatedRecognition->maxConfidence;
            }
            std::cout << "\n";

            results.push_back(std::move(result));
        }
    }

    std::cout << "\nAll pairs completed.\n";
    if (!writeResultsJson(results, args.outputPath, args.commitHash, args.tolerance))
    {
        return 1;
    }

    return 0;
}

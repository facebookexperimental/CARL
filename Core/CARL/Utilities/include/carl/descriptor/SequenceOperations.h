/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include <gsl/span>

namespace carl::descriptor
{
    template<typename DescriptorT, typename CallableT>
    void extendSequence(const InputSample& newSample, std::vector<DescriptorT>& sequence, InputSample& mostRecentSample, gsl::span<const NumberT> tuning, const CallableT& tryCreate)
    {
        constexpr NumberT THRESHOLD{ 0.5 };

        // Handle startup, in which case m_mostRecentSample will be a default value.
        if (sequence.empty())
        {
            auto descriptor = tryCreate(newSample, newSample);
            if (descriptor.has_value())
            {
                sequence.emplace_back(std::move(*descriptor));
                mostRecentSample = newSample;
            }
        }
        else
        {
            // Iteratively create new mostRecentSamples and descriptors until the most recent descriptor is sufficiently close to that of newSample.
            // TODO: Figure out if delta descriptors (rotation, translation, etc.) will play correctly with this, since lerping inputs should not 
            //       change them. Might be necessary to mute such descriptors with tuning.
            constexpr size_t MAX_FILL_ITERATIONS = 12;
            for (size_t fillIter = 0; fillIter < MAX_FILL_ITERATIONS; ++fillIter)
            {
                auto sampleDesc = tryCreate(newSample, mostRecentSample);

                // If we weren't able to create a descriptor, stop iterating.
                if (!sampleDesc.has_value())
                {
                    break;
                }

                auto outerDistance = DescriptorT::Distance(*sampleDesc, sequence.back(), sequence.back(), sequence.back(), tuning);
                if (outerDistance < THRESHOLD)
                {
                    break;
                }

                // Posit new samples, and associated descriptors, to bring the end of sequence closer to newSample.
                NumberT upper = 1;
                NumberT lower = 0;
                NumberT mid{};
                constexpr size_t MAX_BISECT_ITERATIONS = 20;
                for (size_t bisectIter = 0; bisectIter < MAX_BISECT_ITERATIONS; ++bisectIter)
                {
                    mid = (upper + lower) / NumberT{ 2 };
                    auto intermediateDesc = DescriptorT::Lerp(sequence.back(), *sampleDesc, mid);
                    auto distance = DescriptorT::Distance(intermediateDesc, sequence.back(), sequence.back(), sequence.back(), tuning);
                    if (distance > THRESHOLD)
                    {
                        // intermediate sample is too distant, continue searching for a nearer sample
                        upper = mid;
                    }
                    else if (distance < NumberT{ 0.9 } *THRESHOLD)
                    {
                        // intermediate sample is too close, continue searching for a more distant sample
                        lower = mid;
                    }
                    else
                    {
                        // intermediate sample is satisfactory, stop searching
                        sequence.push_back(intermediateDesc);
                        mostRecentSample = InputSample::Lerp(mostRecentSample, newSample, mid);
                        break;
                    }
                }
            }
        }
    }

    template<typename DescriptorT, typename RecordingT, typename CallableT>
    std::vector<DescriptorT> createDescriptorSequenceFromRecording(const RecordingT& recording, double startTimestamp, double endTimestamp, gsl::span<const NumberT> tuning, const CallableT& tryCreate)
    {
        std::vector<DescriptorT> sequence{};

        auto samples = recording.getSamples();
        InputSample mostRecentSample{};

        size_t idx = 0;
        while (idx < samples.size() && samples[idx].Timestamp < startTimestamp)
        {
            ++idx;
        }
        do
        {
            descriptor::extendSequence(samples[idx], sequence, mostRecentSample, tuning, tryCreate);
            ++idx;
        } while (idx < samples.size() && samples[idx].Timestamp < endTimestamp);

        return sequence;
    }

    template<typename DescriptorT, typename ExampleT, typename CallableT>
    std::vector<std::vector<DescriptorT>> createDescriptorSequencesFromExamples(gsl::span<const ExampleT> examples, gsl::span<const NumberT> tuning, const CallableT& tryCreate, double padding = 0.)
    {
        std::vector<std::vector<DescriptorT>> sequences{};
        sequences.reserve(examples.size());
        for (const auto& example : examples)
        {
            sequences.push_back(createDescriptorSequenceFromRecording<DescriptorT>(example.getRecording(), example.getStartTimestamp() - padding, example.getEndTimestamp() + padding, tuning, tryCreate));
        }
        return sequences;
    }
}

/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include <carl/Definition.h>
#include <carl/Session.h>
#include <carl/Signaling.h>

#include <ostream>

namespace carl::action
{
    /// <summary>
    /// Observes a Session and attempts to recognize an action specified by a Definition.
    /// </summary>
    class Recognizer
    {
    public:
        class Impl;

        Recognizer(Session&, const Definition&);
        ~Recognizer();

        double currentScore();

        RecordingInspector getCanonicalRecordingInspector() const;

        void setSensitivity(double sensitivity);

        Example createAutoTrimmedExample(const Recording&) const;

        void analyzeRecording(const Recording&, std::ostream&) const;

    private:
        std::unique_ptr<Impl> m_impl{};

    public:
        Signal<bool> whenRecognitionChangedSignal;
    };
}

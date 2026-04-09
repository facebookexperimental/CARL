#!/usr/bin/env python3
# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

"""
Compare two CARL test harness JSON result files and produce a markdown report.

Usage:
    python compare_results.py baseline.json current.json --output report.md

No external dependencies — stdlib only.
"""

import argparse
import json
from pathlib import Path


def load_results(path):
    with open(path) as f:
        return json.load(f)


def make_pair_key(result):
    return (result["definition"], result["example"])


def fmt_delta(val, invert=False, prec=4):
    """Format a delta value with direction indicator.
    invert=True means positive delta is bad (e.g., frame time increase, false positive increase).
    invert=False means positive delta is good (e.g., true positive confidence increase).
    """
    if val is None:
        return "N/A"
    sign = "+" if val >= 0 else ""
    formatted = f"{sign}{val:.{prec}f}"
    if abs(val) < 1e-6:
        return formatted
    if invert:
        return f"⚠️ {formatted}" if val > 0 else f"✅ {formatted}"
    else:
        return f"⚠️ {formatted}" if val < 0 else f"✅ {formatted}"


def safe_get(d, *keys, default=None):
    for k in keys:
        if isinstance(d, dict) and k in d:
            d = d[k]
        else:
            return default
    return d


def compute_deltas(baseline_results, current_results):
    baseline_map = {make_pair_key(r): r for r in baseline_results}
    current_map = {make_pair_key(r): r for r in current_results}

    all_keys = set(baseline_map.keys()) | set(current_map.keys())

    associated_deltas = []
    non_associated_deltas = []
    performance_deltas = []
    new_pairs = []
    removed_pairs = []

    for key in sorted(all_keys):
        b = baseline_map.get(key)
        c = current_map.get(key)

        if b is None:
            new_pairs.append(c)
            continue
        if c is None:
            removed_pairs.append(b)
            continue

        perf_delta = {
            "definition": key[0],
            "example": key[1],
            "associated": c.get("associated", False),
            "avgFrameTimeDelta": (
                safe_get(c, "performance", "averageFrameTimeUs", default=0)
                - safe_get(b, "performance", "averageFrameTimeUs", default=0)
            ),
            "p95FrameTimeDelta": (
                safe_get(c, "performance", "p95FrameTimeUs", default=0)
                - safe_get(b, "performance", "p95FrameTimeUs", default=0)
            ),
            "maxFrameTimeDelta": (
                safe_get(c, "performance", "maxFrameTimeUs", default=0)
                - safe_get(b, "performance", "maxFrameTimeUs", default=0)
            ),
            "baseline_avg": safe_get(b, "performance", "averageFrameTimeUs", default=0),
            "current_avg": safe_get(c, "performance", "averageFrameTimeUs", default=0),
        }
        performance_deltas.append(perf_delta)

        if c.get("associated"):
            rec_delta = {
                "definition": key[0],
                "example": key[1],
                "maxInWindowDelta": (
                    safe_get(c, "recognition", "maxConfidenceInWindow", default=0)
                    - safe_get(b, "recognition", "maxConfidenceInWindow", default=0)
                ),
                "maxOutsideWindowDelta": (
                    safe_get(c, "recognition", "maxConfidenceOutsideWindow", default=0)
                    - safe_get(b, "recognition", "maxConfidenceOutsideWindow", default=0)
                ),
                "baseline_inWindow": safe_get(b, "recognition", "maxConfidenceInWindow", default=0),
                "current_inWindow": safe_get(c, "recognition", "maxConfidenceInWindow", default=0),
                "baseline_outsideWindow": safe_get(b, "recognition", "maxConfidenceOutsideWindow", default=0),
                "current_outsideWindow": safe_get(c, "recognition", "maxConfidenceOutsideWindow", default=0),
            }
            associated_deltas.append(rec_delta)
        else:
            rec_delta = {
                "definition": key[0],
                "example": key[1],
                "maxConfidenceDelta": (
                    safe_get(c, "recognition", "maxConfidence", default=0)
                    - safe_get(b, "recognition", "maxConfidence", default=0)
                ),
                "baseline_max": safe_get(b, "recognition", "maxConfidence", default=0),
                "current_max": safe_get(c, "recognition", "maxConfidence", default=0),
            }
            non_associated_deltas.append(rec_delta)

    return {
        "associated": associated_deltas,
        "non_associated": non_associated_deltas,
        "performance": performance_deltas,
        "new_pairs": new_pairs,
        "removed_pairs": removed_pairs,
    }


def short_path(path):
    """Shorten a path for display: take last two path components."""
    parts = Path(path).parts
    return "/".join(parts[-2:]) if len(parts) >= 2 else path


def generate_report(baseline_data, current_data, deltas):
    lines = []
    lines.append("# 🧪 CARL Test Harness Comparison Report\n")

    b_commit = baseline_data.get("commitHash", "unknown")
    c_commit = current_data.get("commitHash", "unknown")
    lines.append(f"**Baseline:** `{b_commit}` → **Current:** `{c_commit}`\n")

    # Summary
    lines.append("## Summary\n")
    n_associated = len(deltas["associated"])
    n_non_associated = len(deltas["non_associated"])
    n_perf = len(deltas["performance"])
    n_new = len(deltas["new_pairs"])
    n_removed = len(deltas["removed_pairs"])

    lines.append(f"| Metric | Count |")
    lines.append(f"|--------|-------|")
    lines.append(f"| Total pairs tested | {n_perf} |")
    lines.append(f"| Associated (true positive) pairs | {n_associated} |")
    lines.append(f"| Non-associated (false positive) pairs | {n_non_associated} |")
    if n_new:
        lines.append(f"| New pairs (not in baseline) | {n_new} |")
    if n_removed:
        lines.append(f"| Removed pairs (not in current) | {n_removed} |")
    lines.append("")

    # Aggregate stats
    if deltas["performance"]:
        avg_deltas = [d["avgFrameTimeDelta"] for d in deltas["performance"]]
        mean_avg_delta = sum(avg_deltas) / len(avg_deltas)
        lines.append(f"**Mean Δ avg frame time:** {fmt_delta(mean_avg_delta, invert=True, prec=1)} μs\n")

    if deltas["associated"]:
        in_window_deltas = [d["maxInWindowDelta"] for d in deltas["associated"]]
        mean_in_window = sum(in_window_deltas) / len(in_window_deltas)
        lines.append(f"**Mean Δ true positive confidence:** {fmt_delta(mean_in_window, invert=False, prec=4)}\n")

    if deltas["non_associated"]:
        fp_deltas = [d["maxConfidenceDelta"] for d in deltas["non_associated"]]
        mean_fp = sum(fp_deltas) / len(fp_deltas)
        lines.append(f"**Mean Δ false positive confidence:** {fmt_delta(mean_fp, invert=True, prec=4)}\n")

    # True Positive Detail
    if deltas["associated"]:
        lines.append("## True Positive Recognition (Associated Pairs)\n")
        lines.append("| Definition | Example | Baseline | Current | Δ In-Window | Δ Outside-Window |")
        lines.append("|------------|---------|----------|---------|-------------|------------------|")
        for d in deltas["associated"]:
            lines.append(
                f"| {short_path(d['definition'])} "
                f"| {short_path(d['example'])} "
                f"| {d['baseline_inWindow']:.4f} "
                f"| {d['current_inWindow']:.4f} "
                f"| {fmt_delta(d['maxInWindowDelta'], invert=False)} "
                f"| {fmt_delta(d['maxOutsideWindowDelta'], invert=True)} |"
            )
        lines.append("")

    # False Positive Detail
    if deltas["non_associated"]:
        lines.append("## False Positive Recognition (Non-Associated Pairs)\n")
        lines.append("| Definition | Example | Baseline | Current | Δ Max Confidence |")
        lines.append("|------------|---------|----------|---------|------------------|")
        for d in deltas["non_associated"]:
            lines.append(
                f"| {short_path(d['definition'])} "
                f"| {short_path(d['example'])} "
                f"| {d['baseline_max']:.4f} "
                f"| {d['current_max']:.4f} "
                f"| {fmt_delta(d['maxConfidenceDelta'], invert=True)} |"
            )
        lines.append("")

    # Performance Detail
    if deltas["performance"]:
        lines.append("## Performance\n")
        lines.append("| Definition | Example | Baseline Avg (μs) | Current Avg (μs) | Δ Avg | Δ P95 |")
        lines.append("|------------|---------|--------------------|--------------------|-------|-------|")
        for d in deltas["performance"]:
            lines.append(
                f"| {short_path(d['definition'])} "
                f"| {short_path(d['example'])} "
                f"| {d['baseline_avg']:.1f} "
                f"| {d['current_avg']:.1f} "
                f"| {fmt_delta(d['avgFrameTimeDelta'], invert=True, prec=1)} "
                f"| {fmt_delta(d['p95FrameTimeDelta'], invert=True, prec=1)} |"
            )
        lines.append("")

    if not deltas["associated"] and not deltas["non_associated"] and not deltas["performance"]:
        lines.append("No matching pairs found between baseline and current results.\n")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Compare CARL test harness results")
    parser.add_argument("baseline", help="Path to baseline results JSON")
    parser.add_argument("current", help="Path to current results JSON")
    parser.add_argument("--output", "-o", default="report.md", help="Output markdown path")
    args = parser.parse_args()

    baseline_data = load_results(args.baseline)
    current_data = load_results(args.current)

    deltas = compute_deltas(
        baseline_data.get("results", []),
        current_data.get("results", []),
    )

    report = generate_report(baseline_data, current_data, deltas)

    with open(args.output, "w") as f:
        f.write(report)

    print(f"Report written to: {args.output}")

    # Print summary to stdout
    print(report)


if __name__ == "__main__":
    main()

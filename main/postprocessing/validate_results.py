"""Compare mean drag, lift RMS, and Strouhal number with saved literature ranges."""
from __future__ import annotations

import csv
import io
import json
import math
import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from run.case_config import ROOT

METRICS = ("mean_CD", "CL_rms", "strouhal")
LABELS = {"mean_CD": "Mean CD", "CL_rms": "CL RMS", "strouhal": "St"}
LITERATURE_PATH = ROOT / "literature/reference_values.csv"


def validate_acceptance(acceptance):
    if (not isinstance(acceptance, dict)
            or acceptance.get("method", "literature_min_max") != "literature_min_max"):
        raise ValueError("acceptance.method must be 'literature_min_max'.")
    unknown = set(acceptance) - {"method", "absolute_tolerances"}
    if unknown:
        raise ValueError(f"Unknown acceptance settings: {sorted(unknown)}")
    tolerances = acceptance.get("absolute_tolerances", {})
    if not isinstance(tolerances, dict):
        raise ValueError("acceptance.absolute_tolerances must be an object.")
    unknown = set(tolerances) - set(METRICS)
    if unknown:
        raise ValueError(f"Unknown acceptance metrics: {sorted(unknown)}")
    for metric, value in tolerances.items():
        if (isinstance(value, bool) or not isinstance(value, (int, float))
                or not math.isfinite(value) or value < 0):
            raise ValueError(f"acceptance.absolute_tolerances.{metric} must be finite and nonnegative.")


def literature_ranges(path=None):
    path = Path(path) if path is not None else LITERATURE_PATH
    raw = path.read_bytes()
    values = {name: [] for name in METRICS}
    excluded_lift = []
    for row in csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))):
        if not row.get("Re", "").strip() or float(row["Re"]) != 100.0:
            continue
        for metric, column in (("mean_CD", "mean_CD"), ("CL_rms", "lift_value"), ("strouhal", "St")):
            text = row.get(column, "").strip()
            if not text:
                continue
            if metric == "CL_rms" and not row.get("lift_definition", "").strip().casefold().startswith("rms"):
                excluded_lift.append(row["source_id"])
                continue
            number = float(text)
            if not math.isfinite(number):
                raise ValueError(f"Nonfinite literature {metric} in {row['source_id']}.")
            values[metric].append((number, row["source_id"]))
    ranges = {}
    for metric, entries in values.items():
        if not entries:
            raise ValueError(f"No comparable Re=100 literature values for {metric}.")
        lower = min(value for value, _ in entries)
        upper = max(value for value, _ in entries)
        ranges[metric] = {
            "lower": lower, "upper": upper,
            "source_ids": ";".join(source for _, source in entries),
            "lower_source_ids": ";".join(source for value, source in entries if value == lower),
            "upper_source_ids": ";".join(source for value, source in entries if value == upper),
        }
    source_name = path.relative_to(ROOT).as_posix() if path.is_relative_to(ROOT) else path.name
    metadata = {
        "method": "inclusive minimum/maximum of saved comparable Re=100 literature values",
        "source_file": source_name,
        "inclusion_policy": "All saved finite mean_CD and St values, including unverified entries; lift_value only for rows labelled RMS.",
        "excluded_non_rms_lift_source_ids": excluded_lift,
        "outside_percent_definition": "100 * abs(value - nearest_bound) / abs(nearest_bound) when outside; zero inside. Undefined if the nearest bound is zero.",
    }
    return ranges, metadata


def evaluate(config, statistics, quality, literature_path=None):
    validate_acceptance(config["acceptance"])
    ranges, metadata = literature_ranges(literature_path)
    tolerances = {metric: config["acceptance"].get("absolute_tolerances", {}).get(metric, 0.0)
                  for metric in METRICS}
    metadata.update(
        absolute_tolerances=tolerances,
        acceptance_rule="lower - absolute_tolerance <= value <= upper + absolute_tolerance",
    )
    rows = []
    for metric in METRICS:
        limits = ranges[metric]
        lower, upper = limits["lower"], limits["upper"]
        tolerance = tolerances[metric]
        acceptance_lower, acceptance_upper = lower - tolerance, upper + tolerance
        value = statistics.get(metric)
        finite = value is not None and math.isfinite(value)
        direction, nearest, absolute, percent = "", None, None, None
        passed = None
        if not quality["passed"] or not finite:
            status = "not_evaluated"
            message = "Not compared: run-quality checks failed or the simulation metric is unavailable."
        elif lower <= value <= upper:
            status = "within_literature_range"
            passed = True
            direction, absolute, percent = "within", 0.0, 0.0
            message = "PASS: Within literature range."
        else:
            passed = acceptance_lower <= value <= acceptance_upper
            status = "within_tolerance" if passed else "outside_tolerance"
            direction = "below" if value < lower else "above"
            nearest = lower if value < lower else upper
            absolute = abs(value - nearest)
            percent = 100 * absolute / abs(nearest) if nearest != 0 else None
            side = "lower" if direction == "below" else "upper"
            exceedance = (f"{percent:.6g}%" if percent is not None else f"{absolute:.6g}")
            decision = "PASS: Outside literature range but within" if passed else "FAIL: Outside"
            message = (f"{decision} absolute tolerance {tolerance:g}; "
                       f"{exceedance} {direction} the literature {side} bound.")
        rows.append({"metric": metric, "value": value if finite else None,
                     "lower": lower, "upper": upper, "status": status,
                     "passed": passed, "absolute_tolerance": tolerance,
                     "acceptance_lower": acceptance_lower, "acceptance_upper": acceptance_upper,
                     "direction": direction, "nearest_bound": nearest,
                     "outside_absolute": absolute, "outside_percent": percent,
                     "message": message, "source_ids": limits["source_ids"],
                     "lower_source_ids": limits["lower_source_ids"],
                     "upper_source_ids": limits["upper_source_ids"]})
    if not quality["passed"]:
        status, code = "invalid_run", 2
    elif any(row["status"] == "not_evaluated" for row in rows):
        status, code = "comparison_unavailable", 3
    elif any(row["status"] == "outside_tolerance" for row in rows):
        status, code = "outside_tolerance", 1
    elif any(row["status"] == "within_tolerance" for row in rows):
        status, code = "within_tolerance", 0
    else:
        status, code = "within_literature_range", 0
    return {"status": status, "exit_code": code, "metrics": rows, "quality": quality,
            "case_type": config["case_type"], "preset": config["preset"],
            "literature_comparison": metadata,
            "note": "Passing requires run quality and literature-range membership or an exceedance within the configured absolute tolerance."}


def save_reports(case_dir, config, statistics, quality, details=None):
    case_dir = Path(case_dir)
    report = evaluate(config, statistics, quality)
    report["statistics"] = statistics
    report.update(details or {})
    case_dir.mkdir(parents=True, exist_ok=True)
    (case_dir / "summary.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    table = case_dir / "tables/metrics.csv"
    table.parent.mkdir(parents=True, exist_ok=True)
    with table.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(report["metrics"][0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(report["metrics"])
    print("\nComparison with literature values\n" + "-" * 72, flush=True)
    for row in report["metrics"]:
        value = f"{row['value']:.8g}" if row["value"] is not None else "unavailable"
        print(f"{LABELS[row['metric']]} = {value}; literature [{row['lower']:g}, {row['upper']:g}]; "
              f"absolute tolerance {row['absolute_tolerance']:g}: {row['message']}", flush=True)
    print("-" * 72, flush=True)
    return report


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_dir", type=Path)
    args = parser.parse_args()
    config = json.loads((args.case_dir / "case_config.json").read_text())
    old = json.loads((args.case_dir / "summary.json").read_text())
    regenerated = {"statistics", "quality", "metrics", "status", "exit_code", "note", "preset", "case_type", "literature_comparison",
                   "surface_pressure_coefficient", "additional_diagnostics"}
    details = {key: value for key, value in old.items()
               if key not in regenerated and not key.endswith("_sha256")}
    report = save_reports(args.case_dir, config, old["statistics"], old["quality"], details)
    print(f"Status: {report['status']}")
    return report["exit_code"]


if __name__ == "__main__":
    raise SystemExit(main())

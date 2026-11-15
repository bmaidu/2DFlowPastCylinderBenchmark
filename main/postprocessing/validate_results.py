"""Compare mean drag and Strouhal number; retain lift RMS as a diagnostic."""
from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from run.case_config import ROOT

# Recognized reference metrics include diagnostics that do not determine PASS/FAIL.
LITERATURE_METRICS = ("mean_CD", "CL_rms", "strouhal")
ACCEPTANCE_METRICS = ("mean_CD", "strouhal")
LABELS = {"mean_CD": "Mean CD", "CL_rms": "CL RMS", "strouhal": "St"}
LITERATURE_PATH = ROOT / "literature/reference_values.csv"


def validate_acceptance(acceptance):
    if (not isinstance(acceptance, dict)
            or acceptance.get("method", "literature_min_max") != "literature_min_max"):
        raise ValueError("acceptance.method must be 'literature_min_max'.")
    unknown = set(acceptance) - {"method"}
    if unknown:
        raise ValueError(f"Unknown acceptance settings: {sorted(unknown)}")


def literature_ranges(path=None):
    path = Path(path) if path is not None else LITERATURE_PATH
    observed = {name: [] for name in LITERATURE_METRICS}
    eligible = {name: [] for name in ACCEPTANCE_METRICS}
    excluded_lift = []
    excluded_acceptance = {name: [] for name in LITERATURE_METRICS}
    study_ids, configuration_ids = set(), set()
    required = {"source_id", "study_id", "Re", "dimension", "body_motion", "cylinder_count",
                "mean_CD", "St", "lift_definition", "primary_source_verified",
                "verified_metrics", "observed_eligible", "acceptance_metrics", "eligibility_reason"}
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Literature CSV missing columns: {sorted(required - set(reader.fieldnames or []))}")
        for row in reader:
            source_id = row["source_id"].strip()
            if not source_id or source_id in configuration_ids:
                raise ValueError(f"Empty or duplicate literature source_id: {source_id!r}.")
            configuration_ids.add(source_id)
            if not row["study_id"].strip():
                raise ValueError(f"Missing study_id in {source_id}.")
            study_ids.add(row["study_id"])
            verified = set(filter(None, row["verified_metrics"].split(";")))
            acceptance = set(filter(None, row["acceptance_metrics"].split(";")))
            if (verified - set(LITERATURE_METRICS)
                    or acceptance - set(ACCEPTANCE_METRICS)
                    or not acceptance.issubset(verified)):
                raise ValueError(f"Unknown or unverified acceptance metric in {source_id}.")
            for flag in ("primary_source_verified", "observed_eligible"):
                if row[flag] not in ("true", "false"):
                    raise ValueError(f"{flag} must be true or false in {source_id}.")
            comparable = (float(row["Re"]) == 100 and row["dimension"] == "2"
                          and row["body_motion"] == "stationary" and row["cylinder_count"] == "1"
                          and row["primary_source_verified"] == "true"
                          and row["observed_eligible"] == "true")
            if acceptance and (not comparable or not row["eligibility_reason"].strip()):
                raise ValueError(f"Acceptance entry lacks comparable verified context or eligibility reason: {source_id}.")
            if row.get("lift_value", "").strip() and not row["lift_definition"].casefold().startswith("rms"):
                excluded_lift.append(source_id)
            for metric, column in (("mean_CD", "mean_CD"), ("CL_rms", "CL_RMS"), ("strouhal", "St")):
                # Diagnostic reference columns may be omitted from the literature CSV.
                raw = row.get(column, "").strip()
                if not raw:
                    if metric in verified:
                        raise ValueError(f"Verified {metric} is missing in {source_id}.")
                    continue
                number = float(raw)
                if not math.isfinite(number):
                    raise ValueError(f"Nonfinite literature {metric} in {source_id}.")
                if metric == "CL_rms" and metric in verified and not row["lift_definition"].casefold().startswith("rms"):
                    raise ValueError(f"Verified CL_RMS requires an explicit RMS definition in {source_id}.")
                if not comparable or metric not in verified:
                    continue
                entry = (number, source_id, row["study_id"])
                observed[metric].append(entry)
                if metric in acceptance:
                    eligible[metric].append(entry)
                else:
                    excluded_acceptance[metric].append({"source_id": source_id, "reason": row["eligibility_reason"]})

    def bounds(entries):
        lower = min((value for value, _, _ in entries), default=None)
        upper = max((value for value, _, _ in entries), default=None)
        return {
            "lower": lower, "upper": upper,
            "source_ids": ";".join(source for _, source, _ in entries),
            "lower_source_ids": ";".join(source for value, source, _ in entries if value == lower),
            "upper_source_ids": ";".join(source for value, source, _ in entries if value == upper),
            "configuration_count": len(entries),
            "study_ids": sorted({study for _, _, study in entries}),
            "study_count": len({study for _, _, study in entries}),
        }

    ranges = {metric: bounds(entries) for metric, entries in eligible.items()}
    source_name = path.relative_to(ROOT).as_posix() if path.is_relative_to(ROOT) else path.name
    metadata = {
        "method": "inclusive minimum/maximum of verified, eligible Re=100 reference configurations",
        "source_file": source_name,
        "inclusion_policy": "Primary-verified stationary single-cylinder 2D Re=100 entries; only mean_CD and strouhal in acceptance_metrics determine PASS/FAIL. Explicit CL_RMS values remain diagnostic references; lift amplitude/extrema are never converted. See literature/README.md.",
        "acceptance_metrics": list(ACCEPTANCE_METRICS),
        "database_study_count": len(study_ids),
        "database_configuration_count": len(configuration_ids),
        "verified_numeric_configuration_count": len({source for entries in observed.values() for _, source, _ in entries}),
        "policy_version": "selected_seven_drag_strouhal_v2",
        "observed_ranges": {metric: bounds(entries) for metric, entries in observed.items()},
        "eligible_ranges": ranges,
        "excluded_acceptance": excluded_acceptance,
        "excluded_non_rms_lift_source_ids": excluded_lift,
        "outside_percent_definition": "100 * abs(value - nearest_bound) / abs(nearest_bound) when outside; zero inside. Undefined if the nearest bound is zero.",
    }
    return ranges, metadata


def evaluate(config, statistics, quality, literature_path=None):
    validate_acceptance(config["acceptance"])
    ranges, metadata = literature_ranges(literature_path)
    metadata["acceptance_rule"] = "lower <= value <= upper"
    rows = []
    for metric in ACCEPTANCE_METRICS:
        limits = ranges[metric]
        lower, upper = limits["lower"], limits["upper"]
        available = lower is not None and upper is not None
        value = statistics.get(metric)
        finite = value is not None and math.isfinite(value)
        direction, nearest, absolute, percent = "", None, None, None
        passed = None
        if not quality["passed"] or not finite:
            status = "not_evaluated"
            message = "Not compared: run-quality checks failed or the simulation metric is unavailable."
        elif not available:
            status = "not_evaluated"
            message = "Not compared: no verified, eligible literature values for this metric."
        elif lower <= value <= upper:
            status = "within_literature_range"
            passed = True
            direction, absolute, percent = "within", 0.0, 0.0
            message = "PASS: Within eligible literature range."
        else:
            passed = False
            status = "outside_literature_range"
            direction = "below" if value < lower else "above"
            nearest = lower if value < lower else upper
            absolute = abs(value - nearest)
            percent = 100 * absolute / abs(nearest) if nearest != 0 else None
            side = "lower" if direction == "below" else "upper"
            exceedance = (f"{percent:.6g}%" if percent is not None else f"{absolute:.6g}")
            message = (f"FAIL: Outside eligible literature range; "
                       f"{exceedance} {direction} the literature {side} bound.")
        rows.append({"metric": metric, "value": value if finite else None,
                     "lower": lower, "upper": upper, "status": status,
                     "passed": passed,
                     "observed_lower": metadata["observed_ranges"][metric]["lower"],
                     "observed_upper": metadata["observed_ranges"][metric]["upper"],
                     "reference_study_count": limits["study_count"],
                     "reference_configuration_count": limits["configuration_count"],
                     "direction": direction, "nearest_bound": nearest,
                     "outside_absolute": absolute, "outside_percent": percent,
                     "message": message, "source_ids": limits["source_ids"],
                     "lower_source_ids": limits["lower_source_ids"],
                     "upper_source_ids": limits["upper_source_ids"]})
    if not quality["passed"]:
        status, code = "invalid_run", 2
    elif any(row["status"] == "not_evaluated" for row in rows):
        status, code = "comparison_unavailable", 3
    elif any(row["status"] == "outside_literature_range" for row in rows):
        status, code = "outside_literature_range", 1
    else:
        status, code = "within_literature_range", 0
    return {"status": status, "exit_code": code, "metrics": rows, "quality": quality,
            "case_type": config["case_type"], "preset": config["preset"],
            "literature_comparison": metadata,
            "note": "PASS requires run quality and inclusive eligible literature ranges for mean_CD and strouhal only. CL_rms remains a simulation diagnostic. It is not a solver-baseline regression comparison or a mesh-convergence certificate."}


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
        limits = (f"[{row['lower']:g}, {row['upper']:g}]" if row["lower"] is not None else "unavailable")
        print(f"{LABELS[row['metric']]} = {value}; eligible literature {limits}: {row['message']}", flush=True)
    print("-" * 72, flush=True)
    return report


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_dir", type=Path)
    args = parser.parse_args()
    config = json.loads((args.case_dir / "case_config.json").read_text())
    old = json.loads((args.case_dir / "summary.json").read_text())
    regenerated = {"statistics", "quality", "metrics", "status", "exit_code", "note", "preset", "case_type", "literature_comparison"}
    details = {key: value for key, value in old.items()
               if key not in regenerated and not key.endswith("_sha256")}
    report = save_reports(args.case_dir, config, old["statistics"], old["quality"], details)
    print(f"Status: {report['status']}")
    return report["exit_code"]


if __name__ == "__main__":
    raise SystemExit(main())

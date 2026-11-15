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

ACCEPTANCE_METRICS = ("mean_CD", "strouhal")
LABELS = {"mean_CD": "Mean CD", "strouhal": "St"}
METRIC_COLUMNS = ("metric", "value", "status")
RUN_DETAILS = ("n_snapshots", "missing_steps", "unexpected_steps", "rejected_snapshots",
               "time_warnings", "force_operator", "input_type")
STATUS_MESSAGES = {
    "within_literature_range": "PASS: Within literature range.",
    "outside_literature_range": "FAIL: Outside literature range.",
    "not_evaluated": "Not evaluated.",
}
LITERATURE_PATH = ROOT / "literature/reference_values.csv"


def validate_acceptance(acceptance):
    if (not isinstance(acceptance, dict)
            or acceptance.get("method", "literature_min_max") != "literature_min_max"):
        raise ValueError("acceptance.method must be 'literature_min_max'.")
    unknown = set(acceptance) - {"method"}
    if unknown:
        raise ValueError(f"Unknown acceptance settings: {sorted(unknown)}")


def literature_ranges(path=None):
    """Use every selected paper's reported mean drag and St as reference values."""
    path = Path(path) if path is not None else LITERATURE_PATH
    observed = {name: [] for name in ACCEPTANCE_METRICS}
    papers = set()
    required = {"paper", "mean_CD", "St"}
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Literature CSV missing columns: {sorted(required - set(reader.fieldnames or []))}")
        for row in reader:
            paper = (row["paper"] or "").strip()
            if not paper or paper in papers:
                raise ValueError(f"Empty or duplicate literature paper: {paper!r}.")
            papers.add(paper)
            for metric, column in (("mean_CD", "mean_CD"), ("strouhal", "St")):
                raw = (row[column] or "").strip()
                try:
                    number = float(raw)
                except ValueError as error:
                    raise ValueError(f"Missing or invalid literature {column} in {paper}.") from error
                if not math.isfinite(number):
                    raise ValueError(f"Nonfinite literature {column} in {paper}.")
                observed[metric].append(number)

    return {
        metric: {"lower": min(values, default=None), "upper": max(values, default=None)}
        for metric, values in observed.items()
    }


def evaluate(config, statistics, quality, literature_path=None):
    validate_acceptance(config["acceptance"])
    ranges = literature_ranges(literature_path)
    rows = []
    for metric in ACCEPTANCE_METRICS:
        limits = ranges[metric]
        lower, upper = limits["lower"], limits["upper"]
        available = lower is not None and upper is not None
        value = statistics.get(metric)
        finite = value is not None and math.isfinite(value)
        if not quality["passed"] or not finite or not available:
            status = "not_evaluated"
        elif lower <= value <= upper:
            status = "within_literature_range"
        else:
            status = "outside_literature_range"
        rows.append({"metric": metric, "value": value if finite else None, "status": status})
    if not quality["passed"]:
        status, code = "invalid_run", 2
    elif any(row["status"] == "not_evaluated" for row in rows):
        status, code = "comparison_unavailable", 3
    elif any(row["status"] == "outside_literature_range" for row in rows):
        status, code = "outside_literature_range", 1
    else:
        status, code = "within_literature_range", 0
    return {"status": status, "exit_code": code, "metrics": rows, "quality": quality,
            "case_type": config["case_type"], "preset": config["preset"]}


def save_reports(case_dir, config, statistics, quality, details=None):
    case_dir = Path(case_dir)
    report = evaluate(config, statistics, quality)
    report["statistics"] = statistics
    # Copy only run diagnostics so refreshing old reports cannot restore retired fields.
    report.update({key: value for key, value in (details or {}).items() if key in RUN_DETAILS})
    case_dir.mkdir(parents=True, exist_ok=True)
    (case_dir / "summary.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    table = case_dir / "tables/metrics.csv"
    table.parent.mkdir(parents=True, exist_ok=True)
    with table.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=METRIC_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(report["metrics"])
    print("\nComparison with literature values\n" + "-" * 72, flush=True)
    for row in report["metrics"]:
        value = f"{row['value']:.8g}" if row["value"] is not None else "unavailable"
        print(f"{LABELS[row['metric']]} = {value}; {STATUS_MESSAGES[row['status']]}", flush=True)
    print("-" * 72, flush=True)
    return report


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_dir", type=Path)
    args = parser.parse_args()
    config = json.loads((args.case_dir / "case_config.json").read_text())
    old = json.loads((args.case_dir / "summary.json").read_text())
    details = {key: old[key] for key in RUN_DETAILS if key in old}
    report = save_reports(args.case_dir, config, old["statistics"], old["quality"], details)
    print(f"Status: {report['status']}")
    return report["exit_code"]


if __name__ == "__main__":
    raise SystemExit(main())

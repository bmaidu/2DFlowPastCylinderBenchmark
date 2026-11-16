from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from .validate_results import save_reports

from . import plotting
from .force_integration import (coefficient_scale, direct_viscous_force, element_vorticity,
    pressure_force, validate_force_operator, velocity_gradients)
from .validate_2d_fields import (discover_results, expected_steps, load_case_geometry,
    read_snapshot, read_solver_config)
from .cycle_statistics import cycle_window_statistics, trapezoidal_mean


def write_csv(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def summarize_history(config, rows, quality_failures=None):
    """Shared CSV/field statistics; time is physical time (D=U=1)."""
    t = np.asarray([row["time"] for row in rows])
    drag = np.asarray([row["CD_total"] for row in rows])
    lift = np.asarray([row["CL_total"] for row in rows])
    failures = list(quality_failures or [])
    if len(t) < 3 or not np.all(np.isfinite(np.column_stack((t, drag, lift)))):
        failures.append("Insufficient or nonfinite force-history samples.")
        cycle = {"status": "insufficient cycles", "mean_CD": None, "mean_CL": None,
                 "CL_rms": None, "strouhal": None, "cycles_used": 0, "t_start": None, "t_end": None}
    elif np.any(np.diff(t) <= 0):
        raise ValueError("Physical times must be strictly increasing.")
    else:
        expected_interval = config["solver"]["dt"] * config["solver"]["vtk_increment"]
        if not np.allclose(np.diff(t), expected_interval, rtol=0, atol=1e-7):
            failures.append("Missing or irregular saved samples; gaps cannot be bridged for validation.")
        cycle = cycle_window_statistics(t, drag, lift, diameter=1.0, speed=1.0,
                                         criteria=config["statistics"])
        if cycle["status"] != "developed":
            failures.extend(cycle.get("reasons", []) or ["Developed periodic shedding is not established."])
        if cycle.get("t_start") is not None:
            window = (t >= cycle["t_start"] - 1e-8) & (t <= cycle["t_end"] + 1e-8)
            for column, name in (("CD_pressure", "mean_CD_pressure"),
                                  ("CD_viscous", "mean_CD_viscous")):
                if all(column in row for row in rows):
                    cycle[name] = trapezoidal_mean(t[window], np.asarray([row[column] for row in rows])[window])
    return cycle, {"passed": not failures, "failures": failures}


def field_figures(case_dir, prepared, velocity, pressure, physical_time, step):
    points, triangles = prepared["points"], prepared["triangles"]
    gradients = velocity_gradients(velocity, triangles, prepared["gradients"])
    x0, x1 = float(points[:, 0].min()), float(points[:, 0].max())
    y0, y1 = float(points[:, 1].min()), float(points[:, 1].max())
    zoom = ((max(x0, -2.0), min(x1, 15.0)), (max(y0, -3.0), min(y1, 3.0)))
    centroid = points[triangles].mean(axis=1)
    node_mask = ((points[:, 0] >= zoom[0][0]) & (points[:, 0] <= zoom[0][1])
                 & (points[:, 1] >= zoom[1][0]) & (points[:, 1] <= zoom[1][1]))
    cell_mask = ((centroid[:, 0] >= zoom[0][0]) & (centroid[:, 0] <= zoom[0][1])
                 & (centroid[:, 1] >= zoom[1][0]) & (centroid[:, 1] <= zoom[1][1]))
    vorticity = element_vorticity(gradients)
    items = (("ux", velocity[:, 0], node_mask, "streamwise velocity", "u_x", None),
             ("pressure", pressure, node_mask, "pressure", "p", None),
             ("vorticity", vorticity, cell_mask, "vorticity", "omega_z", (-5.0, 5.0)))
    for name, values, mask, title, label, limits in items:
        if limits is None:
            sample = values[mask] if np.any(mask) else values
            lo, hi = np.percentile(sample, [0.5, 99.5])
            if hi <= lo:
                hi = lo + 1e-12
            limits = (float(lo), float(hi))
        plotting.plot_scalar_field(points, triangles, values,
            case_dir / "figures" / f"{name}_step_{step}_t{physical_time:.3f}.png",
            f"{title}, t={physical_time:.3f}", label, limits, "coolwarm", zoom)


def postprocess(case_dir, config, *, solver_completed=False):
    case_dir = Path(case_dir)
    # An interrupted or failed analysis must not leave an old successful report.
    case_dir.mkdir(parents=True, exist_ok=True)
    (case_dir / "summary.json").write_text(
        json.dumps({"status": "processing", "exit_code": 2}, indent=2, allow_nan=False) + "\n")
    try:
        return _postprocess(case_dir, config, solver_completed=solver_completed)
    except Exception as error:
        save_reports(case_dir, config, {key: None for key in
            ("mean_CD", "mean_CL", "CL_rms", "strouhal")},
            {"passed": False, "failures": [f"Postprocessing failed: {error}"]},
            {"input_type": "raw solver fields"})
        raise


def _postprocess(case_dir, config, *, solver_completed=False):
    case_dir = Path(case_dir)
    for name in ("tables", "figures"):
        (case_dir / name).mkdir(exist_ok=True)
    for folder, pattern in (
        ("figures", "cylinder_2d_summary.png"), ("figures", "cylinder_2d_summary.pdf"),
        ("figures", "speed_step_*.png"), ("figures", "wall_shear_last_snapshot.png"),
        ("figures", "ux_step_*.png"), ("figures", "pressure_step_*.png"),
        ("figures", "vorticity_step_*.png"),
        ("tables", "wall_shear_last_snapshot.csv"),
    ):
        for path in (case_dir / folder).glob(pattern):
            path.unlink()
    stable_seconds = 0.0 if solver_completed else 15.0
    solver = read_solver_config(case_dir / "solver.xml")
    prepared = load_case_geometry(case_dir, solver)
    operator = validate_force_operator(prepared)
    points = np.column_stack((prepared["points"], prepared["coordinate_z"]))
    triangles = prepared["triangles"]
    records, rejected = discover_results(case_dir / solver["results_dir"], solver["vtk_prefix"],
                                         solver["dt"], stable_seconds, None)
    failures = []
    if not operator["passed"]:
        failures.append("Force-operator validation failed.")
    if rejected:
        failures.append("Some snapshots were rejected or are still being written.")
    missing = sorted(set(expected_steps(solver)) - {row["step"] for row in records})
    unexpected = sorted({row["step"] for row in records} - set(expected_steps(solver)))
    if missing:
        failures.append(f"Run is missing {len(missing)} expected snapshots, including any unsaved final step.")
    if unexpected:
        failures.append("Snapshots outside the configured step schedule were found.")
    rows, warnings = [], []
    scale = coefficient_scale(solver["density"], config["solver"]["speed"], 1.0)
    last = None
    for index, record in enumerate(records):
        velocity, pressure, time_value, _ = read_snapshot(record["path"], points, triangles)
        if time_value is None:
            time_value = record["filename_time"]
            warnings.append({"step": record["step"], "action": "time reconstructed from step*dt"})
        elif abs(time_value - record["filename_time"]) > 1e-6:
            failures.append(f"Physical time disagrees with step*dt at step {record['step']}.")
        geometry = prepared["boundaries"]["cylinder"]
        p_force = pressure_force(pressure, geometry) / scale
        v_force = direct_viscous_force(velocity, geometry, solver["viscosity"]) / scale
        rows.append({"step": record["step"], "time": time_value,
            "CD_pressure": float(p_force[0]), "CD_viscous": float(v_force[0]),
            "CD_total": float(p_force[0]+v_force[0]), "CL_pressure": float(p_force[1]),
            "CL_viscous": float(v_force[1]), "CL_total": float(p_force[1]+v_force[1])})
        last = (velocity, pressure, time_value, record["step"])
        if (index+1) % 100 == 0 or index+1 == len(records):
            print(f"{index+1}/{len(records)} snapshots analyzed", flush=True)
    columns = ("step", "time", "CD_pressure", "CD_viscous", "CD_total", "CL_pressure", "CL_viscous", "CL_total")
    write_csv(case_dir / "tables/force_history.csv", rows, columns)
    cycle, quality = summarize_history(config, rows, failures)
    if rows:
        plotting.plot_summary_dashboard(case_dir / "figures/force_history",
            {"time": [row["time"] for row in rows], "cd_total": [row["CD_total"] for row in rows],
             "cl_total": [row["CL_total"] for row in rows], "final_time": rows[-1]["time"],
             "cycle_statistics": cycle})
        field_figures(case_dir, prepared, *last)
    print("\n" + "=" * 72 + "\nPostprocessing is finished.\n" + "=" * 72, flush=True)
    report = save_reports(case_dir, config, cycle, quality,
        {"n_snapshots": len(rows), "missing_steps": missing, "unexpected_steps": unexpected,
         "rejected_snapshots": rejected, "time_warnings": warnings,
         "force_operator": operator,
         "input_type": "raw solver fields"})
    print(f"Reports saved in {case_dir}\n", flush=True)
    return report["exit_code"]

"""Figures for the 2D cylinder postprocessor."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.ticker import MaxNLocator
from matplotlib.tri import Triangulation
import numpy as np


PLOT_CONFIG = {
    "cd_ylim": None,
    "cl_ylim": None,
    "axis_reference_tmin": 0.10,
    "display_tmin": 0.0,
    "display_tmax": None,
    "axis_padding_fraction": 0.10,
}


def axis_limits(
    time,
    series,
    *,
    manual=None,
    reference_tmin=0.10,
    display_tmin=0.0,
    display_tmax=None,
    padding=0.10,
    symmetric=False,
):
    """Y limits from post-reference samples. Manual limits are not padded.

    Samples before reference_tmin are omitted from the automatic range only.
    They remain available to the caller for plotting.
    """
    time = np.asarray(time, dtype=float).reshape(-1)
    columns = [np.asarray(values, dtype=float).reshape(-1) for values in series]
    if any(len(column) != len(time) for column in columns):
        raise ValueError("A plotted series does not match the time array.")
    finite = np.isfinite(time) & np.all([np.isfinite(column) for column in columns], axis=0)
    usable_time = time[finite]
    usable = [column[finite] for column in columns]
    x_max = float(np.nanmax(usable_time)) if len(usable_time) else 0.0
    if display_tmax is None:
        display_tmax = x_max
    x_min = 0.0 if display_tmin is None else float(display_tmin)
    shown = (usable_time >= x_min) & (usable_time <= float(display_tmax))
    reference = shown & (usable_time >= float(reference_tmin))
    if np.count_nonzero(reference) < 2:
        reference = shown
        reference_note = "too few samples after the reference time; displayed samples were used"
    else:
        reference_note = "samples at or after the reference time"
    if np.count_nonzero(reference) == 0:
        reference = finite
        reference_note = "no displayed samples; all finite samples were used"

    if manual is not None:
        low, high = float(manual[0]), float(manual[1])
        mode = "manual"
    else:
        chosen = np.concatenate([column[reference] for column in usable]) if np.count_nonzero(reference) else np.array([0.0])
        if len(chosen) == 0 or not np.any(np.isfinite(chosen)):
            low, high, mode = -1.0, 1.0, "empty"
        else:
            data_low = float(np.min(chosen))
            data_high = float(np.max(chosen))
            if symmetric:
                bound = max(abs(data_low), abs(data_high))
                span = bound if bound > 0.0 else 1.0
                pad = float(padding) * span
                low, high = -(bound + pad), bound + pad
                mode = "symmetric"
            else:
                span = data_high - data_low
                if span <= 0.0:
                    span = max(abs(data_high), 1.0)
                    mode = "constant"
                else:
                    mode = "auto"
                pad = float(padding) * span
                low, high = data_low - pad, data_high + pad
    if not np.isfinite(low) or not np.isfinite(high) or high <= low:
        low, high, mode = -1.0, 1.0, "fallback"

    outside = 0
    outside_early = 0
    for column in usable:
        hidden = shown & ((column < low) | (column > high))
        outside += int(np.count_nonzero(hidden))
        outside_early += int(np.count_nonzero(hidden & (usable_time < float(reference_tmin))))
    return {
        "ylim": (low, high),
        "xlim": (x_min, float(display_tmax)),
        "mode": mode,
        "reference_note": reference_note,
        "reference_count": int(np.count_nonzero(reference)),
        "outside_count": outside,
        "outside_before_reference": outside_early,
    }


def _save(figure, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(figure)


def _triangulation(points, triangles):
    return Triangulation(points[:, 0], points[:, 1], triangles)


def plot_scalar_field(points, triangles, values, path, title, label, limits, cmap, zoom):
    figure, axis = plt.subplots(figsize=(8.6, 4.6))
    triangulation = _triangulation(points, triangles)
    shading = "flat" if len(values) == len(triangles) else "gouraud"
    image = axis.tripcolor(
        triangulation,
        values,
        shading=shading,
        cmap=cmap,
        vmin=limits[0],
        vmax=limits[1],
    )
    axis.add_patch(Circle((0.0, 0.0), 0.5, facecolor="white", edgecolor="black", lw=0.8, zorder=3))
    axis.set_aspect("equal")
    axis.set_xlim(*zoom[0])
    axis.set_ylim(*zoom[1])
    axis.set_xlabel("x/D")
    axis.set_ylabel("y/D")
    axis.set_title(title)
    figure.colorbar(image, ax=axis, label=label, fraction=0.046, pad=0.04)
    figure.tight_layout()
    _save(figure, path)


def _style_time_axis(axis, xlim):
    axis.set_xlim(*xlim)
    axis.set_xlabel(r"Nondimensional time, $t^*$")
    axis.grid(True, color="0.88", linewidth=0.6)
    axis.xaxis.set_major_locator(MaxNLocator(nbins=6, steps=[1, 2, 2.5, 5, 10]))
    axis.yaxis.set_major_locator(MaxNLocator(nbins=5, steps=[1, 2, 2.5, 5, 10]))
    axis.tick_params(labelsize=10)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)


def _format_value(value, digits=3):
    if value is None or not np.isfinite(value):
        return "n/a"
    magnitude = abs(float(value))
    if magnitude != 0.0 and (magnitude < 1.0e-2 or magnitude >= 1.0e3):
        return f"{value:.2e}"
    return f"{value:.{digits}f}"


def plot_summary_dashboard(path_stem, data, config=None):
    """Two time histories and three cycle metrics, saved as PNG."""
    settings = dict(PLOT_CONFIG)
    if config:
        settings.update(config)
    time = np.asarray(data["time"], dtype=float)
    total_drag = np.asarray(data["cd_total"], dtype=float)
    total_lift = np.asarray(data["cl_total"], dtype=float)
    drag_limits = axis_limits(
        time,
        (total_drag,),
        manual=settings["cd_ylim"],
        reference_tmin=settings["axis_reference_tmin"],
        display_tmin=settings["display_tmin"],
        display_tmax=settings["display_tmax"],
        padding=settings["axis_padding_fraction"],
        symmetric=False,
    )
    lift_limits = axis_limits(
        time,
        (total_lift,),
        manual=settings["cl_ylim"],
        reference_tmin=settings["axis_reference_tmin"],
        display_tmin=settings["display_tmin"],
        display_tmax=settings["display_tmax"],
        padding=settings["axis_padding_fraction"],
        symmetric=settings["cl_ylim"] is None,
    )
    cycle = data.get("cycle_statistics") or {"status": "insufficient cycles", "t_start": None, "t_end": None}

    with plt.rc_context(
        {
            "font.size": 11,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.facecolor": "white",
        }
    ):
        figure = plt.figure(figsize=(13.0, 5.6), constrained_layout=True)
        grid = figure.add_gridspec(2, 2, height_ratios=(1.0, 0.30))
        drag_axis = figure.add_subplot(grid[0, 0])
        lift_axis = figure.add_subplot(grid[0, 1])
        strip = figure.add_subplot(grid[1, :])
        _draw_drag(drag_axis, time, total_drag, drag_limits, cycle)
        _draw_lift(lift_axis, time, total_lift, lift_limits, cycle)
        _draw_cycle_strip(strip, cycle)
        figure.suptitle(
            f"2D flow past cylinder\nFinal t={data['final_time']:.3g}",
            fontsize=14,
        )
        stem = Path(path_stem)
        stem.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(stem.with_suffix(".png"), dpi=200)
        plt.close(figure)
    return {"drag": drag_limits, "lift": lift_limits}


def _shade_cycles(axis, cycle):
    start = cycle.get("t_start")
    end = cycle.get("t_end")
    if start is None or end is None or cycle.get("status") == "insufficient cycles":
        return
    axis.axvspan(start, end, facecolor="#d6e4f0", edgecolor="none", alpha=0.85, zorder=0)


def _style_summary_axis(axis, limits, ylabel, title):
    axis.set_ylim(*limits["ylim"])
    _style_time_axis(axis, limits["xlim"])
    axis.set_xlabel(r"$t$")
    axis.set_ylabel(ylabel)
    axis.set_title(title)
    axis.set_frame_on(True)
    for side in ("top", "right", "bottom", "left"):
        axis.spines[side].set_visible(True)


def _draw_drag(axis, time, total, limits, cycle):
    _shade_cycles(axis, cycle)
    axis.plot(time, total, color="black", lw=2.05, zorder=3)
    _style_summary_axis(axis, limits, r"$C_D$", r"Drag coefficient $C_D$")


def _draw_lift(axis, time, total, limits, cycle):
    _shade_cycles(axis, cycle)
    axis.axhline(0.0, color="0.72", lw=0.8, zorder=1)
    axis.plot(time, total, color="#215f9a", lw=2.05, zorder=3)
    _style_summary_axis(axis, limits, r"$C_L$", r"Lift coefficient $C_L$")
    if limits["outside_count"]:
        axis.text(
            0.01,
            0.03,
            f"{limits['outside_count']} plotted samples lie outside this axis",
            transform=axis.transAxes,
            ha="left",
            va="bottom",
            fontsize=8,
            color="0.25",
        )


def _draw_cycle_strip(axis, cycle):
    axis.set_axis_off()
    axis.set_xlim(0.0, 1.0)
    axis.set_ylim(0.0, 1.0)
    axis.plot([0.0, 1.0], [0.96, 0.96], color="0.82", lw=0.8, transform=axis.transAxes, clip_on=False)
    status = cycle.get("status") or "insufficient cycles"
    if status == "insufficient cycles" or cycle.get("t_start") is None:
        found = cycle.get("available_complete_cycles")
        requested = cycle.get("requested_cycles")
        detail = "No complete-cycle averaging interval is available."
        if found is not None and requested is not None:
            detail = (
                f"Found {found} of {requested} requested complete cycles. "
                "No shorter window was substituted."
            )
        axis.text(0.5, 0.62, "Insufficient complete cycles", ha="center", va="center", fontsize=12)
        axis.text(0.5, 0.28, detail, ha="center", va="center", fontsize=10, color="0.25")
        return

    strouhal = _format_value(cycle.get("strouhal"), 3)
    if status == "provisional" and cycle.get("strouhal") is not None:
        strouhal = f"{strouhal} (provisional)"
    elif cycle.get("strouhal") is None:
        strouhal = "n/a"
    values = (
        f"Mean $C_D$ = {_format_value(cycle.get('mean_CD'))}"
        f"      $C_{{L,\\mathrm{{rms}}}}$ = {_format_value(cycle.get('CL_rms'))}"
        f"      St = {strouhal}"
    )
    axis.text(0.5, 0.55, values, ha="center", va="center", fontsize=12)
    if status == "provisional":
        pending = []
        if cycle.get("amplitude_stabilized") is False and not any(
            "increasing" in reason for reason in (cycle.get("reasons") or [])
        ):
            pending.append("lift amplitude")
        if cycle.get("period_stabilized") is False:
            pending.append("period")
        if cycle.get("drag_stabilized") is False:
            pending.append("cycle-mean drag")
        increasing = any("increasing" in reason for reason in (cycle.get("reasons") or []))
        parts = []
        if increasing:
            parts.append("Lift amplitude is increasing")
        if pending:
            parts.append("not stabilized: " + ", ".join(pending))
        if not parts and cycle.get("reasons"):
            parts.append(cycle["reasons"][0])
        if parts:
            axis.text(
                0.5,
                0.08,
                ". ".join(parts),
                ha="center",
                va="center",
                fontsize=9,
                color="0.30",
            )

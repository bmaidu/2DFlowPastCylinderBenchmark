"""Qualified lift cycles and time-averaged force statistics."""

from __future__ import annotations

import numpy as np


def _distance_filter(indices, time, signal, separation, kind):
    kept = []
    for index in indices:
        if not kept or time[index] - time[kept[-1]] >= separation:
            kept.append(index)
            continue
        if kind == "max" and signal[index] > signal[kept[-1]]:
            kept[-1] = index
        if kind == "min" and signal[index] < signal[kept[-1]]:
            kept[-1] = index
    return kept


CYCLE_CRITERIA = {
    "n_cycles": 3,
    "min_separation": 2.0,
    "prominence_fraction": 0.25,
    "minimum_prominence": 0.01,
    "max_period_relative_range": 0.05,
    "max_amplitude_relative_range": 0.10,
    "max_drag_relative_range": 0.02,
    "diameter": 1.0,
    "speed": 1.0,
}

CYCLE_DEFINITION = (
    "A complete cycle is the closed interval between consecutive qualified lift "
    "maxima that contains at least one qualified lift minimum. Samples after the "
    "last qualified maximum are not a complete cycle, and the simulation end time "
    "does not close a cycle by itself."
)


def trapezoidal_integral(time, values):
    """Integrate with the trapezoidal rule. Time may be nonuniform."""
    time = np.asarray(time, dtype=float)
    values = np.asarray(values, dtype=float)
    if len(time) != len(values):
        raise ValueError("Time and value arrays have different lengths.")
    if len(time) < 2:
        raise ValueError("At least two samples are required for integration.")
    if np.any(~np.isfinite(time)) or np.any(~np.isfinite(values)):
        raise ValueError("Integration received a non-finite sample.")
    if np.any(np.diff(time) <= 0.0):
        raise ValueError("Integration time must be strictly increasing.")
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(values, time))
    return float(np.sum(0.5 * (values[1:] + values[:-1]) * np.diff(time)))


def trapezoidal_mean(time, values):
    time = np.asarray(time, dtype=float)
    duration = float(time[-1] - time[0])
    if duration <= 0.0:
        raise ValueError("The averaging interval has zero duration.")
    return trapezoidal_integral(time, values) / duration


def _local_extrema(signal):
    maxima = []
    minima = []
    for index in range(1, len(signal) - 1):
        if signal[index] > signal[index - 1] and signal[index] >= signal[index + 1]:
            maxima.append(index)
        elif signal[index] < signal[index - 1] and signal[index] <= signal[index + 1]:
            minima.append(index)
    return maxima, minima


def _side_excursion(index, signal, opposite, kind):
    """Smaller rise above the available neighboring opposite extrema.

    A peak at the end of the record is judged from the preceding trough only.
    An interior peak must clear both neighboring troughs.
    """
    earlier = [item for item in opposite if item < index]
    later = [item for item in opposite if item > index]
    rises = []
    if earlier:
        neighbor = signal[earlier[-1]]
        rises.append(signal[index] - neighbor if kind == "max" else neighbor - signal[index])
    if later:
        neighbor = signal[later[0]]
        rises.append(signal[index] - neighbor if kind == "max" else neighbor - signal[index])
    if not rises:
        return 0.0
    return float(min(rises))


def _relative_range(values):
    values = np.asarray(values, dtype=float)
    if len(values) < 2 or not np.all(np.isfinite(values)):
        return None
    span = float(np.max(values) - np.min(values))
    scale = abs(float(np.mean(values)))
    if scale == 0.0:
        return 0.0 if span == 0.0 else None
    return span / scale


def _cycle_record(time, drag, lift, start, end):
    segment_time = time[start : end + 1]
    segment_drag = drag[start : end + 1]
    segment_lift = lift[start : end + 1]
    amplitude = 0.5 * (float(np.max(segment_lift)) - float(np.min(segment_lift)))
    return {
        "start_index": int(start),
        "end_index": int(end),
        "t_start": float(segment_time[0]),
        "t_end": float(segment_time[-1]),
        "period": float(segment_time[-1] - segment_time[0]),
        "mean_CD": trapezoidal_mean(segment_time, segment_drag),
        "mean_CL": trapezoidal_mean(segment_time, segment_lift),
        "lift_amplitude": float(amplitude),
        "selected": False,
    }


def cycle_window_statistics(time, drag, lift, n_cycles=None, diameter=None, speed=None, criteria=None):
    """Statistics on the last N complete lift cycles.

    Means and the lift RMS are trapezoidal time integrals over the selected
    interval. They are not averages of every saved sample from startup.
    Fewer cycles are never substituted for a missing requested count.
    """
    settings = dict(CYCLE_CRITERIA)
    if criteria:
        settings.update(criteria)
    if n_cycles is not None:
        settings["n_cycles"] = int(n_cycles)
    if diameter is not None:
        settings["diameter"] = float(diameter)
    if speed is not None:
        settings["speed"] = float(speed)
    requested = int(settings["n_cycles"])
    if requested < 1:
        raise ValueError("The requested cycle count must be at least 1.")

    time = np.asarray(time, dtype=float)
    drag = np.asarray(drag, dtype=float)
    lift = np.asarray(lift, dtype=float)
    if not (len(time) == len(drag) == len(lift)):
        raise ValueError("Time, drag, and lift arrays have different lengths.")

    raw_maxima, raw_minima = _local_extrema(lift)
    maxima = _distance_filter(raw_maxima, time, lift, settings["min_separation"], "max")
    minima = _distance_filter(raw_minima, time, lift, settings["min_separation"], "min")
    maximum_excursions = [_side_excursion(index, lift, minima, "max") for index in maxima]
    positive = [value for value in maximum_excursions if value > 0.0]
    scale = float(np.median(positive)) if positive else 0.0
    prominence = max(float(settings["minimum_prominence"]), float(settings["prominence_fraction"]) * scale)
    qualified_maxima = [
        index for index, height in zip(maxima, maximum_excursions) if height >= prominence
    ]
    minimum_excursions = [_side_excursion(index, lift, maxima, "min") for index in minima]
    qualified_minima = [
        index for index, height in zip(minima, minimum_excursions) if height >= prominence
    ]

    detected = []
    for start, end in zip(qualified_maxima, qualified_maxima[1:]):
        if not any(start < index < end for index in qualified_minima):
            continue
        segment = time[start : end + 1]
        if len(segment) < 2 or segment[-1] <= segment[0]:
            continue
        if not np.all(np.isfinite(segment)) or not np.all(np.isfinite(drag[start : end + 1])):
            continue
        if not np.all(np.isfinite(lift[start : end + 1])):
            continue
        detected.append(_cycle_record(time, drag, lift, start, end))

    result = {
        "definition": CYCLE_DEFINITION,
        "integration": "trapezoidal",
        "requested_cycles": requested,
        "available_complete_cycles": len(detected),
        "cycles_used": 0,
        "fallback_applied": False,
        "fallback_note": None,
        "t_start": None,
        "t_end": None,
        "status": "insufficient cycles",
        "mean_CD": None,
        "mean_CL": None,
        "CL_rms": None,
        "mean_period": None,
        "strouhal": None,
        "strouhal_is_benchmark": False,
        "amplitude_stabilized": False,
        "period_stabilized": False,
        "drag_stabilized": False,
        "period_relative_range": None,
        "amplitude_relative_range": None,
        "drag_relative_range": None,
        "reasons": [],
        "qualified_maxima": [
            {"time": float(time[index]), "CL": float(lift[index])} for index in qualified_maxima
        ],
        "qualified_minima": [
            {"time": float(time[index]), "CL": float(lift[index])} for index in qualified_minima
        ],
        "detected_cycles": detected,
        "selected_cycles": [],
        "criteria": {
            "n_cycles": requested,
            "min_separation": float(settings["min_separation"]),
            "prominence_fraction": float(settings["prominence_fraction"]),
            "minimum_prominence": float(settings["minimum_prominence"]),
            "prominence_used": float(prominence),
            "excursion_scale": scale,
            "max_period_relative_range": float(settings["max_period_relative_range"]),
            "max_amplitude_relative_range": float(settings["max_amplitude_relative_range"]),
            "max_drag_relative_range": float(settings["max_drag_relative_range"]),
            "diameter": float(settings["diameter"]),
            "speed": float(settings["speed"]),
        },
    }
    if len(detected) < requested:
        note = (
            f"Requested {requested} complete cycles, found {len(detected)}. "
            "Fewer cycles were not substituted into the reported average."
        )
        result["fallback_note"] = note
        result["reasons"] = [note]
        return result

    selected = detected[-requested:]
    for item in selected:
        item["selected"] = True
    start = selected[0]["start_index"]
    end = selected[-1]["end_index"]
    window_time = time[start : end + 1]
    window_drag = drag[start : end + 1]
    window_lift = lift[start : end + 1]
    mean_lift = trapezoidal_mean(window_time, window_lift)
    periods = [item["period"] for item in selected]
    amplitudes = [item["lift_amplitude"] for item in selected]
    cycle_drag = [item["mean_CD"] for item in selected]
    period_range = _relative_range(periods)
    amplitude_range = _relative_range(amplitudes)
    drag_range = _relative_range(cycle_drag)
    period_ok = period_range is not None and period_range <= settings["max_period_relative_range"]
    amplitude_ok = (
        amplitude_range is not None and amplitude_range <= settings["max_amplitude_relative_range"]
    )
    drag_ok = drag_range is not None and drag_range <= settings["max_drag_relative_range"]
    reasons = []
    if not amplitude_ok:
        if len(amplitudes) >= 2 and np.all(np.diff(amplitudes) > 0.0):
            reasons.append("lift amplitude is increasing across the selected cycles")
        else:
            reasons.append("lift amplitude is not stable across the selected cycles")
    if not period_ok:
        reasons.append("shedding period is not stable across the selected cycles")
    if not drag_ok:
        reasons.append("cycle-averaged drag is not stable across the selected cycles")
    mean_period = float(np.mean(periods))
    strouhal = None
    if mean_period > 0.0 and settings["speed"] != 0.0:
        strouhal = float(settings["diameter"] / (settings["speed"] * mean_period))
    developed = not reasons and strouhal is not None
    result.update(
        {
            "cycles_used": requested,
            "t_start": float(window_time[0]),
            "t_end": float(window_time[-1]),
            "status": "developed" if developed else "provisional",
            "mean_CD": trapezoidal_mean(window_time, window_drag),
            "mean_CL": mean_lift,
            "CL_rms": float(np.sqrt(trapezoidal_mean(window_time, (window_lift - mean_lift) ** 2))),
            "mean_period": mean_period,
            "strouhal": strouhal,
            "strouhal_is_benchmark": developed,
            "amplitude_stabilized": amplitude_ok,
            "period_stabilized": period_ok,
            "drag_stabilized": drag_ok,
            "period_relative_range": period_range,
            "amplitude_relative_range": amplitude_range,
            "drag_relative_range": drag_range,
            "reasons": reasons,
            "selected_cycles": selected,
        }
    )
    return result

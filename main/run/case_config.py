from __future__ import annotations

import copy
import json
import math
from pathlib import Path

from mesh.mesh_parameters import Domain, MeshParameters, SIZE_FIELDS

ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = ROOT / "benchmark_default_config.json"


def case_name(config):
    x0, x1 = config["geometry"]["domain_x"]
    y0, y1 = config["geometry"]["domain_y"]
    return f"case_{x1-x0:g}D_{y1-y0:g}D_hwall{config['mesh']['wall_size']:.3f}"


def validate(config):
    g, m, s = (config[key] for key in ("geometry", "mesh", "solver"))
    if g["diameter"] != 1.0 or g["center"] != [0.0, 0.0]:
        raise ValueError("This benchmark fixes D=1 and the cylinder center at (0,0).")
    for key in ("domain_x", "domain_y"):
        if not isinstance(g[key], (list, tuple)) or len(g[key]) != 2:
            raise ValueError(f"{key} must contain exactly two bounds.")
    Domain(*g["domain_x"], *g["domain_y"])
    MeshParameters(**m)
    for key in ("density", "viscosity", "speed", "dt"):
        if not isinstance(s[key], (int, float)) or not math.isfinite(s[key]) or s[key] <= 0:
            raise ValueError(f"solver.{key} must be finite and positive.")
    if not math.isclose(s["density"] * s["speed"] / s["viscosity"], 100.0):
        raise ValueError("The configuration must retain Re=100.")
    if any(s[key] != value for key, value in (("density", 1.0), ("speed", 1.0), ("viscosity", 0.01))):
        raise ValueError("This benchmark retains density=1, speed=1, and dynamic viscosity=0.01.")
    for key in ("n_steps", "vtk_increment", "restart_increment"):
        if not isinstance(s[key], int) or s[key] <= 0:
            raise ValueError(f"solver.{key} must be a positive integer.")
    if s["n_steps"] % s["vtk_increment"]:
        raise ValueError("n_steps must be divisible by vtk_increment so the final field is saved.")
    criteria = config["statistics"]
    if not isinstance(criteria["n_cycles"], int) or criteria["n_cycles"] < 3:
        raise ValueError("At least three complete cycles are required.")
    for key, value in criteria.items():
        if not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
            raise ValueError(f"statistics.{key} must be finite and positive.")
    from postprocessing.validate_results import validate_acceptance
    validate_acceptance(config["acceptance"])


def resolve(preset=None, overrides=None, config_file=None):
    base = json.loads(CONFIG_PATH.read_text())
    overrides = {key: value for key, value in (overrides or {}).items() if value is not None}
    custom = bool(overrides) or config_file is not None
    if preset and custom:
        raise ValueError("Preset names cannot be combined with domain/mesh overrides or --config. "
                         "Omit coarse/medium/fine to define a custom case.")
    if preset is not None and preset not in base["presets"]:
        raise ValueError(f"Unknown preset: {preset}")
    selected = preset or base["default_case"]
    result = {"schema_version": 1, "case_type": "custom" if custom else "preset",
              "preset": None if custom else selected,
              "geometry": copy.deepcopy(base["geometry"]),
              "mesh": copy.deepcopy(base["presets"][selected]),
              "solver": copy.deepcopy(base["solver"]),
              "statistics": copy.deepcopy(base["statistics"]),
              "acceptance": copy.deepcopy(base["acceptance"])}
    if config_file:
        patch = json.loads(Path(config_file).read_text())
        allowed = {"geometry", "mesh", "solver", "statistics", "acceptance", "schema_version",
                   "case_type", "preset"}
        unknown = set(patch) - allowed
        if unknown:
            raise ValueError(f"Unknown configuration sections: {sorted(unknown)}")
        for section in ("geometry", "mesh", "solver", "statistics"):
            if section not in patch:
                continue
            if not isinstance(patch[section], dict):
                raise ValueError(f"{section} must be an object.")
            unknown = set(patch[section]) - set(result[section])
            if unknown:
                raise ValueError(f"Unknown {section} settings: {sorted(unknown)}")
            result[section].update(patch[section])
        if "acceptance" in patch:
            result["acceptance"] = patch["acceptance"]
    for key, value in overrides.items():
        if key in ("domain_x", "domain_y"):
            result["geometry"][key] = list(value)
        elif key in result["mesh"]:
            result["mesh"][key] = value
        else:
            raise ValueError(f"Unknown override: {key}")
    validate(result)
    return result

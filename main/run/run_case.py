#!/usr/bin/env python3
"""Prepare, run, and analyze a preset or user-defined case."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
import traceback

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from run.case_config import CONFIG_PATH, ROOT, SIZE_FIELDS, case_name, resolve, validate
from run.create_solver import TEMPLATE, create_solver
from mesh.mesh_parameters import BOUNDARY_NAMES, Domain, MeshParameters
from run.run_svMP import executable, final_step_reached, run_solver, sha256


def parser():
    base = json.loads(CONFIG_PATH.read_text())
    p = argparse.ArgumentParser(description=__doc__,
        epilog=f"For custom cases, omit the preset name; unspecified settings inherit "
               f"{base['default_case']} values. Use --prepare-only to generate inputs without running CFD.")
    p.add_argument("preset", nargs="?", choices=tuple(base["presets"]),
                   help=f"Preset case; without arguments {base['default_case']} is selected.")
    for axis, bounds in (("x", ("XMIN", "XMAX")), ("y", ("YMIN", "YMAX"))):
        limits = ", ".join(f"{value:g}" for value in base["geometry"][f"domain_{axis}"])
        p.add_argument(f"--domain-{axis}", nargs=2, type=float, metavar=bounds,
                       help=f"{axis.upper()} bounds in D. All presets: [{limits}].")
    for field in (*SIZE_FIELDS, "distance_sampling"):
        values = "; ".join(f"{name}={mesh[field]:g}" for name, mesh in base["presets"].items())
        description = ("Wall-distance sampling points" if field == "distance_sampling"
                       else "Custom mesh setting in D")
        p.add_argument("--" + field.replace("_", "-"),
                       type=int if field == "distance_sampling" else float, default=None,
                       help=f"{description}. Presets: {values}.")
    p.add_argument("--config", type=Path, help="Custom JSON settings; CLI geometry/mesh options take precedence.")
    p.add_argument("--results-dir", type=Path, default=ROOT / "results",
                   help="Case output directory (default: %(default)s).")
    location = p.add_mutually_exclusive_group()
    location.add_argument("--case-name", help="Results subfolder name; omitted names use domain dimensions and wall size.")
    location.add_argument("--case-dir", type=Path, help="Existing prepared case for --postprocess-only.")
    mode = p.add_mutually_exclusive_group()
    mode.add_argument("--prepare-only", action="store_true", help="Generate mesh and solver input, without CFD.")
    mode.add_argument("--postprocess-only", action="store_true", help="Read an existing prepared case; never run CFD.")
    p.add_argument("--solver", default="svmultiphysics",
                   help="Executable name on PATH, or executable path (default: %(default)s).")
    p.add_argument("--mpi-exec", default="mpirun",
                   help="MPI executable name or path (default: %(default)s).")
    p.add_argument("--np", type=int, default=8, help="MPI process count (default: %(default)s).")
    return p


def input_files(case_dir):
    return [case_dir / "solver.xml", case_dir / "mesh/mesh-complete.mesh.vtu"] + [
        case_dir / "mesh/mesh-surfaces" / f"{name}.vtp"
        for name in BOUNDARY_NAMES]


def prepare(case_dir, config):
    """Reuse existing inputs and refresh mesh figures; create missing inputs as needed."""
    from mesh.generate_mesh import generate
    from mesh.plot_mesh import plot_mesh_figures
    import gmsh
    import numpy
    import vtk

    case_dir = Path(case_dir)
    config_path = case_dir / "case_config.json"
    if config_path.is_file():
        config = json.loads(config_path.read_text())
        print(f"Using existing case inputs and saved settings in {case_dir}", flush=True)
    case_dir.mkdir(parents=True, exist_ok=True)
    if not config_path.exists():
        config_path.write_text(json.dumps(config, indent=2, allow_nan=False) + "\n")
    mesh_files = input_files(case_dir)[1:]
    generated = not all(path.is_file() for path in mesh_files)
    provenance_path = case_dir / "provenance.json"
    provenance = json.loads(provenance_path.read_text()) if provenance_path.is_file() else {}
    try:
        if generated:
            g = config["geometry"]
            print(f"Generating mesh in {case_dir / 'mesh'}; see mesh/mesh_generation.log", flush=True)
            generate(MeshParameters(**config["mesh"]), case_dir, Domain(*g["domain_x"], *g["domain_y"]))
            provenance.update(
                created_utc=datetime.now(timezone.utc).isoformat(), python=Path(sys.executable).name,
                python_version=sys.version,
                dependencies={"gmsh": gmsh.__version__, "numpy": numpy.__version__,
                              "vtk": vtk.vtkVersion.GetVTKVersion()},
                source_sha256={str(path.relative_to(ROOT)): sha256(path) for path in
                    (ROOT / "main/mesh/generate_mesh.py", ROOT / "main/mesh/mesh_parameters.py",
                     ROOT / "main/mesh/plot_mesh.py", TEMPLATE, ROOT / "main/run/create_solver.py")})
        else:
            plot_mesh_figures(case_dir / "mesh")
        if not (case_dir / "solver.xml").is_file():
            create_solver(case_dir, config)
        for name in ("solver_output", "figures", "tables"):
            (case_dir / name).mkdir(exist_ok=True)
        if generated:
            provenance["input_sha256"] = {
                str(path.relative_to(case_dir)): sha256(path) for path in input_files(case_dir)}
            provenance["preparation_status"] = "completed"
            provenance_path.write_text(json.dumps(provenance, indent=2, allow_nan=False) + "\n")
        return config
    except BaseException as error:
        if generated:
            provenance.update(preparation_status="failed", error=str(error))
            provenance_path.write_text(json.dumps(provenance, indent=2, allow_nan=False) + "\n")
        raise


def main(argv=None):
    p = parser()
    args = p.parse_args(argv)
    overrides = {name: getattr(args, name) for name in
                 ("domain_x", "domain_y", *SIZE_FIELDS, "distance_sampling")}
    try:
        if args.np < 1:
            raise ValueError("--np must be positive.")
        if args.case_name is not None and (
                not args.case_name.strip() or args.case_name in (".", "..")
                or any(character in args.case_name for character in ("/", "\\", "\0"))):
            raise ValueError("--case-name must be a nonempty folder name without path separators.")
        if args.case_dir:
            if not args.postprocess_only:
                raise ValueError("--case-dir is only valid with --postprocess-only.")
            if args.preset or args.config or any(value is not None for value in overrides.values()):
                raise ValueError("--case-dir uses its saved settings; omit preset/config/mesh overrides.")
            case_dir = args.case_dir.expanduser().resolve()
            config = json.loads((case_dir / "case_config.json").read_text())
            validate(config)
        else:
            config = resolve(args.preset, overrides, args.config)
            name = args.case_name if args.case_name is not None else case_name(config)
            case_dir = args.results_dir.expanduser().resolve() / name
        if not args.prepare_only:
            from postprocessing.postprocess_case import postprocess
        if args.postprocess_only:
            saved = case_dir / "case_config.json"
            if saved.is_file():
                config = json.loads(saved.read_text())
            from mesh.plot_mesh import plot_mesh_figures
            plot_mesh_figures(case_dir / "mesh")
            message = "Starting postprocessing of saved CFD results."
        else:
            if not args.prepare_only and not final_step_reached(case_dir):
                executable(args.solver)
                executable(args.mpi_exec)
            config = prepare(case_dir, config)
            if args.prepare_only:
                print(f"Prepared {case_dir}")
                return 0
            ran_solver = run_solver(case_dir, args.solver, args.mpi_exec, args.np)
            message = ("CFD run is finished. Now starting postprocessing." if ran_solver else
                       "Using existing CFD results. Now starting postprocessing.")
        print("\n" + "=" * 72 + f"\n{message}\n" + "=" * 72 + "\n", flush=True)
        return postprocess(case_dir, config)
    except (ValueError, OSError, RuntimeError, ImportError) as error:
        print(f"Error: {error}", file=sys.stderr)
        if isinstance(error, ImportError):
            print("Install requirements.txt into the Python environment used for this command.", file=sys.stderr)
        return 2
    except Exception:
        traceback.print_exc()
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

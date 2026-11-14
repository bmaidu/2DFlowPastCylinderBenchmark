# 2D flow past cylinder at Re=100: svMultiPhysics benchmark

This benchmark generates a triangular mesh, runs 2D incompressible flow past cylinder at Re = 100, and compares mean drag, lift RMS, and Strouhal number with literature values. The default configuration supports coarse, medium, and fine presets, plus custom domains and regional mesh sizes.

## Problem setup

The default domain size and boundary conditions are shown below. The cylinder is centered at the origin, and lengths are normalized by its diameter $D$. Here, $u$ and $v$ are the streamwise ($x$) and transverse velocity ($y$) components, respectively.

![Default cylinder-flow domain and boundary conditions](docs/figures/problem_setup.svg)

The meshes use linear triangular (TRI3) elements generated with Gmsh's Frontal-Delaunay algorithm (6). The default simulations use a time step of $0.005$ for $30,000$ steps, reaching $t = 150$. 

## Requirements

- Python 3.10 or newer and the packages in [requirements.txt](requirements.txt).
- svMultiPhysics available as `svmultiphysics` on `PATH`, or supplied with `--solver`.
- An MPI launcher available as `mpirun` on `PATH`, or supplied with `--mpi-exec`.

Install dependencies in your Python environment:

```bash
python -m pip install -r requirements.txt
```

Run all commands below from the repository root using that environment's Python.

## Default presets

| Preset | Nodes | Elements | Cylinder wall | Near-wall outer | Far upstream | Upstream | Approach | Initial wake | Middle wake | Extended wake | Transition | Outlet | Far field |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| coarse | 23,837 | 47,327 | 0.012375 | 0.18 | 1.53 | 0.81 | 0.18 | 0.063 | 0.0855 | 0.12375 | 0.2475 | 0.495 | 1.575 |
| medium | 54,050 | 107,582 | 0.00825 | 0.12 | 1.02 | 0.54 | 0.12 | 0.042 | 0.057 | 0.0825 | 0.165 | 0.33 | 1.05 |
| fine | 120,353 | 239,933 | 0.0055 | 0.08 | 0.68 | 0.36 | 0.08 | 0.028 | 0.038 | 0.055 | 0.11 | 0.22 | 0.70 |

All regional sizes are target element sizes $h/D$. Preset cases are saved as `preset_coarse`, `preset_medium`, and `preset_fine`.

Mesh targets are defined in [benchmark_default_config.json](benchmark_default_config.json); solver settings are in [solver_template.xml](main/run/solver_template.xml).

## Running cases

Run a preset; coarse is the default if no preset name is provided:

```bash
python main/run/run_case.py coarse
python main/run/run_case.py medium
python main/run/run_case.py fine
```

Specify executable paths and the MPI process count when needed (`--np` defaults to 8):

```bash
python main/run/run_case.py coarse \
  --solver /path/to/svmultiphysics --mpi-exec /path/to/mpirun --np 8
```

For a custom case, omit the preset name and specify target element size in a region. Unspecified settings inherit coarse defaults:

```bash
python main/run/run_case.py \
  --domain-x -40 20 --domain-y -30 30 \
  --wall-size 0.01 --initial-wake-size 0.05 --case-name custom_case
```

Lengths are in cylinder diameters. Expanding the domain keeps refinement locations fixed and extends the endpoint size targets. `--far-size` also sets the global size maximum; other element-size targets must not exceed it. Use `python main/run/run_case.py --help` for all options and preset values.

A custom JSON file can also provide settings:

```json
{
  "geometry": {"domain_x": [-40, 20], "domain_y": [-30, 30]},
  "mesh": {"wall_size": 0.01, "initial_wake_size": 0.05}
}
```

```bash
python main/run/run_case.py --config /path/to/custom.json
```

Settings take precedence in this order: coarse defaults, custom JSON, then explicit command-line domain/mesh values.

Generate inputs without CFD, or postprocess saved results:

```bash
python main/run/run_case.py fine --prepare-only
python main/run/run_case.py --postprocess-only --case-name custom_case
python main/run/run_case.py --postprocess-only --case-dir results/<case-directory>
```

Existing cases reuse their saved settings and inputs. Choose a new `--case-name` or `--results-dir` when changing settings for an independent run. Automatic custom names encode domain dimensions and wall size, for example `custom_40D_24D_hwall0.012`.

The solver is skipped when a nonempty .vtu reaches the configured final step. Otherwise, the saved `solver.xml` controls restarting: `Continue_previous_simulation=true` resumes from restart data; `false` starts from scratch.

## Outputs

```text
benchmark_default_config.json
main/
  run/
  mesh/
  postprocessing/
literature/
results/<case-directory>/
  case_config.json
  solver.xml
  mesh/
    figures/
  solver_output/             # .vtu, restart data, solver history
  figures/
  tables/
  summary.json
  run.log
```

Paths below are relative to each case directory:

| Output | Contents |
|---|---|
| `mesh/figures/` | Full-domain, near-cylinder, and near/middle-wake mesh views |
| `figures/force_history.png` | Drag and lift histories with mean drag, lift RMS, and St |
| `figures/ux_step_*.png`, `figures/pressure_step_*.png`, `figures/vorticity_step_*.png` | Final flow fields; vorticity color limits [-5, 5] |
| `tables/force_history.csv` | Drag and lift coefficients, including pressure and viscous contributions |
| `tables/cycle_statistics.json` | Selected cycles, statistics, and stability checks, including mean lift |
| `tables/metrics.csv` | Observed/eligible literature bounds, passing intervals, source counts and percentage exceedances |
| `summary.json` | Combined metrics and run-quality status |
| `run.log` | Solver output, also printed live in the terminal |

## Validation

Statistics use the last three qualified complete lift cycles and trapezoidal time averages. Drag and lift forces include pressure and viscous contributions per unit span.

| Metric | Definition |
|---|---|
| Mean drag | Mean of `CD = 2*Fx/(rho*U^2*D)` |
| Lift RMS | `sqrt(mean((CL - mean(CL))^2))`, with `CL = 2*Fy/(rho*U^2*D)` |
| Strouhal number | `St = D/(U*mean_period)` |

Literature comparisons require complete, consistent snapshots and developed periodic shedding. Statistical criteria are defined in `benchmark_default_config.json`.

Each metric is compared with a primary-verified literature envelope in [literature/reference_values.csv](literature/reference_values.csv). It supports both preset and user-defined meshes. Reports distinguish the full observed span from acceptance: author-identified confinement and numerical sensitivity tests can remain in the database without determining passing limits. Lift amplitude, extrema and peak-to-peak values are not treated as RMS. Unverified claims are retained separately.

| Metric | Eligible literature envelope | Additional absolute tolerance | Passing interval |
|---|---|---:|---|
| Mean CD | 1.300–1.3917 | 0.02 | 1.280–1.4117 |
| CL RMS | 0.215–0.237 | 0.02 | 0.195–0.257 |
| St | 0.162–0.171 | 0.002 | 0.160–0.173 |

A metric passes when `lower - absolute_tolerance <= value <= upper + absolute_tolerance`. Eligibility and source locations are documented in [literature/README.md](literature/README.md); bounds are recalculated from the CSV. Margins are set in `acceptance.absolute_tolerances` and saved per case. Existing cases retain their saved margins; an omitted margin is zero.

PASS means meeting a broad literature validation criterion and the existing run-quality checks. It does not certify mesh convergence or agreement with a previous solver baseline. Reports identify values within the eligible envelope, outside it but within tolerance, or failing. Outside percentages use `100 * abs(value - nearest_bound) / abs(nearest_bound)`. Exit codes are 0 for passing, 1 for outside tolerance, 2 for invalid run quality and 3 for an unavailable comparison.

Refresh literature comparisons from saved statistics:

```bash
python main/postprocessing/validate_results.py results/<case-directory>
```

This refresh updates only `summary.json` and `tables/metrics.csv`, preserving saved metrics and analysis windows. It does not rerun CFD or rewrite force histories.

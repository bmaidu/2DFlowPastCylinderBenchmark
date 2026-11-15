# 2D flow past cylinder at Re=100: svMultiPhysics benchmark

This benchmark generates a triangular mesh, runs 2D incompressible flow past cylinder at Re = 100, and compares mean drag and Strouhal number with literature values. The default configuration supports coarse, medium, and fine presets, plus custom domains and regional mesh sizes.

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
| `figures/force_history.png` | Drag and lift histories with mean drag and lift coefficients, and Stouhal number |
| `figures/ux_step_*.png`, `figures/pressure_step_*.png`, `figures/vorticity_step_*.png` | Flow fields visualization at final time step |
| `tables/force_history.csv` | Drag and lift coefficients |
| `tables/cycle_statistics.json` | Selected cycles, statistics, and stability checks |
| `tables/metrics.csv` | literature bounds and other metrics |
| `summary.json` | Combined metrics and run status |
| `run.log` | Solver output, also printed live in the terminal |

## Validation

2D flow past a circular cylinder at $Re=100$ is a well-established benchmark in CFD. Here, we assess svMultiPhysics against published studies with comparable domain sizes, mesh resolutions, and numerical methods. Pointwise agreement between flow fields is not expected, the mean drag coefficient ($\overline{C_D}$) and Strouhal number ($St$) serve as the primary comparison metrics and determine the PASS/FAIL outcome in svbenchmark.

| Paper | Domain size (D) | $\overline{C_D}$ | $St$ | Numerical method |
| --- | --- | --- | --- | --- |
| [Tumkur et al. (2013)](https://doi.org/10.1016/j.jfluidstructs.2013.03.008) | 48 × 24 | 1.360 | 0.168 | RB-VMS FEM |
| [Kang & Masud (2021)](https://doi.org/10.1007/s11012-020-01227-w) | 48 × 24 | 1.374 | 0.168 | Immersed-boundary VMS FEM |
| [Kamensky et al. (2015)](https://doi.org/10.1016/j.cma.2014.10.040) | 40 × 20 | 1.381 | 0.171 | Immersogeometric VMS FEM |
| [Behr et al. (1995)](https://doi.org/10.1016/0045-7825(94)00736-7) | 30.5 × 24 | 1.3721 | 0.1624 | Space-time FEM |
| [Tsetoglou et al. (2025)](https://doi.org/10.1002/fld.5334) | 65 × 30 | 1.335 | 0.167 | Body-fitted FVM |
| [Deng et al. (2019)](https://doi.org/10.1103/PhysRevE.100.053107) | 46 × 32 | 1.374 | 0.164 | Immersed-boundary FVM |
| [Alam et al. (2017)](https://doi.org/10.1016/j.jfluidstructs.2017.01.014) | 35 × 20 | 1.372 | 0.165 | FVM |

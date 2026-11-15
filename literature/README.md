# Re=100 cylinder literature benchmark

These seven publications were selected primarily for comparable computational
domain sizes and/or numerical methods, particularly the variational multiscale
finite-element methods related to svMultiPhysics's residual-based VMS (RB-VMS FEM).
Only stationary, isolated-cylinder validation results for two-dimensional,
incompressible flow at Re=100 are used, including those within broader FSI,
moving-body or multiple-cylinder studies. One configuration is retained per paper.

Domain size below is total streamwise length × total height, normalized by diameter D.

| Paper | Domain size (D) | Mean CD | St | Numerical method |
| --- | --- | --- | --- | --- |
| [Tumkur et al. (2013)](https://doi.org/10.1016/j.jfluidstructs.2013.03.008) | 48 × 24 | 1.360 | 0.168 | RB-VMS FEM |
| [Kang & Masud (2021)](https://doi.org/10.1007/s11012-020-01227-w) | 48 × 24 | 1.374 | 0.168 | Immersed-boundary VMS FEM |
| [Kamensky et al. (2015)](https://doi.org/10.1016/j.cma.2014.10.040) | 40 × 20 | 1.381 | 0.171 | Immersogeometric VMS FEM |
| [Behr et al. (1995)](https://doi.org/10.1016/0045-7825(94)00736-7) | 30.5 × 24 | 1.3721 | 0.1624 | Space-time FEM |
| [Tsetoglou et al. (2025)](https://doi.org/10.1002/fld.5334) | 65 × 30 | 1.335 | 0.167 | Body-fitted FVM |
| [Deng et al. (2019)](https://doi.org/10.1103/PhysRevE.100.053107) | 46 × 32 | 1.374 | 0.164 | Immersed-boundary FVM |
| [Alam et al. (2017)](https://doi.org/10.1016/j.jfluidstructs.2017.01.014) | 35 × 20 | 1.372 | 0.165 | FVM |

[reference_values.csv](reference_values.csv) preserves configuration identifiers,
bibliographic metadata, original-source locations and lift definitions.
Tumkur's 48D length comprises 12D upstream and 36D downstream; its two-dimensional
flow is represented using one spanwise layer of brick elements.

PASS requires the existing run-quality checks to pass and **both mean drag and St**
to lie within the inclusive minimum/maximum bounds derived dynamically from
verified, eligible CSV records: **mean CD [1.335, 1.381]** and **St [0.1624, 0.171]**.
All seven CSV rows use `mean_CD;strouhal` for both `verified_metrics` and
`acceptance_metrics`.

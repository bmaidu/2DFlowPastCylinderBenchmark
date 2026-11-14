# Re=100 cylinder-flow literature

[reference_values.csv](reference_values.csv) records original-source results for stationary, two-dimensional, single-cylinder flow. It contains **71 numerical configurations from 15 studies**, plus one verified setup/period-only record from Liu et al. (2023). [unverified_candidates.csv](unverified_candidates.csv) preserves 13 unresolved numerical-claim or retrieval records, including 10 additional studies. Across both files, there are **26 distinct publications**. These counts describe the documented search, not an exhaustive review of all Re=100 publications.

The three comparison metrics are mean drag, RMS of lift fluctuations, and Strouhal number. Mean lift remains in simulation statistics. Different domains, boundary conditions, discretizations, time steps and averaging procedures can produce different coefficients; a literature envelope is not an experimental uncertainty interval or evidence of mesh convergence.

## Data and verification

Each row identifies one reported configuration. The table records authors, year, title, journal, DOI/URL, physical context, method, geometry in cylinder diameters, boundary conditions, mesh, time step and units, coefficient values, original location, verification and eligibility. Empty cells mean unreported or unconfirmed, never zero. Distances are measured from the cylinder center unless the row notes state otherwise. Approximate counts remain approximate in the notes.

`source_id` identifies a configuration; `study_id` groups a publication. Repeated appearances of a result in another table or another paper do not create another configuration. In particular:

- Qu D10 appears in Tables 2 and 3; it retains the existing ID `qu2013_table3`. D9's more precise lift RMS is taken from Table 2.
- Alam's selected 69,000-node, dt=0.01 configuration appears in Tables 1 and 2; it retains `alam2017_table2`.
- Kamensky's repeated l=3 result and boundary-fitted reference are each counted once across Tables 1, 2 and 4.
- Tsetoglou's finest VOS-IVP case appears in Tables 2 and 3 and is counted once. The 2024 online publication belongs to the 2025 journal issue.
- Deng's former `deng2019_candidate` is now the verified bare-cylinder record `deng2019_bare`.

`primary_source_verified=true` means the original article or author manuscript was inspected. `verified_metrics` lists only independently checked principal metrics. Deng and Tumkur numerical passages were checked in indexed **original PDF text**; direct PDF downloads were blocked/rate-limited. The other numerical records were checked in accessible original text/PDFs. This distinction and any unresolved setup details are recorded explicitly. No plotted force curve was digitized to invent precise drag or RMS values.

`CL_RMS` is used only when RMS is explicitly reported. `lift_value` and `lift_definition` preserve amplitude, extrema, maximum or peak-to-peak lift separately. They are never divided by sqrt(2), halved, or otherwise silently converted. Park's original Table 1, Stalberg's inspected author version and TF-73 report amplitude/maxima; their older secondary RMS claims remain in the candidate table. Posdziech and Mittal (2005) also remain there until their original numerical tables and definitions can be checked. Their previously saved values have not been discarded or represented as newly verified RMS measurements.

## Observed range and acceptance policy

The **observed range** spans all finite, directly verified comparable metric values in the recorded stationary single-cylinder Re=100 tests. It includes documented domain, mesh, time-step and quadrature sensitivity points. Some are intentionally constrained or demonstrate numerical error; their inclusion describes what was published, not what a well-resolved unbounded-flow calculation should reproduce.

The **eligible literature envelope** uses metric-specific `acceptance_metrics` with the following rules:

1. Require primary verification, Re=100, dimension=2, one stationary cylinder, a compatible coefficient definition and an explicit eligibility reason. Other Reynolds numbers, interacting cylinders, moving-body cases and secondary comparison rows do not contribute.
2. Numerical method, mesh count or domain differences from svMultiPhysics do not by themselves exclude a reference. Accept reported validation configurations and the resolved/refined subsets supported by each study's own sensitivity discussion.
3. Retain author-identified confinement, under-resolution and under-integration tests in the observed data, but exclude their affected metrics from acceptance. The per-study rules below apply uniformly to their configurations, irrespective of the current svMultiPhysics values.
4. Use only directly reported principal metrics. Missing values remain unavailable. Liu (2023) reports an approximate lift period of 5.8 D/U; its implied reciprocal is about 0.172, but no precise St is inserted into the acceptance data. This is distinct from the explicitly reported St=0.171 values in Shim, Kamensky and Tsetoglou.
5. Compute inclusive minima/maxima from the eligible rows. No svMultiPhysics result contributes to those bounds. Each generated report names all contributing configurations, endpoint configurations and distinct-study counts.

The **passing interval is the eligible literature envelope**, with both endpoints included. A result outside that envelope fails. Passing also requires the existing run-quality checks.

| Metric | Previous saved-table range | Observed span | Passing literature range |
|---|---|---|---|
| Mean CD | 1.317–1.374 | 1.297–1.477 | 1.300–1.3917 |
| CL RMS | 0.222–0.235 | 0.2031–0.237 | 0.215–0.237 |
| St | 0.164–0.166 | 0.157–0.178 | 0.162–0.171 |

The previous saved table had already reached mean CD=1.374; the earlier approximately 1.330 upper bound described a smaller historical subset.

| Metric | Eligible studies | Eligible configurations | Lower-bound source | Upper-bound source |
|---|---:|---:|---|---|
| Mean CD | 11 | 37 | Tsetoglou VOS-IVP D/dx=50, 1.300 | Behr M240 velocity-pressure-stress, 1.3917 |
| CL RMS | 2 | 4 | Tsetoglou VOS-IVP D/dx=50, 0.215 | Tsetoglou body-fitted, 0.237 |
| St | 13 | 46 | Kravchenko TF-73 case 4, 0.162 | Shim fixed-domain formulations, Kamensky M1 tau=1000, Tsetoglou VOS50: 0.171 |

The observed drag maximum, Kamensky l=0 at 1.477, does not set acceptance because the authors explicitly demonstrate quadrature error. The observed St extremes include Alam's coarse-time-step 0.157 and Tsetoglou's coarsest-grid 0.178. The lift observed minimum, Qu D11=0.2031, reflects time-step damping. These are evaluated as configurations, not omitted merely because they are numerical extremes.

## Study-specific audit

| Original source | Recorded configurations | Comparable metrics and eligibility |
|---|---:|---|
| [Qu et al. (2013)](https://www.tfd.chalmers.se/~lada/postscript_files/lixia_paper_CPL.pdf), Tables 2–3, Sections 3.1/4.1 | 11 | All tabulated CD/St tests retained: source reports small changes. RMS acceptance uses the time-resolved reference D9 and selected production step D10. Larger steps damp lift RMS. D10 is not counted twice. |
| [Park et al. (1998)](https://doi.org/10.1007/BF02942594), Table 1, p.1203 | 1 | Mean CD=1.33 verified; lift .3321 is amplitude. Precise St/RMS secondary claims are separately retained, not treated as explicit primary numbers. |
| [Stalberg et al. (2006), author version](https://www.mech.kth.se/~eriks/stalberg_050420.pdf), Paper 2 Table 2, printed p.85 | 1 | Mean CD=1.32 verified; lift .330 is amplitude. Final-version setup discrepancies and secondary RMS/St remain unresolved. |
| [Kravchenko and Moin (1998), TF-73](https://web.stanford.edu/group/tfsa/TF_reports/TF-073_Kravchenko.pdf), Table 5.1.1, pp.64–66 | 9 | St verified; frequency/domain independence stated for Rd>=60, so cases 4–9 contribute St. Earlier-domain cases observed only. Table lift and drag maxima are not assigned as mean drag/RMS. Rd is normalized by cylinder radius: total domain diameter/D=Rd. |
| [Meneghini et al. (2001)](https://doi.org/10.1006/jfls.2000.0343), Section 3.1/Table 1 | 1 | Isolated-cylinder validation only: mean CD=1.37 and St=.165. No interacting-cylinder coefficients. |
| [Alam et al. (2017)](https://eng.monash.edu/lbe/Publications/2017/AlZhHo-JFS-2017.pdf), Table 1/Section 2.2 | 20 | Nine configurations at Nn>=69,000 and dt<=.01 contribute CD/St: the author's selected adequate setup and its mesh/time refinements. Remaining 11 sensitivity points observed only. Lift is maximum, not RMS. |
| [Kang and Masud (2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10470433/), Section 4.4/Figure 12/Table 2 | 1 | Mean CD=1.374, St=.168; lift extrema +/- .332. Domain x=[-12,36], y=[-12,12] in D units. |
| [Behr et al. (1995)](https://www.tafsm.org/PUB_PRE/jALL/j56-CMAME-LateBoun.pdf), Tables 1–4 | 10 | Five domains and two formulations kept separately. M240/M320 enter acceptance for both formulations, using the domain-insensitive plateau beyond A=16 discussed by the authors. M090/M125/M160 remain observed confinement tests. Lift is amplitude. |
| [Deng, Mao and Xie (2019)](https://doi.org/10.1103/PhysRevE.100.053107), p.053107-4 | 1 | Bare-cylinder validation: CD=1.374, St=.164, peak-to-peak lift .668. Selected H=.0125D, dt=.0005. Filament-attached cases excluded. |
| [Mittal and Tezduyar (1992)](https://home.iitk.ac.in/~smittal/publi_%26_present/sm_journals/notes_on_the_space_time.pdf), Section 5.3.1 | 1 | Fixed cylinder: directly reported St=.167. No invented precise drag/RMS from force-history plots; later forced/VIV cases excluded. |
| [Shim et al. (2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8299810/), Section 3.2 | 3 | CFD and explicitly fixed-domain FSI/BFSI reductions each report St=.171; these are not moving cylinders. Missing exact drag/RMS remain blank. Geometry/dt/BC details depend on the unrecovered earlier supplement. |
| [Kamensky et al. (2015)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4274080/), Section 3.4/Tables 1,2,4 | 6 | Four l=3 immersed configurations and one shared boundary-fitted reference contribute CD/St. The l=0 under-integration test is observed only. M1 tau=1000 explicitly gives St=.171. Lift extrema are not RMS. |
| [Tumkur et al. (2013)](https://doi.org/10.1016/j.jfluidstructs.2013.03.008), Table 1, p.221 | 1 | Stationary validation: CD=1.36, St=.168, maximum lift=.328. Excludes sprung-cylinder/nonlinear-attachment simulations. Full six-author metadata retained. |
| [Singh and Mittal (2005), drag-crisis paper](https://home.iitk.ac.in/~smittal/publi_%26_present/sm_journals/drag_crisis.pdf), Section IV.D | 1 | Re=100 St=.1643 explicitly stated, lift amplitude=.319. Precise drag not assigned from a graph. Same-year VIV publication is a distinct unresolved candidate. |
| [Liu et al. (2023)](https://link.springer.com/article/10.1007/s00466-022-02230-6), Sections 3.1–3.4 | 1 setup-only | Short-wake cylinder domain, mesh and approximate period verified. No directly reported principal metric inserted. The downstream 350D carrier calculation is not a separate cylinder-force validation case. |
| [Tsetoglou et al. (2025)](https://doi.org/10.1002/fld.5334), Section 3.2/Tables 2–3 | 4 | Body-fitted reference and selected finest VOS-IVP case contribute CD/RMS/St; two coarsenings remain observed sensitivity tests. RMS explicitly identified in Section 3.2.2. Source is independent of the svMultiPhysics results. |

### Behr coordinate normalization

Figure 1 uses native cylinder diameter 2, inlet distance 16, outlet distance 45 and lateral half-height A. Thus upstream/downstream distances are **8D/22.5D**, lateral half-height is **A/2 D**, and total height is **A D**. M240 is y=+/-12D, not +/-24D. Native dt=.125 with U=1 becomes dt*U/D=.0625. Both formulations retain their own coefficients for each domain.

### Unresolved sources and legacy claims

The candidate CSV records full identified bibliographic metadata, the claim, its origin, and the verification problem. None enters either verified range.

- **Posdziech and Grundmann (2007); Mittal (2005):** older secondary numbers retained; original numerical tables/definitions were not recovered. Mittal's main investigation constrains wake symmetry, so any unconstrained validation case must be identified explicitly before using those claims.
- **Park, Kravchenko, Stalberg legacy claims:** previously saved secondary RMS values remain preserved, alongside any original amplitude/maxima now verified. This policy does not assert an RMS value from a sinusoidal approximation.
- **Wang, Fan and Cen (2009):** DOI 10.1016/j.jcp.2008.10.038, correct pages **1504–1520**, not the abbreviated/incorrect page range in Kang's bibliography. CD=1.379, St=.170 remain secondary claims.
- **Wang and Zhang (2011):** DOI 10.1016/j.jcp.2011.01.045; publisher metadata verified, Re=100 numerical table unavailable. CD=1.330/St=.166 unresolved.
- **Prabhakar, Pontaza and Reddy:** publisher confirms **2008**, DOI 10.1016/j.cma.2007.06.013; Kang's reference year 2012 is inconsistent. CD=1.356/St=.167 remain unverified numerical claims.
- **Rajani, Kandasamy and Majumdar (2009):** DOI 10.1016/j.apm.2008.01.017. Original repository PDF retrieval failed. Secondary CD is variously 1.340 or 1.335; St=.157 is preserved but does not set acceptance until the original configuration is verified.
- **Main and Scovazzi (2018), Part II:** DOI 10.1016/j.jcp.2018.01.023. CD=1.350/St=.167 remain secondary; Part I concerns Poisson/Stokes and must not be substituted.
- **Ateshian et al. (2018):** main article DOI 10.1115/1.4038716 checked; relevant cylinder test lies in Supplement S3, whose download was not recovered. No numbers copied from Shim's later study.
- **Liu et al. (2022):** DOI 10.1142/S0218202522500579. Author/institutional abstract confirms the Re=100 parameter study; full numerical tables not recovered. No values or artificial rows generated from its parameter count.
- **Singh and Mittal (2005), VIV paper:** DOI 10.1016/j.jfluidstructs.2005.05.011. Distinct from the verified drag-crisis paper; no moving-cylinder coefficients used as a fixed-cylinder reference.

## Configuration and reports

Eligibility is configurable by editing `acceptance_metrics` and the documented reason in the reference CSV. Its metric names are `mean_CD`, `CL_rms`, `strouhal`, separated by semicolons. Only independently verified metrics can be made eligible. `observed_eligible` and physical-context fields control the observed dataset. Configuration IDs must be unique. Unverified claims belong in the candidate CSV.

The configuration selects `acceptance.method = "literature_min_max"`. The passing bounds are calculated directly from eligible reference values; no fixed bounds or svMultiPhysics baseline numbers are embedded in the validator.

Refresh only the literature decisions from an existing report:

```bash
python main/postprocessing/validate_results.py results/<case-directory>
```

This command preserves computed statistics, selected time windows, mean lift, quality checks and other saved diagnostics, updating `summary.json` and `tables/metrics.csv`. It does not read/reintegrate VTUs or regenerate meshes/figures/force histories. The general `run_case.py --postprocess-only` command performs the fuller existing postprocessing workflow and is unnecessary for a literature-only refresh.

The CSV report includes observed bounds, eligible passing bounds, source IDs and contributing-study/configuration counts. The summary also includes exclusion reasons and policy version. Missing literature metrics are reported as unavailable individually; they cannot pass by omission. Exit codes remain 0 for PASS, 1 for outside the literature range, 2 for invalid run quality and 3 for unavailable comparison.

PASS means meeting this broad literature criterion and existing quality checks, for either a preset or user-defined mesh. It is not a comparison to an earlier svMultiPhysics baseline; no such baseline gate currently exists in this repository. A future strict regression gate should be reported separately. A changed acceptance decision changes the interpretation of the saved result, not its numerical accuracy.

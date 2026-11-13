# Literature reference values

`reference_values.csv` lists published comparison values for cylinder flow at Re=100. The benchmark compares mean drag, lift RMS, and Strouhal number against the inclusive minimum and maximum of the saved compatible Re=100 values.

The `verification` column distinguishes values checked in source publications from tabulated values and unverified candidates. Empty cells mean unconfirmed or unreported, not zero. Domain distances are measured from the cylinder center and normalized by D.

- Qu et al. (2013), Table 3 and case D9: RMS explicitly means standard deviation of fluctuations. The author manuscript is available at https://www.tfd.chalmers.se/~lada/postscript_files/lixia_paper_CPL.pdf .
- Alam et al. (2017): the single-cylinder comparison in Table 2 is `CD=1.372`, `St=0.165`, with the selected `dt=0.01` case. Table 1 also gives `CD=1.372`, `St=0.166` at `dt=0.005` on 82,000 nodes. The lift column is maximum lift, **not RMS**. The single-cylinder height of 20D follows from its stated 5% blockage; do not substitute the four-cylinder domain height of 80D.
- Meneghini et al. (2001), Section 3.1 and Table 1: single-cylinder Re=100 drag 1.37 and St 0.165; dt=0.005, 27,062 triangles and 13,696 nodes. Exact outer-domain dimensions remain unconfirmed. An author-uploaded original is at https://www.researchgate.net/publication/257280881_Numerical_simulation_of_flow_interference_between_two_circular_cylinders_in_tandem_and_side-by-side_arrangement .
- Deng, Mao and Xie (2019): the paper identity was confirmed, but the supplied bare-cylinder values, domain dimensions, and time step were not independently checked in accessible full text. This candidate row is included in the current saved-table min/max comparison; its verification status remains explicitly recorded.
- Other tabulated values retain their verification notes. Rows labelled RMS are used for lift RMS even when the original definition has not been independently audited; maximum-lift entries are excluded.

Mean lift is calculated and retained in the numerical statistics, but it is not included in this three-metric literature comparison.

Bounds are derived dynamically from the CSV, with missing values omitted. The current ranges are mean drag [1.317, 1.374], lift RMS [0.222, 0.235], and Strouhal number [0.164, 0.166]. Outside values are reported as a percentage above or below the nearest bound: `100 * abs(simulation - bound) / abs(bound)`. Source IDs accompany each comparison.

The range of published drag values is not a statistical uncertainty interval. Method, mesh, domain, boundary conditions, averaging window, and time step differ. Similar drag values do not identify the cause of differences from other studies.

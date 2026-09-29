# Replication: the Monte Carlo validation of equaids's standard errors

This folder lets anyone rebuild the data-generating processes (DGP) used to
validate the analytic standard errors of `equaids`, check that they obtain
exactly the same populations, and rerun the Monte Carlo.

## What is validated

The analytic standard errors of the aggregate elasticities and of the free
coefficients, under simple random sampling (`vce(robust)`) and under a
stratified cluster design (`vce(svy)` with the finite-population correction),
for AIDS and QUAIDS with Ray's demographic scaling. The DGP is the **perfect
case**: the model is correctly specified, consistent with consumer theory,
and every parameter is clearly identified. Weak identification is the
business of `equaidsdiag`, not of this validation.

## The DGP

- **Population**: 100,000 households, 10 strata x 400 PSUs x 25 households,
  generated once per model x error type with a fixed seed (six populations).
- **Regressors**: log prices set at the PSU level (spatial variation) plus a
  small household noise; a log-normal expenditure,
  ln x = 5.5 + 0.12 (stratum - 5.5) + PSU effect N(0, 0.1) + N(0, 0.6);
  demographics z1 (indicator, p = 0.4) and z2 (household size, 1 to 8).
- **Shares**: the model equations (equaids's own `_eq_parts`) at known
  parameters, plus errors drawn once: normal, standard deviations 0.2 x the
  mean shares with the correlations of a multinomial; heteroskedastic
  (scaled by exp((ln x - mean)/2)) or clustered (30% of the variance at the
  PSU level). The last share follows by adding-up.
- **alpha_0**: the smallest ln x of the population minus 0.1, fixed in the
  census and in every sample (the default rule of `equaids` was validated
  separately, by the bootstrap).
- **Parameters**: mean shares 40/30/20/10; gamma = 0.3 (wbar wbar' - diag(wbar));
  eta = (0.03, -0.012, -0.010, -0.008) for z1, (0.01, -0.004, -0.003, -0.003)
  for z2.
  - QUAIDS: beta = (-0.15, 0.05, 0.06, 0.04), lambda = (0.02, -0.008,
    -0.007, -0.005), rho = (1.2, 0.3);
  - AIDS: linear Engel curves, beta = (-0.06, 0.02, 0.025, 0.015),
    rho = (1.6, 0.3).
  The file `mcs_<cell>_dgp.csv` (in `audit/mcs/`) lists the full parameter vector of each cell.
- **Truth**: the census, `equaids` on the whole population.

## Checks passed before any replication

The script refuses to draw samples unless the DGP passes:

| Check | Criterion |
|---|---|
| Slutsky | s_ij = w_i e*_ij symmetric and negative semidefinite for every household |
| Shares | systematic shares in [0, 1] for every household |
| Engel curvature | turning points -beta_i/(2 lambda_i) inside the range of l; condition index of (1, l, l^2) below 100 (equaidsdiag D4) |
| Identification | every free parameter with abs(theta)/SE >= 4 at n = 5,000 (census SE scaled by sqrt(N/n)) |

The parameters were tuned in pilots (`tune*.do` in `audit/mcs/`) until every
check passed; the smallest abs(theta)/SE is 4.2 to 5.3 over the six cells.

## Designs

- `srs`: 5,000 households drawn independently (with replacement), no
  weights, `vce(robust)`.
- `svy`: one-stage stratified cluster sample, 200 PSUs drawn without
  replacement within strata, poorer strata oversampled (32 28 24 24 20 20
  16 16 10 10), all their households; weights = inverse probabilities;
  `svyset psu [pw=pw], strata(stratum) fpc(NPSU)` and `vce(svy)`; t with
  190 design degrees of freedom.

## How to rerun

```stata
* one cell: model (quaids|aids), errors (homo|hetero|clus), design (srs|svy),
* n, replications, part
do mcs_equaids.do quaids homo srs 5000 1000 1
```

The twelve cells (2 models x 3 error types x 2 designs), 1,000 replications
each, take about 45 minutes on four Stata processes. `tools/stata_queue.sh`
runs a list of do-files (one wrapper per cell, which sets the paths and calls
`mcs_equaids.do`) in at most K parallel processes, each taking the next cell
as soon as it is free:

```bash
bash tools/stata_queue.sh 4 cell1.do cell2.do ...
```

The seed of the samples depends only on the cell and the part, not on the
process: a rerun draws exactly the same samples, so that a change of the
variance formula is measured on the same estimates. The two designs of a
given model and error type share the same population. `mcs_summary.py` writes
`mcs_summary.md` (coverage of the 95% intervals, bias, standard error over
standard deviation; t with the design degrees of freedom for `vce(svy)`).

## Results

Coverage of the 95% intervals and mean analytic SE / Monte Carlo SD of the
estimates (1,000 replications; Monte Carlo bands: coverage 0.933 to 0.967,
ratio 0.944 to 1.056), aggregate elasticities and free coefficients:

| Design | Coverage | SE / SD |
|---|---|---|
| `srs`, `vce(robust)` | 0.945 to 0.954 | 0.98 to 1.02 |
| `svy`, `vce(svy)` | 0.932 to 0.948 | 0.97 to 1.01 |

All replications converged (9 to 11 iterations at the median, at most 26),
and no estimate lies more than 10 robust standard deviations from the median.
Details by cell and by family of statistics: `audit/mcs/mcs_summary.md`.

*Information (not part of the replication).* With 100 PSUs (5 to 16 per
stratum, PSUs of 50 households) the design coverage was 0.925 to 0.944 and
the ratio 0.94 to 0.99; doubling the number of PSUs brought every cell
toward 0.95 while the `srs` cells, on the same populations, did not move.
The linearized design variance needs enough PSUs per stratum -- a
small-sample limit of the linearization, not of the formula.

## Checking the populations

Each cell records the Stata `datasignature` of its population. Running the
script with 0 replications rebuilds the population, runs the checks and
prints the signature, to compare with the one stored in
`mcs_<cell>_truth.dta`:

```stata
do mcs_equaids.do quaids homo srs 5000 0 1
```

Signatures (`datasignature`) of the six populations:

| Model | Errors | Signature |
|---|---|---|
| QUAIDS | homo | 100000:15(125485):3538290244:1742659787 |
| QUAIDS | hetero | 100000:15(125485):4114593653:299916656 |
| QUAIDS | clus | 100000:15(125485):716497240:1513381787 |
| AIDS | homo | 100000:15(125485):3884304284:3257601964 |
| AIDS | hetero | 100000:15(125485):4110030658:458639672 |
| AIDS | clus | 100000:15(125485):1635707945:1177785054 |

The populations themselves (about 10 MB each) are archived on Zenodo with the
technical note, for those who prefer not to regenerate them.

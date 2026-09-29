# equaids — replication

Scripts that reproduce the numbers of the technical note. Run them from this
folder (`cd replication`); each writes its results into `out/` (or into the
`out/` of its subfolder). `master.do` runs the fast ones in one go.

Data: Poi's data come from Stata (`webuse food`); the Mexican cereals are in
`../examples/`; the data of Lecocq and Robin (2015) are downloaded once from
the Stata Journal by `get_lr_data.do` (package st0393_3) into `data/`.

## Scripts by section of the note

| script | in the note | time |
|---|---|---|
| `types_elasticities.do` | §4.1: the five types of elasticities on the Mexican cereals (table after the remark) | seconds |
| `remark_alpha0.do` | §4.1, remark: corn and rice, market vs mean of the household elasticities, α0 = 6.31 and 9.5, shares below 0.001, design SE of corn | seconds |
| `pimpute_effect.do` | §6.1 and §9.3: weight of the imputation of the prices in the robust SE (Mexico, with and without the correction; the design of §6.4) | 1 min |
| `selection.do` | §6.6: Mexican cereals without and with the correction for the non-buyers (`vce(bootstrap)`, 100 replications), `e(m0_t)`, linearized vs bootstrap SE of δ | minutes |
| `dgp_nonbuyers.do` | §6.4: generator of the population of the Monte Carlo of the correction | — |
| `dgp_describe.do` | §6.4: buyers and pseudo-R² of that population | seconds |
| `make_example_data.do` | §8.1: the simulated data of example 8 of `help equaids`, installed with the package (`../examples/equaids_nonbuyers.dta`, 3,000 households of that population) | seconds |
| `mc_nonbuyers.do` | §6.4: Monte Carlo of the correction (700 samples); raw draws shipped in `out/mc_nonbuyers_raw.dta`, `do mc_nonbuyers.do combine 4` rebuilds the summary | hours (combine: seconds) |
| `diagnostics.do` | §7.2 and §8: `equaidsdiag` (condition indexes 23, 36, 31 and 331; the Mexican report) | seconds |
| `stability.do` | §7.3, Table `tab:stab` and the paragraph after it | 4 min |
| `example_design.do` | §8, Table `tab:mex` and the estimation facts; design/robust ratios of §9 | seconds |
| `example_engel.do` | §8, Figure `fig:engel` (`out/engel_mexico.png`), turning points, corn below zero | seconds |
| `bootstrap/` | §9.4, Table `tab:boot`, the Gauss–Newton ranges and Poi's α0 = 10 (see below) | summary: 1 min |
| `bruteforce/` | §9.3, Table `tab:numif` and the decomposition (see below) | summary: seconds |
| `montecarlo/` | §9.5, Table `tab:mc` (see its README) | summary: seconds |
| `timing.do` | §9.6, Table `tab:time` (needs Poi's `quaids`, st0268_1, and Stata 18+ for `demandsys`) | more than 1 hour; **not re-run** for this version |

### `bootstrap/`
`se_boot.do` draws the replications (shipped in `raw/`); `boot_summary.do`
recomputes the analytic standard errors of the current `equaids` and of
`equaids` 1.1.0 (Gauss–Newton, downloaded from the GitHub release v1.1.0) and
scores both against the same replications.

### `bruteforce/`
The influence function by brute force: `bf_worker.do` (one block of
households), `bf_parallel.sh` (K blocks in parallel), `bf_combine.do`; the
specifications `spec_poi.do`, `spec_sel.do`, `spec_selimp.do`; the blocks U
are shipped in `raw/`. `gn_equaids110.do` computes the standard errors of
`equaids` 1.1.0 (the Gauss–Newton column), `numif_summary.do` writes
Table `tab:numif`, `decomposition.do` the paragraph that follows it (the
pieces of the influence function added one at a time). `demo_mean.do`:
the brute force on a mean, every step printed.

## Numbers from the test suite

The numerical checks quoted in the note (analytic Jacobian vs numerical,
distance to the optimum, equality with Poi and `demandsys`, aggregation
identities, the design oracle against `svy: total`, the oracles of the
correction, the 293 checks) are printed by the test suite of the development
repository (`tests/test_step*.do`), which is not distributed with the
package.

## Not re-run for this version

- `timing.do` (Table `tab:time`): the times depend on the computer; the
  table reports one run.
- The failure of Poi's `quaids` on the data of Lecocq and Robin with
  α0 = 10 (§7.1: twelve minutes, log likelihood 155,890.8 against 156,634.1).
- The Monte Carlo cells with 100 sampling units (§9.5), computed with the
  variance of version 1.1.0.

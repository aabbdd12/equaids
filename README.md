# equaids

AIDS and QUAIDS demand systems for household survey data, in Stata.

`equaids` estimates the almost ideal demand system (Deaton and Muellbauer
1980) and its quadratic extension (Banks, Blundell and Lewbel 1997), with
demographic variables entering by Ray's (1983) scaling, by iterated feasible
generalized nonlinear least squares. It reports the elasticities of the
**households** (the mean of the household elasticities, the default), of the
**individuals** (each household weighted by its size), of the **market** (the
aggregate elasticities), of a **reference** household (at the means) or the
plain mean of the household elasticities, all with analytic standard errors
built from the exact influence function of the whole procedure: robust, by
cluster, or by survey design (`svyset`:
strata, primary sampling units, finite-population correction), or by a
bootstrap of the whole procedure (`vce(bootstrap)`). For survey data with zero
shares and unit values, `pimpute()` fills the missing prices of the
non-buyers from the households of the same group, and `selection` corrects
the system for the non-buyers (Shonkwiler and Yen 1999), the standard errors
including the estimation of the probits.

- `equaids`: estimation, elasticities, standard errors, dialog box
- `equaidsdiag`: diagnostics of a specification before estimating it
- `estat diagnostics`, `estat engel`: after estimation
- every example of `help equaids` runs from its links, in the command window,
  in the dialog box or as a do-file, without losing the data in memory
  (`equaids_examples`); the Mexican example data are installed with the
  package (`sysuse mexico_2014_cereals`)

Version 1.2.0. Requires Stata 14.2 or later.

## Installation

```stata
net install equaids, from("https://raw.githubusercontent.com/aabbdd12/equaids/main") replace
help equaids
help equaidsdiag
```

## Replication

The folder `replication/` reproduces the numbers of the technical note (its
README maps each table and figure to a script); it is not installed by
`net install`. Run its scripts from that folder, in a copy of this repository.

## Citation

Araar, A. 2026. *Estimating AIDS and QUAIDS demand systems with survey data:
the equaids Stata module*. Technical note, Zenodo.
https://doi.org/10.5281/zenodo.22959991

## Author

Abdelkrim Araar, Universite Laval and Partnership for Economic Policy (PEP).

## License

GPL-3.0-or-later (see `LICENSE`).

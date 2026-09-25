# equaids

AIDS and QUAIDS demand systems for household survey data, in Stata.

`equaids` estimates the almost ideal demand system (Deaton and Muellbauer
1980) and its quadratic extension (Banks, Blundell and Lewbel 1997), with
demographic variables entering by Ray's (1983) scaling, by iterated feasible
generalized nonlinear least squares. It reports the **aggregate elasticities**
(the elasticities of total demand) with analytic standard errors built from
influence functions: robust, by cluster, or by survey design (`svyset`:
strata, primary sampling units, finite-population correction).

- `equaids`: estimation, elasticities, standard errors, dialog box
- `equaidsdiag`: diagnostics of a specification before estimating it
- `estat diagnostics`, `estat engel`: after estimation

Version 1.0.0. Requires Stata 14.2 or later.

## Installation

```stata
net install equaids, from("https://raw.githubusercontent.com/aabbdd12/equaids/main") replace
help equaids
help equaidsdiag
```

## Citation

Araar, A. 2026. *Estimating AIDS and QUAIDS demand systems with survey data:
the equaids Stata module*. Technical note, Zenodo.
https://doi.org/10.5281/zenodo.22959991

## Author

Abdelkrim Araar, Universite Laval and Partnership for Economic Policy (PEP).

## License

GPL-3.0-or-later (see `LICENSE`).

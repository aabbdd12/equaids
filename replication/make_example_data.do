* make_example_data.do -- the simulated data of example 8 of help equaids,
* an ancillary file of the package (equaids_nonbuyers.dta): 3,000 households
* drawn from the population of the Monte Carlo of the correction for the
* non-buyers (Section 6.4, dgp_nonbuyers.do: rho = 0.3, beta = (-0.10,
* 0.08), every parameter strongly identified).  Writes
* ../examples/equaids_nonbuyers.dta (Stata 14 format).
*
*   cd <repository>/replication
*   do make_example_data.do
version 14.2
clear all
set more off
local ROOT = subinstr("`c(pwd)'", "\", "/", .) + "/.."
capture confirm file "`ROOT'/src/equaids.ado"
if _rc {
    display as error "run this script from the replication/ directory"
    exit 601
}
run dgp_nonbuyers.do
_dgp_nonbuyers 3000 0.3 -0.10 0.08 20260929
keep w1 w2 w3 p1 p2 p3 x hs q1 q2 grp
order w1 w2 w3 p1 p2 p3 x hs q1 q2 grp
label variable w1  "budget share of good 1 (0 for the non-buyers)"
label variable w2  "budget share of good 2 (0 for the non-buyers)"
label variable w3  "budget share of good 3 (bought by every household)"
label variable p1  "price of good 1 (missing for the non-buyers)"
label variable p2  "price of good 2 (missing for the non-buyers)"
label variable p3  "price of good 3"
label variable x   "total expenditure"
label variable hs  "household size (Ray's scaling)"
label variable q1  "variable of the probit of good 1 only"
label variable q2  "variable of the probit of good 2 only"
label variable grp "group of 25 households sharing a price shock (sampling unit)"
label data "equaids: simulated non-buyers, 3,000 households (technical note, Section 6.4)"
compress
saveold "`ROOT'/examples/equaids_nonbuyers.dta", version(14) replace

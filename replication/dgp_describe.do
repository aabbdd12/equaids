* dgp_describe.do -- Section 6.4: the population of the Monte Carlo of the
* correction for the non-buyers (dgp_nonbuyers.do, rho = 0.3, beta =
* (-0.10, 0.08)): the share of buyers of the two corrected goods and the
* pseudo-R2 of their probits, on one large draw (100,000 households), from
* the selection table of equaids (e(sel_diag): % buyers, McFadden pseudo-R2,
* households predicted with probability 0 or 1).
*
*   cd <repository>/replication
*   do dgp_describe.do
version 14.2
clear all
set more off
local ROOT = subinstr("`c(pwd)'", "\", "/", .) + "/.."
capture confirm file "`ROOT'/src/equaids.ado"
if _rc {
    display as error "run this script from the replication/ directory"
    exit 601
}
adopath ++ "`ROOT'/src"
run dgp_nonbuyers.do
_dgp_nonbuyers 100000 0.3 -0.10 0.08 20260929
quietly equaids w1 w2 w3, prices(p1 p2 p3) expenditure(x) noquadratic demographics(hs) anot(0) ///
    pimpute(grp) selection selvars(w1: q1 ; w2: q2) noelastse notable nolog
tempname D
matrix `D' = e(sel_diag)
display as text _n "100,000 households: % buyers, pseudo-R2, predicted 0 or 1"
matrix list `D', format(%9.3f)

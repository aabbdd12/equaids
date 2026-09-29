* selection.do -- the non-buyers in the note (equaids 1.2.0): Mexican cereals,
* QUAIDS with the demographics hhsize and isMale, sampling weights, robust
* standard errors, the prices of the non-buyers filled by pimpute(psu rururb).
* For each good: the market expenditure and own-price elasticities, with
* standard errors, without and with the correction for the non-buyers
* (selection, perc_ocupa in the probits only), and the selection table
* (buyers, pseudo-R2, VIF of beta, delta and its bootstrap and linearized
* standard errors). Under the
* correction, household size in Ray's scaling brings the estimate near the
* boundary m0 > 0 (e(m0_t) below 3): its standard errors are those of
* vce(bootstrap) (100 replications of the whole procedure, seed 2026); the
* distance e(m0_t) of both estimations is displayed.
* Writes out/selection.csv and out/selection_diag.csv. Run from replication/:
*
*   cd <repository>/replication
*   do selection.do
version 14.2
clear all
set more off
local ROOT = subinstr("`c(pwd)'", "\", "/", .) + "/.."
adopath ++ "`ROOT'/src"
capture mkdir out
use "`ROOT'/examples/mexico_2014_cereals.dta", clear
local W "wcorn wwheat wrice wother wcomp"
local P "pcorn pwheat price pother pcomp"
local O "prices(`P') expenditure(hh_current_inc) demographics(hhsize isMale) pimpute(psu rururb) nolog notable"

tempname T D
matrix `T' = J(5, 8, .)
local c 0
foreach s in "" "selection selvars(perc_ocupa) vce(bootstrap, reps(100) seed(2026))" {
    quietly equaids `W' [pw=sweight], `O' `s'
    display as text "`s': N = " e(N) ", converged " e(converged) ", vce " e(vce) ", m0_t " %6.2f e(m0_t)
    tempname X VX U VU
    matrix `X'  = e(elas_x)
    matrix `VX' = e(V_elas_x)
    matrix `U'  = e(elas_u)
    matrix `VU' = e(V_elas_u)
    forvalues j = 1/5 {
        matrix `T'[`j', `c' + 1] = `X'[1, `j']
        matrix `T'[`j', `c' + 2] = sqrt(`VX'[`j', `j'])
        matrix `T'[`j', `c' + 3] = `U'[`j', `j']
        local d = (`j' - 1) * 5 + `j'
        matrix `T'[`j', `c' + 4] = sqrt(`VU'[`d', `d'])
    }
    local c = `c' + 4
}
matrix `D' = e(sel_diag), e(sel_delta)', e(se_sel_delta)'
* the linearized (robust) standard errors of delta, against the bootstrap
quietly equaids `W' [pw=sweight], `O' selection selvars(perc_ocupa)
matrix `D' = `D', e(se_sel_delta)'
mata: st_local("rq", invtokens(strofreal(st_matrix("`D'")[., 6] :/ st_matrix("`D'")[., 7], "%5.2f")'))
display as text "delta: bootstrap standard error / linearized one, by good: " as result "`rq'"
matrix colnames `T' = ex se_ex eu se_eu ex_sel se_ex_sel eu_sel se_eu_sel
matrix rownames `T' = corn wheat rice other composite
matrix colnames `D' = buy_pct pseudo_r2 perfect vif_beta delta se_delta se_delta_lin
matrix rownames `D' = corn wheat rice other
matrix list `T', format(%9.4f)
matrix list `D', format(%9.4f)

preserve
clear
svmat double `T', names(col)
gen str10 good = ""
local g 0
foreach n in corn wheat rice other composite {
    local ++g
    replace good = "`n'" in `g'
}
order good
export delimited using "out/selection.csv", replace
restore
preserve
clear
svmat double `D', names(col)
gen str10 good = ""
local g 0
foreach n in corn wheat rice other {
    local ++g
    replace good = "`n'" in `g'
}
order good
export delimited using "out/selection_diag.csv", replace
restore

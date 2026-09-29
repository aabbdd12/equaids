* mc_nonbuyers.do -- the Monte Carlo of the correction for the non-buyers in
* the note (equaids 1.2.0): 700 samples of N = 3,000 households from the
* perfect case of dgp_nonbuyers.do (rho = 0.3, beta = (-0.10, 0.08)), each
* estimated by
*
*   equaids w1 w2 w3, prices(p1 p2 p3) expenditure(x) noquadratic
*       demographics(hs) anot(0) pimpute(grp) selection
*       selvars(w1: q1 ; w2: q2) vce(cluster grp)
*
* (the prices share a shock by group: the population is clustered by grp, and
* so is the variance). For rho, delta and the market elasticities (expenditure
* and own price): the mean analytic standard error against the standard
* deviation of the estimates, with the 95% band of that ratio allowing for
* the kurtosis of the estimates; the bias of rho and delta against the truth.
* Seeds 100001-100200 and 700001-700500, one per sample: the samples do not
* depend on how the replications are split.  Writes the estimates and
* standard errors of every sample to out/mc_nonbuyers_raw.dta, then the
* summary to out/mc_nonbuyers.csv.  Run from replication/ (about 5 minutes):
*
*   cd <repository>/replication
*   do mc_nonbuyers.do
*
* or in parallel Stata processes, each a block of samples, then combined:
*
*   do mc_nonbuyers.do 1 175 1        (samples 1-175, out/mc_nonbuyers_raw_1.dta)
*   ...                                (other blocks, other processes)
*   do mc_nonbuyers.do combine 4      (appends the 4 blocks, then the summary)
args first last part
version 14.2
clear all
set more off
local ROOT = subinstr("`c(pwd)'", "\", "/", .) + "/.."
adopath ++ "`ROOT'/src"
capture mkdir out
run dgp_nonbuyers.do

local R 700
local names r d1 d2 x1 x2 x3 u1 u2 u3
local senames ""
foreach v of local names {
    local senames `senames' se_`v'
}
if "`first'" == "combine" {
    * the blocks of the parallel run, appended
    clear
    forvalues k = 1/`last' {
        append using "out/mc_nonbuyers_raw_`k'.dta"
    }
    sort rep
    save "out/mc_nonbuyers_raw.dta", replace
    quietly count
    display as text "samples in the blocks: " r(N) " of `R'"
}
else {
if "`first'" == "" {
    local first 1
    local last `R'
}
local raw = cond("`part'" == "", "out/mc_nonbuyers_raw.dta", "out/mc_nonbuyers_raw_`part'.dta")
tempname res
postfile `res' int rep `names' `senames' using "`raw'", replace
local nf 0
forvalues r = `first'/`last' {
    local seed = cond(`r' <= 200, 100000 + `r', 700000 + `r' - 200)
    _dgp_nonbuyers 3000 0.3 -0.10 0.08 `seed'
    capture quietly equaids w1 w2 w3, prices(p1 p2 p3) expenditure(x) noquadratic demographics(hs) ///
        anot(0) pimpute(grp) selection selvars(w1: q1 ; w2: q2) vce(cluster grp) notable nolog
    if _rc | !e(converged) {
        local ++nf
        continue
    }
    tempname D SD X VX U VU
    matrix `D' = e(sel_delta)
    matrix `SD' = e(se_sel_delta)
    matrix `X' = e(elas_x)
    matrix `VX' = e(V_elas_x)
    matrix `U' = e(elas_u)
    matrix `VU' = e(V_elas_u)
    post `res' (`r') (_b[rho:rho_hs]) (`D'[1,1]) (`D'[1,2]) (`X'[1,1]) (`X'[1,2]) (`X'[1,3]) ///
        (`U'[1,1]) (`U'[2,2]) (`U'[3,3]) ///
        (_se[rho:rho_hs]) (`SD'[1,1]) (`SD'[1,2]) (sqrt(`VX'[1,1])) (sqrt(`VX'[2,2])) (sqrt(`VX'[3,3])) ///
        (sqrt(`VU'[1,1])) (sqrt(`VU'[5,5])) (sqrt(`VU'[9,9]))
    if mod(`r', 100) == 0 display as text "  `r' samples"
}
postclose `res'
display as text "failed samples: `nf'"
* a block of the parallel run stops here (the summary needs all the blocks)
if "`part'" != "" exit
}
use "out/mc_nonbuyers_raw.dta", clear
tempname T
matrix `T' = J(9, 7, .)
local k 0
foreach v of local names {
    local ++k
    quietly summarize `v', detail
    local sd = r(sd)
    matrix `T'[`k', 1] = r(mean)
    matrix `T'[`k', 2] = `sd'
    * the kurtosis of the estimates widens the sampling error of an sd
    matrix `T'[`k', 5] = r(kurtosis)
    matrix `T'[`k', 6] = 1.96 * sqrt((r(kurtosis) - 1) / (4 * r(N)))
    quietly summarize se_`v'
    matrix `T'[`k', 3] = r(mean)
    matrix `T'[`k', 4] = r(mean) / `sd'
}
* bias of rho and delta, in standard errors of the Monte Carlo mean
local k 0
foreach p in "r 0.3" "d1 -0.012" "d2 0.008" {
    local ++k
    gettoken v t : p
    quietly summarize `v'
    matrix `T'[`k', 7] = (r(mean) - `t') / (r(sd) / sqrt(r(N)))
}
matrix colnames `T' = mean sd mean_se ratio kurtosis band95 bias_z
matrix rownames `T' = `names'
matrix list `T', format(%9.4f)
clear
svmat double `T', names(col)
gen str4 param = ""
local k 0
foreach v of local names {
    local ++k
    replace param = "`v'" in `k'
}
order param
export delimited using "out/mc_nonbuyers.csv", replace

* boot_summary.do -- Table tab:boot of the note (analytic standard errors /
* bootstrap standard deviations) and the numbers of the paragraph on Poi's
* alpha_0 = 10.  The replications are shipped in raw/ (se_boot.do recomputes
* them); the analytic standard errors are recomputed here with the current
* equaids (one estimation per case, seconds), so that an improvement of the
* variance formula is scored against the same replications.  Writes
* out/boot_summary.csv.
*
*   cd <repository>/replication/bootstrap
*   do boot_summary.do
*
* Ratio = analytic standard error / standard deviation of the converged
* replications, for the market expenditure, uncompensated and compensated
* elasticities and for the free coefficients theta.  Monte Carlo band of a
* ratio: 1 +/- 2.5/sqrt(2B).  The same ratios for theta with the
* Gauss-Newton variance of equaids 1.1.0 (the note's comparison; for Poi's
* alpha_0 = 10 it is also the variance of Poi's quaids and of demandsys),
* from the ado-files of the release v1.1.0 in ../bruteforce/v110/
* (downloaded from GitHub if absent).
version 14.2
clear all
set more off
* the analytic side, with the current equaids (the full-sample estimation of
* each case)
global SE_ANALYTIC_ONLY 1
quietly do se_boot.do 500 food 1 auto
quietly do se_boot.do 125 mex 1 auto
quietly do se_boot.do 125 mex 1 auto quaids rw
quietly do se_boot.do 500 food 1
* the same with equaids 1.1.0 (Gauss-Newton), into gn/
local V110 = subinstr("`c(pwd)'", "\", "/", .) + "/../bruteforce/v110"
capture mkdir "`V110'"
foreach f in equaids.ado equaids_estat.ado _equaids_tabstars.ado _equaids_pimpute.ado {
    capture confirm file "`V110'/`f'"
    if _rc copy "https://raw.githubusercontent.com/aabbdd12/equaids/v1.1.0/src/`f'" "`V110'/`f'"
}
global SE_SRC "`V110'"
global SE_OUT "gn"
quietly do se_boot.do 500 food 1 auto
quietly do se_boot.do 125 mex 1 auto
quietly do se_boot.do 125 mex 1 auto quaids rw
quietly do se_boot.do 500 food 1
capture adopath - "`V110'"
global SE_SRC
global SE_OUT
global SE_ANALYTIC_ONLY

capture mkdir out
mata:
real scalar _med(real rowvector x)
{
    real colvector v
    real scalar m
    v = sort(x', 1) ; m = rows(v)
    return(mod(m, 2) ? v[(m + 1) / 2] : (v[m / 2] + v[m / 2 + 1]) / 2)
}
// one case: ratios of the elasticities and of theta
void _case(string scalar tag, real scalar M, real matrix OUT)
{
    real matrix    D, E, T
    real rowvector e0, s0, t0, st0, re, rt, bias, sg, rg
    real scalar    K, P, B, k
    string rowvector v
    stata(`"quietly use "raw/boot_eq_"' + tag + `"_analytic.dta", clear"')
    v = st_varname(1..st_nvar())
    K = sum(regexm(v, "^e[0-9]+$")) ; P = sum(regexm(v, "^t[0-9]+$"))
    e0 = st_data(1, "e" :+ strofreal(1..K)) ; s0 = st_data(1, "s" :+ strofreal(1..K))
    t0 = st_data(1, "t" :+ strofreal(1..P)) ; st0 = st_data(1, "st" :+ strofreal(1..P))
    // the Gauss-Newton standard errors of theta (equaids 1.1.0)
    sg = J(1, P, .)
    if (fileexists("gn/boot_eq_" + tag + "_analytic.dta")) {
        stata(`"quietly use "gn/boot_eq_"' + tag + `"_analytic.dta", clear"')
        sg = st_data(1, "st" :+ strofreal(1..P))
    }
    stata("clear")
    for (k = 1; k <= 9; k++) {
        if (fileexists("raw/boot_eq_" + tag + "_" + strofreal(k) + ".dta"))
            stata(`"quietly append using "raw/boot_eq_"' + tag + "_" + strofreal(k) + `".dta""')
    }
    stata("quietly keep if conv == 1")
    B = st_nobs()
    E = st_data(., "e" :+ strofreal(1..K)) ; T = st_data(., "t" :+ strofreal(1..P))
    re = s0 :/ sqrt(diagonal(variance(E)))'
    rt = st0 :/ sqrt(diagonal(variance(T)))'
    bias = (colsum(T) :/ B :- t0) :/ st0
    rg = sg :/ sqrt(diagonal(variance(T)))'
    // B; elasticities min max; expenditure min max; theta median min max;
    // theta of the smallest ratio: its ratio and bias in standard errors
    k = order(rt', 1)[1]
    OUT = (B, min(re), max(re), min(re[1..M]), max(re[1..M]), _med(rt), min(rt), max(rt),
           rt[k], bias[k], k, sum(abs(re :- 1) :> 2.5 / sqrt(2 * B)) + sum(abs(rt :- 1) :> 2.5 / sqrt(2 * B)),
           _med(rg), min(rg), max(rg))
}
end


tempname fh
file open `fh' using "out/boot_summary.csv", write text replace
file write `fh' "case,draws,elas_min,elas_max,exp_min,exp_max,theta_median,theta_min,theta_max,theta_worst,bias_worst,worst_index,outside_band,gn_theta_median,gn_theta_min,gn_theta_max" _n
mata: OUT = J(0, 0, .)
display as text _n "Table tab:boot (ratio = analytic SE / bootstrap SD)"
display as text "{hline 78}"
foreach c in food_aauto:4 mex_aauto:5 mex_aauto_rw:5 food:4 {
    gettoken tag M : c, parse(":")
    local M : subinstr local M ":" ""
    mata: _case("`tag'", `M', OUT) ; st_matrix("OUT", OUT)
    display as text %-14s "`tag'" " draws " as result OUT[1,1] ///
        as text "  elasticities " as result %5.3f OUT[1,2] "--" %5.3f OUT[1,3] ///
        as text "  theta median " as result %5.3f OUT[1,6] ///
        as text ", range " as result %5.3f OUT[1,7] "--" %5.3f OUT[1,8] ///
        as text "  outside the band " as result OUT[1,12]
    display as text "   Gauss-Newton (equaids 1.1.0): theta median " as result %5.3f OUT[1,13] ///
        as text ", range " as result %5.3f OUT[1,14] "--" %5.3f OUT[1,15]
    if "`tag'" == "food" {
        display as text "   alpha_0 = 10: expenditure elasticities " as result %5.3f OUT[1,4] "--" %5.3f OUT[1,5] ///
            as text "; smallest theta ratio (t" as result OUT[1,11] as text ", lambda_1) " ///
            as result %5.3f OUT[1,9] as text ", bootstrap mean " as result %4.2f OUT[1,10] as text " SE from the estimate"
    }
    file write `fh' "`tag'"
    forvalues k = 1/15 {
        file write `fh' "," %12.6f (OUT[1, `k'])
    }
    file write `fh' _n
}
file close `fh'

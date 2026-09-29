* remark_alpha0.do -- the remark of Section 4.1 ("Weighted by the shares, not
* a plain mean"): Mexican cereals, QUAIDS, demographics hhsize and isMale,
* sampling weights, robust standard errors.  At the default alpha_0 and at
* alpha_0 = 9.5: the market expenditure elasticities of corn and rice and
* the means of the household elasticities, with standard errors, and the
* households with a predicted share below 0.001 (e(n_fsmall)); at the
* default alpha_0, the standard error with the survey design (svyset of the
* data, single-unit strata centered as in Section 8) and those of the
* households type.  Writes out/remark_alpha0.csv.
*
*   cd <repository>/replication
*   do remark_alpha0.do
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
capture mkdir out
use "`ROOT'/examples/mexico_2014_cereals.dta", clear
local S "wcorn wwheat wrice wother wcomp"
local O "prices(pcorn pwheat price pother pcomp) expenditure(hh_current_inc) demographics(hhsize isMale) nolog notable"

tempname fh
file open `fh' using "out/remark_alpha0.csv", write text replace
file write `fh' "alpha0,good,market,se_market,hhmean,se_hhmean,households,se_households,n_small" _n
foreach a in default 9.5 {
    local ao = cond("`a'" == "default", "", "anot(`a')")
    quietly equaids `S' [pw=sweight], `O' `ao'
    tempname X VX XH VH XW VW NS
    matrix `X'  = e(elas_x)
    matrix `VX' = e(V_elas_x)
    matrix `XH' = e(elas_xh)
    matrix `VH' = e(V_elas_xh)
    matrix `XW' = e(elas_xw)
    matrix `VW' = e(V_elas_xw)
    matrix `NS' = e(n_fsmall)
    display as text _n "alpha_0 = " as result %6.2f e(anot) as text " (`a'), N = " as result e(N)
    display as text %-8s "" %10s "market" %8s "(se)" %10s "hhmean" %9s "(se)" %10s "househ." %8s "(se)" %10s "share<.001"
    local j 0
    foreach g in corn wheat rice other comp {
        local ++j
        display as text %-8s "`g'" as result %10.3f `X'[1, `j'] %8.3f sqrt(`VX'[`j', `j']) ///
            %10.2f `XH'[1, `j'] %10.4g sqrt(`VH'[`j', `j']) %10.3f `XW'[1, `j'] %8.3f sqrt(`VW'[`j', `j']) ///
            %10.0f `NS'[1, `j']
        file write `fh' "`=e(anot)',`g'," %9.6f (`X'[1, `j']) "," %9.6f (sqrt(`VX'[`j', `j'])) "," ///
            %12.6f (`XH'[1, `j']) "," %12.6f (sqrt(`VH'[`j', `j'])) "," %9.6f (`XW'[1, `j']) "," ///
            %9.6f (sqrt(`VW'[`j', `j'])) "," %8.0f (`NS'[1, `j']) _n
    }
}
file close `fh'

* the market elasticity of corn with the survey design (default alpha_0)
svyset psu [pweight=sweight], strata(strata) singleunit(centered)
quietly equaids `S', `O' vce(svy)
tempname VS
matrix `VS' = e(V_elas_x)
display as text _n "corn, market expenditure elasticity, standard error with the survey design: " ///
    as result %6.3f sqrt(`VS'[1, 1])

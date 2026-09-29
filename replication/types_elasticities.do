* types_elasticities.do -- the comparison of the types of elasticities in
* Section 4 of the note (equaids 1.2.0): Mexican cereals, QUAIDS with the
* demographics hhsize and isMale, sampling weights, robust standard errors.
* For each good: the expenditure and own-price elasticities, with standard
* errors, of the households (the default: the mean of the household
* elasticities, each household by its weight), the individuals (weights times
* household size), the market (weights times expenditure), the reference
* household (at the means) and the unweighted mean of the household
* elasticities (hhmean). Writes out/types_elasticities.csv. Run from
* replication/:
*
*   cd <repository>/replication
*   do types_elasticities.do
version 14.2
clear all
set more off
local ROOT = subinstr("`c(pwd)'", "\", "/", .) + "/.."
adopath ++ "`ROOT'/src"
capture mkdir out
use "`ROOT'/examples/mexico_2014_cereals.dta", clear
local W "wcorn wwheat wrice wother wcomp"
local P "pcorn pwheat price pother pcomp"
local O "prices(`P') expenditure(hh_current_inc) demographics(hhsize isMale) nolog notable"

tempname T
matrix `T' = J(4, 20, .)
local c 0
foreach t in households individuals market reference hhmean {
    if "`t'" == "individuals" quietly equaids `W' [aw=sweight], `O' hhsize(hhsize)
    else if "`t'" == "households" quietly equaids `W' [aw=sweight], `O'
    else quietly equaids, elasticities(`t')
    local f = cond("`t'" == "market", "", cond("`t'" == "hhmean", "h", cond("`t'" == "reference", "m", "w")))
    tempname X VX U VU
    matrix `X'  = e(elas_x`f')
    matrix `VX' = e(V_elas_x`f')
    matrix `U'  = e(elas_u`f')
    matrix `VU' = e(V_elas_u`f')
    forvalues j = 1/4 {
        matrix `T'[`j', `c' + 1] = `X'[1, `j']
        matrix `T'[`j', `c' + 2] = sqrt(`VX'[`j', `j'])
        matrix `T'[`j', `c' + 3] = `U'[`j', `j']
        local d = (`j' - 1) * 5 + `j'
        matrix `T'[`j', `c' + 4] = sqrt(`VU'[`d', `d'])
    }
    local c = `c' + 4
    * back to the household estimation after individuals
    if "`t'" == "individuals" quietly equaids `W' [aw=sweight], `O'
}
matrix rownames `T' = corn wheat rice other
matrix colnames `T' = hh_exp hh_exp_se hh_own hh_own_se ind_exp ind_exp_se ind_own ind_own_se ///
    mkt_exp mkt_exp_se mkt_own mkt_own_se ref_exp ref_exp_se ref_own ref_own_se ///
    hhm_exp hhm_exp_se hhm_own hhm_own_se
matlist `T', format(%8.3f) title("Expenditure and own-price elasticities: households, individuals, market, reference, hhmean")
* CSV
tempname fh
file open `fh' using "out/types_elasticities.csv", write text replace
local cn : colnames `T'
file write `fh' "good"
foreach n of local cn {
    file write `fh' ",`n'"
}
file write `fh' _n
local rn : rownames `T'
forvalues i = 1/4 {
    file write `fh' "`: word `i' of `rn''"
    forvalues j = 1/20 {
        file write `fh' "," %21.12g (`T'[`i', `j'])
    }
    file write `fh' _n
}
file close `fh'

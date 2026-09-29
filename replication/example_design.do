* example_design.do -- the example of the note with the survey design
* (Section "Example: Mexican cereals with the survey design", Table tab:mex):
* Mexican cereals, QUAIDS with the demographics hhsize and isMale, the default
* alpha_0; the aggregate (market) elasticities with the design standard
* errors (vce(svy): strata, sampling units, weights, singleunit(centered)) and
* the robust ones (vce(robust), same weights); the ratios design / robust of
* the price elasticities quoted in the validation section.  Writes
* out/example_design.csv.  Run from replication/:
*
*   cd <repository>/replication
*   do example_design.do
version 14.2
clear all
set more off
local ROOT = subinstr("`c(pwd)'", "\", "/", .) + "/.."
adopath ++ "`ROOT'/src"
capture mkdir out
use "`ROOT'/examples/mexico_2014_cereals.dta", clear
svyset psu [pweight=sweight], strata(strata) singleunit(centered)
local W "wcorn wwheat wrice wother wcomp"
local P "pcorn pwheat price pother pcomp"
local O "prices(`P') expenditure(hh_current_inc) demographics(hhsize isMale) nolog notable"
tempname S X VXd VUd VXr VUr U
* warm-up (the first call compiles the Mata of equaids), then the timed run
capture quietly equaids `W' in 1/2000, `O' iterate(1)
timer clear 1
timer on 1
quietly equaids `W', `O' vce(svy)
timer off 1
quietly timer list 1
matrix `S' = e(aggshare)
matrix `X' = e(elas_x)
matrix `U' = e(elas_u)
matrix `VXd' = e(V_elas_x)
matrix `VUd' = e(V_elas_u)
display as text "vce(svy): N = " e(N) ", alpha_0 = " %5.2f e(anot) ", strata " e(N_strata) ///
    ", sampling units " e(N_psu) ", single-unit strata " e(N_single) ", design df " e(df_r) ///
    ", iterations " e(iter) ", " %4.1f r(t1) " s"
quietly equaids `W' [pw=sweight], `O' vce(robust)
matrix `VXr' = e(V_elas_x)
matrix `VUr' = e(V_elas_u)
mata:
S = st_matrix("`S'") ; X = st_matrix("`X'") ; U = st_matrix("`U'")
xd = sqrt(diagonal(st_matrix("`VXd'")))' ; xr = sqrt(diagonal(st_matrix("`VXr'")))'
ud = sqrt(diagonal(st_matrix("`VUd'")))' ; ur = sqrt(diagonal(st_matrix("`VUr'")))'
M  = cols(X) ; own = (0..M-1) :* M :+ (1..M)
T  = (100 :* S)', X', xd', xr', diagonal(U), ud[own]', ur[own]'
st_matrix("T", T)
r  = ud :/ ur
printf("{txt}price elasticities (all %g): design SE / robust SE, median %6.3f, max %6.3f\n",
    cols(r), sort(r', 1)[ceil(cols(r) / 2)], max(r))
end
matrix colnames T = share ex se_ex_design se_ex_robust eu se_eu_design se_eu_robust
matrix rownames T = corn wheat rice other composite
matlist T, format(%9.3f)
preserve
clear
svmat double T, names(col)
gen str10 good = ""
local k 0
foreach g in corn wheat rice other composite {
    local ++k
    replace good = "`g'" in `k'
}
order good
export delimited using "out/example_design.csv", replace
restore

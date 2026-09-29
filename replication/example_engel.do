* example_engel.do -- the Engel curves of the example (Section 8, Figure
* fig:engel): Mexican cereals, QUAIDS with hhsize and isMale, the survey
* design (singleunit(centered)), estat engel at the means of the log prices
* and demographics.  Reports the turning points (r(turn): log expenditure and
* its percentile), the percentile of expenditure from which the predicted
* share of corn is below zero and from which its 95% band includes zero, and
* the households with a predicted share outside [0, 1] (e(n_shout)).
* Writes out/engel_mexico.png (the figure of the note) and out/engel.dta
* (the curves).
*
*   cd <repository>/replication
*   do example_engel.do
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
svyset psu [pweight=sweight], strata(strata) singleunit(centered)
quietly equaids wcorn wwheat wrice wother wcomp, prices(pcorn pwheat price pother pcomp) ///
    expenditure(hh_current_inc) demographics(hhsize isMale) snames(corn wheat rice other composite) ///
    vce(svy) notable nolog
local nout = e(n_shout)
estat engel, data("out/engel.dta", replace)
* PNG: exported to PDF, the band (a rarea) shows the seams of its segments
graph export "out/engel_mexico.png", replace width(3000)
tempname T
matrix `T' = r(turn)
display as text _n "turning points (ln x and percentile of expenditure):"
matrix list `T', format(%9.3f)
display as text "households with a predicted share outside [0, 1]: " as result `nout'
use "out/engel.dta", clear
quietly summarize pctile if _w1 < 0
display as text "the predicted share of corn is below zero from percentile " as result %5.1f r(min)
quietly summarize pctile if _lo1 < 0
display as text "its 95% band includes zero from percentile " as result %5.1f r(min)

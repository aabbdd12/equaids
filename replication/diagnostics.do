* diagnostics.do -- the numbers of equaidsdiag in the note, without any
* estimation: the condition index of (1, l, l^2) at the start (D4, Section
* 7.2) on the three data sets at the default alpha_0, and on Poi's data at
* alpha_0 = 10; the full report on the Mexican cereals (Section 8,
* "Diagnostics": households lost, zeros, small goods, correlations of the
* demographics, number of warnings).  Writes out/diagnostics.csv.
*
*   cd <repository>/replication
*   do diagnostics.do
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
tempname fh
file open `fh' using "out/diagnostics.csv", write text replace
file write `fh' "data,alpha0,N,N_lost,warnings,cond_quad,cond_max" _n

capture program drop _dg_line
program define _dg_line
    args fh lab
    display as text _n "`lab': alpha_0 " as result %6.2f r(anot) as text ", N " as result r(N) ///
        as text ", lost " as result r(N_lost) as text ", warnings " as result r(N_warn) ///
        as text ", condition index of (1, l, l^2) " as result %6.1f r(cond_quad) ///
        as text ", largest " as result %8.1f r(cond_max)
    file write `fh' "`lab'," %6.3f (r(anot)) "," (r(N)) "," (r(N_lost)) "," (r(N_warn)) "," ///
        %8.2f (r(cond_quad)) "," %10.2f (r(cond_max)) _n
end

* Poi's data: default alpha_0, then Poi's alpha_0 = 10
webuse food, clear
quietly equaidsdiag w1-w4, prices(p1-p4) expenditure(expfd)
_dg_line `fh' "Poi"
quietly equaidsdiag w1-w4, prices(p1-p4) expenditure(expfd) anot(10)
_dg_line `fh' "Poi, alpha_0 = 10"

* the Mexican cereals: the full report of Section 8
use "`ROOT'/examples/mexico_2014_cereals.dta", clear
equaidsdiag wcorn wwheat wrice wother wcomp [pw=sweight], prices(pcorn pwheat price pother pcomp) ///
    expenditure(hh_current_inc) demographics(hhsize isMale)
_dg_line `fh' "Mexico"

* Lecocq and Robin
do get_lr_data.do
quietly equaidsdiag s1-s7, prices(p1-p7) expenditure(somtot) demographics(nbpers rural)
_dg_line `fh' "Lecocq-Robin"
file close `fh'

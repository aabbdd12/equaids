* timing.do -- Table tab:time (Section 9.6): execution times of equaids, of
* Poi's quaids (Poi 2012, Stata Journal st0268_1) and of Stata's demandsys
* on the same model and alpha_0, no weights, and their log likelihoods.
* Times depend on the computer (the note: one Windows computer, Stata 19);
* the ratios and the log likelihoods do not.  Requires Stata 18 or later
* (demandsys) and Poi's quaids:
*     net install st0268_1, from(http://www.stata-journal.com/software/sj13-3)
*
*   cd <repository>/replication
*   do timing.do            the four rows (more than one hour: Poi's quaids
*                           takes 14 and 27 minutes on the data of Lecocq and
*                           Robin, demandsys 13 and 25)
*   do timing.do case       one row: food, mexico, lr_quaids or lr_aids
*
* food: Poi's data; mexico: the Mexican cereals without the households whose
* price is missing, hhsize and isMale; lr: the data of Lecocq and Robin
* (get_lr_data.do), nbpers and rural, alpha_0 = 2.3.  alpha_0 of food and
* mexico: equaids's default rule rounded to 4 decimals (demandsys refuses the
* full precision, r(481)).  Appends one line per case to out/timing.csv.
args case
version 14.2
set more off
local ROOT = subinstr("`c(pwd)'", "\", "/", .) + "/.."
capture confirm file "`ROOT'/src/equaids.ado"
if _rc {
    display as error "run this script from the replication/ directory"
    exit 601
}
capture mkdir out
if "`case'" == "" {
    capture erase "out/timing.csv"
    foreach c in food mexico lr_quaids lr_aids {
        do timing.do `c'
    }
    exit
}
clear all
adopath ++ "`ROOT'/src"
capture which quaids
if _rc {
    display as error "Poi's quaids is not installed: net install st0268_1, from(http://www.stata-journal.com/software/sj13-3)"
    exit 111
}
local nq ""
local dm "quaids"
if "`case'" == "food" {
    webuse food, clear
    local SH w1 w2 w3 w4
    local PR p1 p2 p3 p4
    local X expfd
    local Z ""
}
else if "`case'" == "mexico" {
    use "`ROOT'/examples/mexico_2014_cereals.dta", clear
    local SH wcorn wwheat wrice wother wcomp
    local PR pcorn pwheat price pother pcomp
    local X hh_current_inc
    local Z hhsize isMale
    * the estimation sample of equaids (missing prices dropped)
    foreach v of local PR {
        quietly drop if missing(`v')
    }
}
else {
    do get_lr_data.do
    local SH s1 s2 s3 s4 s5 s6 s7
    local PR p1 p2 p3 p4 p5 p6 p7
    local X somtot
    local Z nbpers rural
    if "`case'" == "lr_aids" {
        local nq "noquadratic"
        local dm "aids"
    }
}
local zopt = cond("`Z'" != "", "demographics(`Z')", "")
local zds  = cond("`Z'" != "", "demographics(`Z', scaling)", "")
if substr("`case'", 1, 3) == "lr_" local a0 2.3
else {
    tempvar lx
    gen double `lx' = ln(`X')
    quietly summarize `lx', meanonly
    local a0 = round(r(min) - 0.1, 0.0001)
}

* warm-up: the first call of equaids compiles its Mata (a fraction of a second)
capture quietly equaids `SH' in 1/500, anot(`a0') prices(`PR') expenditure(`X') `zopt' `nq' iterate(1) nolog notable
timer clear
timer on 1
quietly equaids `SH', anot(`a0') prices(`PR') expenditure(`X') `zopt' `nq' nolog notable
timer off 1
local it1 = e(iter)
local P = colsof(e(b_free))
timer on 2
quietly equaids `SH', anot(`a0') prices(`PR') expenditure(`X') `zopt' `nq' tolerance(1e-10) nolog notable
timer off 2
local it2 = e(iter)
local ll2 = e(ll)
local N = e(N)
timer on 3
quietly quaids `SH', anot(`a0') prices(`PR') expenditure(`X') `zopt' `nq' nolog
timer off 3
local llp = e(ll)
timer on 4
quietly demandsys `dm' `SH', prices(`PR') expenditures(`X') piconstant(`a0') `zds' nolog
timer off 4
local lld = e(ll)
quietly timer list
display as text _n "`case': N = " as result `N' as text ", " as result `P' ///
    as text " free parameters, alpha_0 = " as result `a0'
display as text "  equaids, default tolerance " as result %9.1f r(t1) as text " s, iterations " as result `it1'
display as text "  equaids, tolerance 1e-10   " as result %9.1f r(t2) as text " s, iterations " as result `it2' ///
    as text ", ll " as result %16.6f `ll2'
display as text "  Poi's quaids               " as result %9.1f r(t3) as text " s, ll " as result %16.6f `llp' ///
    as text " (difference " as result %9.2e (`ll2' - `llp') / abs(`ll2') as text ")"
display as text "  demandsys                  " as result %9.1f r(t4) as text " s, ll " as result %16.6f `lld' ///
    as text " (difference " as result %9.2e (`ll2' - `lld') / abs(`ll2') as text ")"
capture confirm file "out/timing.csv"
local new = _rc
tempname fh
file open `fh' using "out/timing.csv", write text append
if `new' file write `fh' "case,N,parameters,alpha0,t_equaids,t_equaids_1e10,t_poi,t_demandsys,iter,iter_1e10,ll,ll_poi,ll_demandsys" _n
file write `fh' "`case',`N',`P',`a0'," %9.2f (r(t1)) "," %9.2f (r(t2)) "," %9.2f (r(t3)) "," %9.2f (r(t4)) ///
    ",`it1',`it2'," %18.6f (`ll2') "," %18.6f (`llp') "," %18.6f (`lld') _n
file close `fh'

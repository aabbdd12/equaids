* se_boot.do (replication/bootstrap; audit/se/se_boot_equaids.do of the
* development) -- are the analytic standard errors of equaids's
* aggregate elasticities right?  Reference: the bootstrap of the households
* (drawn with replacement; with sampling weights, weight x count), the order
* of se-validation SRS, then pweight (the design comes at step 3).
*
*   data food  Poi's data, 4,048 households, QUAIDS, alpha_0 = 10, no weights
*   data mex   Mexican cereals, 10,804 households, QUAIDS, alpha_0 = 9.5,
*              demographics hhsize isMale, [pw=sweight]
*
* Statistics: the aggregate expenditure elasticities (M), uncompensated and
* compensated price elasticities (M x M each).  Each replication starts from
* the full-sample estimate (from()) and skips the analytic standard errors
* (noelastse); a replication that does not converge is recorded and dropped.
*
* Also the free coefficients theta (e(b_free)) and their analytic robust
* standard errors: the elasticity standard errors are J V_theta J' plus the
* sampling term, so a gap in the elasticities is traced to theta.
* Optional: alpha_0 (default 10 for food, 9.5 for mex) and the model (quaids,
* aids); the output files are tagged when they differ from the defaults.
* anot = auto: equaids's default rule (min ln x - 0.1) is applied to the data
* of each replication, as to the full sample: the whole procedure, alpha_0
* included, is replicated; the alpha_0 of each replication is recorded.
*
* scheme = rw (mex only): Rao-Wu bootstrap of the survey design, the
* reference for vce(svy): in each stratum n_h - 1 PSUs drawn with
* replacement, weight x count x n_h / (n_h - 1); the analytic standard errors
* are vce(svy) (svyset psu [pw=sweight], strata(stratum), no fpc).  The
* estimation sample is fixed first and the strata with a single PSU in it
* are dropped (Rao-Wu needs n_h >= 2).
*
*   cd <repository>/replication/bootstrap
*   do se_boot.do  B  data  part  [anot  model  scheme]
*
* The replications go to raw/boot_eq_<tag>_<part>.dta ($SE_OUT: another
* folder), the full-sample estimate and its analytic standard errors to
* raw/boot_eq_<tag>_analytic.dta ($SE_ANALYTIC_ONLY = 1: that part only).
args Barg data part anot model scheme
version 14.2
clear all
set more off
local B    = cond("`Barg'" == "", 200, real("`Barg'"))
local data = cond("`data'" == "", "food", "`data'")
local part = cond("`part'" == "", 1, real("`part'"))
local model = cond("`model'" == "", "quaids", "`model'")
local a0def = cond("`data'" == "food", "10", "9.5")
if "`anot'" == "" local anot `a0def'
local tag "`data'"
if "`anot'" != "`a0def'" local tag "`tag'_a`=subinstr("`anot'", ".", "p", .)'"
if "`model'" == "aids" local tag "`tag'_aids"
local rw = ("`scheme'" == "rw")
if `rw' local tag "`tag'_rw"
local nq = cond("`model'" == "aids", "noquadratic", "")
local aopt = cond("`anot'" == "auto", "", "anot(`anot')")
local ROOT = subinstr("`c(pwd)'", "\", "/", .) + "/../.."
capture confirm file "`ROOT'/src/equaids.ado"
if _rc {
    display as error "run this script from the replication/bootstrap/ directory"
    exit 601
}
adopath ++ "`ROOT'/src"
* $SE_SRC: another version of equaids, first on the ado-path
if "$SE_SRC" != "" adopath ++ "$SE_SRC"
local OUT "raw"
if "$SE_OUT" != "" local OUT "$SE_OUT"
capture mkdir "`OUT'"

if "`data'" == "food" {
    webuse food, clear
    local CMD "equaids w1-w4 [iw=wt], `aopt' prices(p1-p4) expenditure(expfd) `nq' nolog notable"
    gen double sw = 1
}
else {
    use "`ROOT'/examples/mexico_2014_cereals.dta", clear
    local CMD "equaids wcorn wwheat wrice wother wcomp [pw=wt], `aopt' `nq' prices(pcorn pwheat price pother pcomp) expenditure(hh_current_inc) demographics(hhsize isMale) nolog notable"
    gen double sw = sweight
    if `rw' {
        local CMD0 : subinstr local CMD "[pw=wt]" "[pw=sweight]"
        quietly `CMD0'
        keep if e(sample)
        capture confirm string variable strata
        if !_rc encode strata, gen(stratum)
        else gen long stratum = strata
        bysort stratum psu: gen byte _first = (_n == 1)
        bysort stratum: egen long nh = total(_first)
        drop if nh == 1
        gen long nh1 = nh - 1
        * Rao-Wu rescaling of the weights
        replace sw = sweight * nh / (nh - 1)
    }
}
* food: iweights equal to one reproduce the unweighted estimate and let the
* replications carry the bootstrap counts
gen double wt = sw
gen long cnt = .

* ---- full sample: point values and analytic standard errors ----
if `rw' {
    svyset psu [pw=sweight], strata(stratum)
    local CMDA : subinstr local CMD "[pw=wt]" ""
    quietly `CMDA' vce(svy)
}
else quietly `CMD'
matrix TH = e(b_free)
local M = e(ngoods)
matrix E0 = e(elas_x), vec(e(elas_u)')', vec(e(elas_c)')'
mata: st_matrix("SE0", (sqrt(diagonal(st_matrix("e(V_elas_x)")))', ///
    sqrt(diagonal(st_matrix("e(V_elas_u)")))', sqrt(diagonal(st_matrix("e(V_elas_c)")))'))
local K = colsof(E0)
local P = colsof(TH)
matrix SEt = vecdiag(e(V_free))
mata: st_matrix("SEt", sqrt(st_matrix("SEt")))
if `part' == 1 {
    preserve
    clear
    svmat double E0, names(e)
    svmat double SE0, names(s)
    svmat double TH, names(t)
    svmat double SEt, names(st)
    gen double a0 = e(anot)
    gen str8 data = "`data'"
    gen long N = e(N)
    save "`OUT'/boot_eq_`tag'_analytic.dta", replace
    restore
}
* $SE_ANALYTIC_ONLY = 1: the analytic standard errors only.  A change of the
* variance formula leaves the point estimates, hence the replications,
* unchanged: they are kept, and compared with the new standard errors.
if "$SE_ANALYTIC_ONLY" == "1" exit

* ---- replications ----
tempname pf
local vl ""
forvalues k = 1/`K' {
    local vl `vl' double e`k'
}
forvalues k = 1/`P' {
    local vl `vl' double t`k'
}
postfile `pf' int rep byte conv int iter double a0 `vl' using "`OUT'/boot_eq_`tag'_`part'.dta", replace
set seed `=20260924 + 1000 * `part' + ("`data'" == "mex") * 7'
timer clear 1
timer on 1
forvalues r = 1/`B' {
    preserve
    if `rw' bsample nh1, strata(stratum) cluster(psu) weight(cnt)
    else    bsample, weight(cnt)
    quietly keep if cnt > 0
    quietly replace wt = sw * cnt
    capture quietly `CMD' from(TH) noelastse
    if _rc {
        local line "(`r') (0) (0) (.)"
        forvalues k = 1/`=`K' + `P'' {
            local line "`line' (.)"
        }
    }
    else {
        matrix E = e(elas_x), vec(e(elas_u)')', vec(e(elas_c)')'
        local line "(`r') (`=e(converged)') (`=e(iter)') (`=e(anot)')"
        forvalues k = 1/`K' {
            local line "`line' (`=E[1, `k']')"
        }
        matrix T = e(b_free)
        forvalues k = 1/`P' {
            local line "`line' (`=T[1, `k']')"
        }
    }
    post `pf' `line'
    restore
    if mod(`r', 25) == 0 {
        timer off 1
        quietly timer list 1
        display as text "part `part': " as result `r' as text " of `B' replications, " ///
            as result %6.0f r(t1) as text " s"
        timer on 1
    }
}
postclose `pf'

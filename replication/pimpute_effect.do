* pimpute_effect.do -- Section 6.1 (and the paragraph after Table tab:numif):
* how much the imputation of the missing prices weighs in the robust standard
* errors.  equaids includes it (the influence of each donor on the group means
* it enters); here the same estimation is run a second time on the prices
* already filled, as data: the estimates are identical, the standard errors
* lack the term of the imputation.  Ratio = SE with / SE without, by block of
* estimates:
*   1. Mexican cereals, QUAIDS, hhsize and isMale, pimpute(psu rururb),
*      sampling weights, vce(robust);
*   2. the same with the correction for the non-buyers;
*   3. the perfect design of Section 6.4 (dgp_nonbuyers.do, 1,000 households,
*      prices shared by groups of 25 and filled by pimpute(grp)), with the
*      correction.
* Writes out/pimpute_effect.csv.
*
*   cd <repository>/replication
*   do pimpute_effect.do
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

mata:
// standard errors and block labels: coefficients e(b), probits (the tail of
// e(sel_psi)), then the four families of elasticities
real rowvector _pe_se()
{
    real rowvector t
    string scalar  f
    real scalar    k, nb, np
    t = sqrt(diagonal(st_matrix("e(V)")))'
    if (st_global("e(selection)") != "") {
        nb = cols(st_matrix("e(b_free)"))
        np = cols(st_matrix("e(sel_psi)"))
        t = t, sqrt(diagonal(st_matrix("e(V_sel_psi)")))'[|nb + 1 \ np|]
    }
    for (k = 1; k <= 4; k++) {
        f = ("", "w", "m", "h")[k]
        t = t, sqrt(diagonal(st_matrix("e(V_elas_x" + f + ")")))',
            sqrt(diagonal(st_matrix("e(V_elas_u" + f + ")")))',
            sqrt(diagonal(st_matrix("e(V_elas_c" + f + ")")))'
    }
    return(t)
}
string rowvector _pe_bl()
{
    string rowvector t, fn
    real scalar k, M, nb, np
    M = cols(st_matrix("e(elas_x)"))
    t = J(1, cols(st_matrix("e(b)")), "coefficients")
    if (st_global("e(selection)") != "") {
        nb = cols(st_matrix("e(b_free)"))
        np = cols(st_matrix("e(sel_psi)"))
        t = t, J(1, np - nb, "probits")
    }
    fn = ("market", "households", "reference", "household mean")
    for (k = 1; k <= 4; k++) t = t, J(1, M + 2 * M * M, "elasticities, " + fn[k])
    return(t)
}
void _pe_report(real rowvector s1, real rowvector s0, string rowvector bl,
    string scalar lab, string scalar csv)
{
    real rowvector r, sel
    string rowvector u
    real scalar b, j, fh
    r = s1 :/ s0
    printf("\n{txt}%s: SE with the imputation / SE with the prices as data\n", lab)
    printf("{txt}  %-34s %6s %8s %8s\n", "", "n", "min", "max")
    u = J(1, 0, "")
    for (j = 1; j <= cols(bl); j++) if (!anyof(u, bl[j])) u = u, bl[j]
    fh = fopen(csv, "a")
    for (b = 1; b <= cols(u) + 1; b++) {
        if (b <= cols(u)) sel = (bl :== u[b]) :& (s0 :> 0) :& (s0 :< .)
        else              sel = (s0 :> 0) :& (s0 :< .)
        if (!sum(sel)) continue
        printf("{txt}  %-34s {res}%6.0f %8.4f %8.4f\n", (b <= cols(u) ? u[b] : "all"), sum(sel),
            min(select(r, sel)), max(select(r, sel)))
        fput(fh, sprintf(`""%s","%s",%g,%12.6f,%12.6f"', lab, (b <= cols(u) ? u[b] : "all"),
            sum(sel), min(select(r, sel)), max(select(r, sel))))
    }
    fclose(fh)
}
end

tempname fh
file open `fh' using "out/pimpute_effect.csv", write text replace
file write `fh' "case,block,estimates,min,max" _n
file close `fh'

* one case: the estimation with pimpute(), then the same on the filled prices
capture program drop _pe_case
program define _pe_case
    args lab wvar grp cmd
    * cmd: the estimation with the placeholders @W (weight clause) and @P
    * (pimpute option)
    local c1 : subinstr local cmd "@P" "pimpute(`grp')"
    local c1 : subinstr local c1 "@W" "[pw=`wvar']"
    quietly `c1'
    mata: S1 = _pe_se(); BL = _pe_bl(); B1 = st_matrix("e(b)")
    local prices "`e(prices)'"
    tempvar es
    quietly gen byte `es' = e(sample)
    * the prices as filled at estimation
    local lps ""
    foreach p of local prices {
        tempvar l`p'
        quietly gen double `l`p'' = ln(`p')
        local lps "`lps' `l`p''"
    }
    quietly _equaids_pimpute `lps', touse(`es') wt(`wvar') groups(`grp')
    preserve
    quietly keep if `es'
    foreach p of local prices {
        quietly replace `p' = exp(`l`p'')
    }
    local c0 : subinstr local cmd "@P" ""
    local c0 : subinstr local c0 "@W" "[pw=`wvar']"
    quietly `c0'
    mata: S0 = _pe_se(); st_numscalar("__d", mreldif(st_matrix("e(b)"), B1))
    restore
    display as text _n "`lab'" as text ": N = " as result e(N) ///
        as text ", estimates with pimpute() and on the filled prices: max rel. diff. " as result %9.2e __d
    mata: _pe_report(S1, S0, BL, "`lab'", "out/pimpute_effect.csv")
end

* ---- 1 and 2: the Mexican cereals ----
use "`ROOT'/examples/mexico_2014_cereals.dta", clear
local W wcorn wwheat wrice wother wcomp
local P pcorn pwheat price pother pcomp
local O prices(`P') expenditure(hh_current_inc) demographics(hhsize isMale) vce(robust) notable nolog
_pe_case "Mexico, pimpute(psu rururb)" sweight "psu rururb" "equaids `W' @W, `O' @P"
_pe_case "Mexico, pimpute(psu rururb), selection" sweight "psu rururb" "equaids `W' @W, `O' @P selection"

* ---- 3: the perfect design of Section 6.4 ----
run "dgp_nonbuyers.do"
_dgp_nonbuyers 1000 0.3 -0.10 0.08 20260928
quietly gen double one = 1
_pe_case "design of Section 6.4, pimpute(grp), selection" one grp ///
    "equaids w1 w2 w3 @W, prices(p1 p2 p3) expenditure(x) noquadratic demographics(hs) anot(0) @P selection selvars(w1: q1 ; w2: q2) vce(robust) notable nolog"

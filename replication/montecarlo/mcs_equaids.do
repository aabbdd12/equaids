* mcs_equaids.do -- synthetic Monte Carlo of the analytic standard errors of
* equaids (design-based validation, finite population).
*
* The PERFECT case (the validation of a variance formula needs every
* assumption to hold; weak identification is the business of equaidsdiag).
*
* A finite population of 100,000 households in 10 strata x 400 PSUs x 25
* households, generated ONCE per cell (fixed seed): prices set at the PSU
* level (with small household noise), a log-normal expenditure with stratum
* and PSU effects (the strata differ in income), two demographics (an
* indicator, p = 0.4, and a household size 1 to 8, Ray scaling), budget
* shares from the model equations (QUAIDS or AIDS, equaids's own _eq_parts)
* at known parameters, plus errors drawn once (the last share by adding-up).
* The truth is the census: equaids on the whole population.  Each
* replication draws a sample by the design and estimates; the randomness is
* the sampling only.
*
* Before any replication the DGP is checked (and refused otherwise):
*   theory:  Slutsky s_ij = w_i e*_ij symmetric and negative semidefinite for
*            every household (the quadratic form of consumer theory);
*            systematic shares in [0, 1] for every household;
*   Engel:   the quadratic term is real -- turning points -beta_i/(2 lambda_i)
*            inside the range of l; the condition index of (1, l, l^2) on the
*            population below D4's warning level (100);
*   identification: every free parameter clearly identified at the sample
*            size, |theta|/SE >= 4 (census SE scaled by sqrt(N/n)).
* The parameters were tuned in pilots (tune*.do) until every check passes:
* smallest |theta|/SE 4.2 to 5.3 over the six model x error cells.
*
* alpha_0 is fixed at the population value (min ln x - 0.1) in the census and
* in every sample: the check is of the standard-error formula; the alpha_0 rule
* was validated separately (bootstrap, replication/bootstrap/).
*
* Cells: model (quaids|aids), errors (homo|hetero|clus), design
*   srs  n households drawn independently (with replacement), no weights,
*        vce(robust) -- what the robust variance assumes
*   svy  one-stage stratified cluster sample: n/25 PSUs drawn without
*        replacement within strata (poorer strata oversampled), all their
*        households; weights = inverse probabilities; vce(svy) with the fpc
*        (exact, no second stage to approximate)
* n = 5000.
*
*   do mcs_equaids.do  model  err  design  n  R  part  [small]
args model err design n R part small
version 14.2
clear all
set more off
local R     = cond("`R'" == "", 50, real("`R'"))
local part  = cond("`part'" == "", 1, real("`part'"))
local n     = cond("`n'" == "", 5000, real("`n'"))
local small = cond("`small'" == "", 0, real("`small'"))
* ---- run this from replication/montecarlo/ ----
local ROOT = subinstr("`c(pwd)'", "\", "/", .) + "/../.."
capture confirm file "`ROOT'/src/equaids.ado"
if _rc {
    display as error "run this script from the replication/montecarlo/ directory"
    exit 601
}
adopath ++ "`ROOT'/src"
run "`ROOT'/src/equaids.ado"
mata:
// free parameters of equaids from the full vector: least squares on the
// affine map full = c0 + th * D' (exact, the full vector satisfies the
// restrictions)
real rowvector _mcs_free(real rowvector f, real scalar M, real scalar K, real scalar qd)
{
    real scalar    P, r
    real rowvector c0, e
    real matrix    D
    P  = 2 * (M - 1) + (M - 1) * M / 2 + qd * (M - 1) + K * (M - 1) + K
    c0 = _eq_full(J(1, P, 0), M, K, qd)
    D  = J(cols(c0), P, 0)
    for (r = 1; r <= P; r++) {
        e = J(1, P, 0) ; e[r] = 1
        D[., r] = (_eq_full(e, M, K, qd) - c0)'
    }
    return(((f - c0) * D) * invsym(D' * D))
}
// the full parameter vector of the DGP written as name,value lines
void _mcs_csv(string scalar fn, real rowvector f, real scalar M, real scalar K, real scalar qd)
{
    real scalar    fh, i, j, k
    string rowvector nm
    nm = ("alpha" :+ strofreal(1..M)), ("beta" :+ strofreal(1..M))
    for (j = 1; j <= M; j++) for (i = j; i <= M; i++) nm = nm, sprintf("gamma_%g_%g", i, j)
    if (qd) nm = nm, ("lambda" :+ strofreal(1..M))
    for (k = 1; k <= K; k++) nm = nm, (sprintf("eta_z%g_", k) :+ strofreal(1..M))
    nm = nm, ("rho_z" :+ strofreal(1..K))
    if (fileexists(fn)) unlink(fn)
    fh = fopen(fn, "w")
    fput(fh, "parameter,value")
    for (j = 1; j <= cols(f); j++) fput(fh, nm[j] + "," + strofreal(f[j], "%21.0g"))
    fclose(fh)
}
end
local OUT "raw"
* $MCS_OUT: another folder, e.g. for a short check against raw/
if "$MCS_OUT" != "" local OUT "$MCS_OUT"
capture mkdir "`OUT'"
local qd = ("`model'" == "quaids")
local nq = cond(`qd', "", "noquadratic")
* households per PSU: 25 in the final DGP (200 PSUs sampled, 10 to 32 per
* stratum); with 50 (100 PSUs, 5 to 16 per stratum) the linearized design
* variance is a little too small -- the known small-sample limit of the
* linearization with few PSUs per stratum (mcs_summary.do)
local psz = cond("$MCS_PSUSIZE" == "", 25, real("$MCS_PSUSIZE"))
local tag = "`model'_`err'_`design'_n`n'" + cond(`small', "_small", "") + cond(`psz' != 50, "_psu`psz'", "")

* ---- true parameters (full vector, Poi's order) ----
* mean shares 40/30/20/10 (small: 50/30/18/2); gamma proportional to
* (wbar wbar' - diag(wbar)): symmetric, rows summing to zero, negative
* semidefinite -- the Cobb-Douglas Slutsky shape
mata:
    M = 4 ; K = 2
    wb = (`small' ? (.50, .30, .18, .02) : (.40, .30, .20, .10))
    al = wb
    // values tuned in pilots (tune.do) so that every free parameter is
    // clearly identified at n = 5000 (|theta|/SE >= 4.5): Ray's rho acts
    // only through ln m0 inside l, its effect on the shares scales with beta
    // and lambda, hence Engel curves with real slopes and curvature
    be = (`small' ? (-.06, .03, .028, .002) : (-.15, .05, .06, .04))
    la = (`small' ? (.008, -.004, -.0038, -.0002) : (.02, -.008, -.007, -.005))
    gs = .3
    et = (.03, -.012, -.010, -.008) \ (.01, -.004, -.003, -.003)
    rh = (1.2, .3)
    // AIDS: linear Engel curves, so gentler slopes keep the shares in [0,1]
    // and Slutsky negative semidefinite over the range of l; rho, then
    // identified through beta alone, needs a larger first element
    if (!`qd') be = (-.06, .02, .025, .015)
    if (!`qd') rh = (1.6, .3)
    // overrides for the pilots that tune the DGP (globals MCS_*): each is a
    // list of numbers; eta rows sum to zero, beta and lambda too
    if (st_global("MCS_BETA") != "") be = strtoreal(tokens(st_global("MCS_BETA")))
    if (st_global("MCS_LAMBDA") != "") la = strtoreal(tokens(st_global("MCS_LAMBDA")))
    if (st_global("MCS_GSCALE") != "") gs = strtoreal(st_global("MCS_GSCALE"))
    if (st_global("MCS_ETA1") != "") et[1, .] = strtoreal(tokens(st_global("MCS_ETA1")))
    if (st_global("MCS_ETA2") != "") et[2, .] = strtoreal(tokens(st_global("MCS_ETA2")))
    if (st_global("MCS_RHO") != "") rh = strtoreal(tokens(st_global("MCS_RHO")))
    Ga = gs :* (wb' * wb - diag(wb))
    if (!`qd') la = J(1, M, 0)
    // Poi's full vector: alpha, beta, gamma (lower triangle by column),
    // lambda (QUAIDS), eta (by demographic), rho
    f = al, be
    for (j = 1; j <= M; j++) for (i = j; i <= M; i++) f = f, Ga[i, j]
    if (`qd') f = f, la
    f = f, vec(et')', rh
    st_matrix("F0", f)
end

* ---- finite population ----
set seed `=20260926 + 100 * `qd' + 10 * `small''
set obs 100000
gen long hh = _n
gen int stratum = ceil(hh / 10000)
gen long psu = ceil(hh / `psz')
* income differs by stratum and PSU; prices vary by PSU (spatial variation)
* with small household noise
bysort psu: gen double _pe = rnormal(0, .1) if _n == 1
bysort psu (_pe): replace _pe = _pe[1]
* a large group of poor households concentrated at the bottom (50% overall,
* 70% in the poorest stratum to 30% in the richest): ln x lower by 1.7, with
* a smaller dispersion (0.15, against 0.6 for the others; tuned so that the
* condition index of (1, l, l^2) is about that of Poi's data with the default
* alpha_0, 23).  With a single normal of 200,000 draws the smallest ln x
* lies 4.5 sd below the mean, isolated from the mass, and alpha_0 (set just
* below it) leaves l far from 0 relative to its spread: a conditioning of
* (1, l, l^2) worse than on real data (23 to 36 on three data sets,
* the diagnostics of equaidsdiag).  A dense bottom of the distribution brings the minimum close
* to the mass, as in the surveys of low-income countries.
* by default a single log-normal expenditure (the classical assumption),
* ln x = 5.5 + stratum + PSU effects + N(0, 0.6); MCS_INCOME = mixture gives
* the mixture of two log-normals above instead
gen byte poor = runiform() < .7 - .4 * (stratum - 1) / 9
if "$MCS_INCOME" != "mixture" gen double lx = 5.5 + .12 * (stratum - 5.5) + _pe + rnormal(0, .6)
else gen double lx = 5.5 + .12 * (stratum - 5.5) + _pe + cond(poor, -1.7 + rnormal(0, .15), rnormal(0, .6))
forvalues k = 1/4 {
    bysort psu: gen double _pk = rnormal(0, .2) if _n == 1
    bysort psu (_pk): replace _pk = _pk[1]
    gen double lp`k' = _pk + rnormal(0, .05)
    drop _pk
}
drop _pe
sort hh
gen byte z1 = runiform() < .4
gen byte z2 = 1 + floor(8 * runiform())
quietly summarize lx
local a0 = r(min) - .1
display as text "population: N = 100,000, alpha_0 = " %8.4f `a0' ", model `model', small `small'"

* systematic shares, elasticities of every household
forvalues k = 1/4 {
    gen double f`k' = .
}
mata:
    LP = st_data(., "lp1 lp2 lp3 lp4") ; LX = st_data(., "lx") ; Z = st_data(., "z1 z2")
    // free parameters from the full vector: equaids's own map
    th = _mcs_free(st_matrix("F0"), M, K, `qd')
    Fh = . ; MU = . ; MUJ = .
    _eq_parts(th, LP, LX, Z, `a0', `qd', M, Fh, MU, MUJ)
    st_store(., ("f1", "f2", "f3", "f4"), Fh)
    // theory: Slutsky s_ij = w_i e*_ij = mu_ij + w_i w_j + mu_i w_j - delta_ij w_i
    // (with e_ij = -delta_ij + mu_ij/w_i, e_i = 1 + mu_i/w_i), for every household
    N = rows(LP) ; nbad = 0 ; asym = 0 ; maxev = -.
    for (h = 1; h <= N; h++) {
        w = Fh[h, .]
        Sl = J(M, M, .)
        for (i = 1; i <= M; i++) for (j = 1; j <= M; j++)
            Sl[i, j] = MUJ[h, (i - 1) * M + j] + w[i] * w[j] + MU[h, i] * w[j] - (i == j) * w[i]
        asym = max((asym, max(abs(Sl - Sl'))))
        ev = symeigenvalues((Sl + Sl') :/ 2)
        if (max(ev) > 1e-10) nbad++
        maxev = max((maxev, max(ev)))
    }
    printf("{txt}theory: Slutsky negative semidefinite for {res}%9.5f{txt} of the households (largest eigenvalue %9.2e; asymmetry %9.2e)\n", 1 - nbad / N, maxev, asym)
    st_numscalar("nsd", 1 - nbad / N)
    st_numscalar("shin", mean(rowmin(Fh) :>= 0 :& rowmax(Fh) :<= 1))
end
* Engel: turning points of the quadratic term and conditioning of (1, l, l^2)
quietly equaidsdiag f1 f2 f3 f4, lnprices(lp1-lp4) lnexpenditure(lx) demographics(z1 z2) anot(`a0') `nq'
local cq = r(cond_quad)
mata:
    fl = st_matrix("F0")
    ell = LX :- `a0' :- LP * wb'
    so = sort(ell, 1)
    printf("{txt}Engel: l from %6.2f to %6.2f (5%%-95%%: %6.2f to %6.2f)\n", so[1], so[N], so[ceil(.05 * N)], so[ceil(.95 * N)])
    // (no if-block: interactive Mata then waits for an else and reads -end-)
    for (i = 1; i <= M * `qd'; i++) printf("{txt}   good %g: beta %7.4f lambda %8.5f, turning point of the share at l = %8.2f\n", i, be[i], la[i], -be[i] / (2 * la[i]))
end
display as text "Engel: condition index of (1, l, l^2) on the population: " as result %6.1f `cq'
display as text "shares: systematic shares in [0,1] for " as result %9.5f shin as text " of the households"
* the quadratic block is refused only above D4's warning level (100): the
* binding identification criterion is |theta|/SE >= 4 below (lambda
* included); a log-normal expenditure with a large population puts the
* smallest ln x about 4.4 sd below the mean and the index near 50
if (nsd < 1 | shin < 1 | (`qd' & `cq' > 100)) {
    display as error "the DGP violates theory or identification: no replication"
    exit 459
}

* errors drawn once (the population's shares); scale 0.3 x mean share, the
* correlations of a multinomial; heteroskedastic: sd times exp((ln x - mean)/2)
* (mean variance kept); clustered: 30% of the variance at the PSU level
mata:
    es = (st_global("MCS_ESCALE") != "" ? strtoreal(st_global("MCS_ESCALE")) : .2)
    D = diag(es :* wb[|1 \ M - 1|])
    Cm = diag(wb) - wb' * wb
    Cm = Cm[|1, 1 \ M - 1, M - 1|]
    Cr = Cm :/ sqrt(diagonal(Cm) * diagonal(Cm)')
    C = cholesky(D * Cr * D)
    E = rnormal(N, M - 1, 0, 1) * C'
    if ("`err'" == "hetero") {
        zz = LX :- mean(LX)
        E = E :* (exp(zz :/ 2) :/ sqrt(mean(exp(zz))))
    }
    if ("`err'" == "clus") {
        ps = st_data(., "psu")
        Ec = rnormal(max(ps), M - 1, 0, 1) * C'
        E = sqrt(.7) :* E + sqrt(.3) :* Ec[ps, .]
    }
    Wo = Fh[|1, 1 \ N, M - 1|] + E
    Wo = Wo, 1 :- rowsum(Wo)
    st_store(., ("f1", "f2", "f3", "f4"), Wo)
    st_numscalar("win", mean(rowmin(Wo) :>= 0 :& rowmax(Wo) :<= 1))
end
rename (f1 f2 f3 f4) (w1 w2 w3 w4)
* the population's signature: whoever regenerates it checks it is the same
quietly datasignature
local dsig "`r(datasignature)'"
display as text "population signature: " as result "`dsig'"
display as text "observed shares in [0,1] for " as result %9.5f win as text " of the households"
local EQ "w1 w2 w3 w4, lnprices(lp1-lp4) lnexpenditure(lx) demographics(z1 z2) anot(`a0') `nq' nolog notable"

* ---- the census: the truth ----
quietly equaids `EQ' tolerance(1e-10) noelastse
matrix T0 = e(elas_x), vec(e(elas_u)')', vec(e(elas_c)')'
matrix TH0 = e(b_free)
local K = colsof(T0)
local P = colsof(TH0)
display as text "census: converged " e(converged) ", " e(iter) " iterations"
* identification: every free parameter clearly identified at the sample
* size, |theta| / SE_n >= 4, with SE_n the robust census standard error
* scaled by sqrt(N/n) (a perfect-case DGP: weak identification is the
* business of equaidsdiag, not of this validation)
mata:
    bt = st_matrix("TH0")
    se = sqrt(diagonal(st_matrix("e(V_free)")))' :* sqrt(rows(LP) / `n')
    tr = abs(bt) :/ se
    nm = ("alpha" :+ strofreal(1..M-1)), ("beta" :+ strofreal(1..M-1))
    for (j = 1; j <= M - 1; j++) for (i = j; i <= M - 1; i++) nm = nm, sprintf("gamma%g%g", i, j)
    if (`qd') nm = nm, ("lambda" :+ strofreal(1..M-1))
    nm = nm, ("eta1_" :+ strofreal(1..M-1)), ("eta2_" :+ strofreal(1..M-1)), ("rho1", "rho2")
    printf("{txt}identification at n = `n': |theta| / SE\n")
    for (j = 1; j <= cols(bt); j++) printf("{txt}   %-10s %10.5f  SE %9.5f  {res}%7.1f%s\n", nm[j], bt[j], se[j], tr[j], (tr[j] < 4 ? "  <-- weak" : ""))
    st_numscalar("tmin", min(tr))
end
display as text "identification: smallest |theta|/SE at n = `n': " as result %6.1f tmin
if (tmin < 4 | `R' == 0) {
    if (tmin < 4) display as error "the DGP is not a perfect case (a parameter is weakly identified): no replication"
    exit
}
if `part' == 1 {
    * the true parameters of the DGP (full vector, Poi's order), for replication
    mata: _mcs_csv("`OUT'/mcs_`tag'_dgp.csv", st_matrix("F0"), M, K, `qd')
    preserve
    clear
    svmat double T0, names(e)
    svmat double TH0, names(t)
    gen double nsd = nsd
    gen double cq = `cq'
    gen str80 datasignature = "`dsig'"
    save "`OUT'/mcs_`tag'_truth.dta", replace
    restore
}

* ---- design ----
* one-stage stratified cluster sample: n/50 PSUs drawn without replacement
* within strata (the poorer strata oversampled), all their households; the
* first-stage fpc is then exact (no second stage to approximate)
local nps = `n' / `psz'
local alloc "16 14 12 12 10 10 8 8 5 5"
gen int npsu_h = .
forvalues h = 1/10 {
    quietly replace npsu_h = round(`: word `h' of `alloc'' * `nps' / 100) if stratum == `h'
}
gen int NPSU = 10000 / `psz'

tempname pf
local vl ""
forvalues k = 1/`K' {
    local vl `vl' double e`k' double s`k'
}
forvalues k = 1/`P' {
    local vl `vl' double t`k' double st`k'
}
postfile `pf' int rep byte conv int iter double time `vl' using "`OUT'/mcs_`tag'_`part'.dta", replace
set seed `=20261000 + 1000 * `part' + 7 * `small''
timer clear 1
timer on 1
forvalues r = 1/`R' {
    preserve
    if "`design'" == "srs" {
        * independent draws (with replacement): what the robust variance assumes
        bsample `n'
        local V ""
    }
    else {
        * first stage: PSUs without replacement within strata
        bysort psu: gen double _u = runiform() if _n == 1
        * rank of each PSU within its stratum, by a uniform draw
        bysort stratum (_u): gen long _rk = _n if !missing(_u)
        bysort psu (_rk): replace _rk = _rk[1]
        quietly keep if _rk <= npsu_h
        gen double pw = NPSU / npsu_h
        quietly svyset psu [pw=pw], strata(stratum) fpc(NPSU)
        local V "vce(svy)"
    }
    capture quietly equaids `EQ' `V'
    if _rc {
        local line "(`r') (0) (0) (.)"
        forvalues k = 1/`=2 * (`K' + `P')' {
            local line "`line' (.)"
        }
    }
    else {
        matrix E = e(elas_x), vec(e(elas_u)')', vec(e(elas_c)')'
        mata: st_matrix("S", (sqrt(diagonal(st_matrix("e(V_elas_x)")))', ///
            sqrt(diagonal(st_matrix("e(V_elas_u)")))', sqrt(diagonal(st_matrix("e(V_elas_c)")))'))
        matrix T = e(b_free)
        mata: st_matrix("ST", sqrt(diagonal(st_matrix("e(V_free)")))')
        local line "(`r') (`=e(converged)') (`=e(iter)') (`=e(time)')"
        forvalues k = 1/`K' {
            local line "`line' (`=E[1, `k']') (`=S[1, `k']')"
        }
        forvalues k = 1/`P' {
            local line "`line' (`=T[1, `k']') (`=ST[1, `k']')"
        }
    }
    post `pf' `line'
    restore
    if mod(`r', 10) == 0 {
        timer off 1
        quietly timer list 1
        display as text "`tag' part `part': " as result `r' as text " of `R' replications, " ///
            as result %6.0f r(t1) as text " s"
        timer on 1
    }
}
postclose `pf'


*! equaidsdiag 1.0.0  2026-09-25  Abdelkrim Araar
*! Diagnostics of an AIDS/QUAIDS specification before estimating it: data,
*! prices, expenditure and alpha_0, demographics, conditioning of the design
*! at the starting point, small goods; optionally the sensitivity of the
*! estimates to alpha_0 (D6) and their stability to each demographic (D7),
*! the only parts that estimate.
*!
*! Principle (as easidiag): everything in D0-D5 is computed from the data and
*! from the regressors at the starting point -- alpha at the mean shares, the
*! other parameters at zero, so that ln a(p) is the Stone index plus alpha_0
*! and l = ln x - alpha_0 - sum_k wbar_k ln p_k.  A diagnostic that needed the
*! model to converge would be silent exactly when it does not.
program define equaidsdiag, rclass
    version 14.2
    syntax varlist(min=3 numeric) [if] [in] [aweight fweight pweight iweight], ///
        [ PRices(varlist numeric) LNPRices(varlist numeric)                    ///
          EXPenditure(varname numeric) LNEXPenditure(varname numeric)          ///
          DEMOgraphics(varlist numeric) ANOT(string) noQUadratic               ///
          SENSitivity A0list(numlist) STABility * ]

    local shares `varlist'
    local M : word count `shares'
    if ("`prices'" != "") + ("`lnprices'" != "") != 1 {
        di as err "specify one of prices() or lnprices()"
        exit 198
    }
    if ("`expenditure'" != "") + ("`lnexpenditure'" != "") != 1 {
        di as err "specify one of expenditure() or lnexpenditure()"
        exit 198
    }
    local pv `prices'`lnprices'
    if `: word count `pv'' != `M' {
        di as err "`M' goods but `: word count `pv'' price variables"
        exit 198
    }
    local quad = ("`quadratic'" == "")
    local ignored `"`options'"'

    * ---- households lost to missing values (D0) ----
    marksample all, novarlist
    quietly count if `all'
    local N0 = r(N)
    local lost ""
    foreach v in `shares' `pv' `expenditure' `lnexpenditure' `demographics' {
        quietly count if `all' & missing(`v')
        if r(N) local lost `"`lost' `v' `r(N)'"'
    }
    marksample touse
    markout `touse' `pv' `expenditure' `lnexpenditure' `demographics'
    if "`weight'" != "" {
        tempvar wt
        quietly gen double `wt' `exp' if `touse'
        quietly replace `touse' = 0 if `wt' <= 0 | missing(`wt')
    }
    quietly count if `touse'
    if r(N) == 0 error 2000
    local N = r(N)

    * ---- logs, weights (normalized as in equaids) ----
    local lnp ""
    foreach v of local pv {
        if "`prices'" != "" {
            quietly count if `touse' & `v' <= 0
            if r(N) {
                di as err "`v' has nonpositive values"
                exit 411
            }
            tempvar l`v'
            quietly gen double `l`v'' = ln(`v') if `touse'
            local lnp `lnp' `l`v''
        }
        else local lnp `lnp' `v'
    }
    if "`expenditure'" != "" {
        tempvar lnx
        quietly gen double `lnx' = ln(`expenditure') if `touse'
    }
    else local lnx `lnexpenditure'
    tempvar om
    if "`weight'" == "" quietly gen double `om' = 1 if `touse'
    else {
        quietly gen double `om' = `wt' if `touse'
        if inlist("`weight'", "aweight", "pweight") {
            quietly summarize `om' if `touse', meanonly
            quietly replace `om' = `om' * `N' / r(sum) if `touse'
        }
    }
    quietly summarize `lnx' if `touse', meanonly
    local lnxmin = r(min)
    if `"`anot'"' == "" {
        local a0 = `lnxmin' - 0.1
        local a0rule "min ln x - 0.1 (equaids's default)"
    }
    else {
        capture confirm number `anot'
        if _rc {
            di as err "anot() must be a number"
            exit 198
        }
        local a0 = `anot'
        local a0rule "anot()"
    }

    * ---- header ----
    di _n as txt "{hline 78}"
    di as txt "equaidsdiag: " as res cond(`quad', "QUAIDS", "AIDS") as txt " specification, " ///
        as res `M' as txt " goods, " as res `: word count `demographics'' as txt " demographics"
    * N on the line of alpha_0: the line above passes column 56 with demographics
    di as txt "alpha_0 = " as res %8.4f `a0' as txt "  (`a0rule')" ///
        _col(60) as txt "N = " as res %12.0fc `N'
    di as txt "Nothing is estimated in D0-D5: the data and the regressors at the starting"
    di as txt "point (Stone index) only."
    if `"`ignored'"' != "" di as txt "Options of the estimator ignored: " as res `"`ignored'"'
    di as txt "{hline 78}"
    if `N0' > `N' {
        di _n as txt "Households lost to missing values or nonpositive weights: " ///
            as res `N0' - `N' as txt " of " as res `N0' ///
            as txt " (" as res %4.1f 100 * (`N0' - `N') / `N0' as txt "%)"
        if `"`lost'"' != "" {
            local l `"`lost'"'
            while `"`l'"' != "" {
                gettoken v l : l
                gettoken n l : l
                di as txt "   `v'" _col(30) "missing " as res %10.0fc `n'
            }
        }
    }

    * ---- D0-D5 ----
    tempname SH BKW DS
    mata: _eqd_main("`shares'", "`lnp'", "`lnx'", "`demographics'", "`om'", ///
        "`touse'", `a0', `quad', "`SH'", "`BKW'", "`DS'")
    local nwarn = `nwarn' + (`N0' - `N' > 0.05 * `N0')

    * the model's options, for D6 and D7 (the only sections that estimate)
    local wopt = cond("`weight'" != "", "[`weight'`exp']", "")
    local popt = cond("`prices'" != "", "prices(`prices')", "lnprices(`lnprices')")
    local xopt = cond("`expenditure'" != "", "expenditure(`expenditure')", "lnexpenditure(`lnexpenditure')")

    * ---- D6: sensitivity to alpha_0 (estimates) ----
    tempname SENS
    if "`sensitivity'" != "" | "`a0list'" != "" {
        if "`a0list'" == "" {
            local d = `lnxmin' - 0.1
            local a0list "`d' `=`d' - 1' `=`d' - 2' `=`d' - 4'"
        }
        local na : word count `a0list'
        matrix `SENS' = J(`na', 5 + `M', .)
        local enames ""
        foreach s of local shares {
            local enames `enames' E_`s'
        }
        matrix colnames `SENS' = alpha_0 converged iter ll rcond `enames'
        local dopt = cond("`demographics'" != "", "demographics(`demographics')", "")
        di _n as txt "D6. Sensitivity to alpha_0 (each row estimates the model)"
        di as txt "{hline 78}"
        di as txt %10s "alpha_0" %6s "conv" %6s "iter" %14s "log lik." %11s "rcond" ///
            "   aggregate expenditure elasticities"
        local r 0
        foreach a of local a0list {
            local ++r
            capture quietly equaids `shares' `wopt' if `touse', `popt' `xopt' `dopt' ///
                anot(`a') `quadratic' noelastse notable nolog
            if _rc {
                matrix `SENS'[`r', 1] = `a'
                di as txt %10.4f `a' as err "   estimation failed (r(" _rc "))"
                continue
            }
            tempname EX
            matrix `EX' = e(elas_x)
            matrix `SENS'[`r', 1] = `a'
            matrix `SENS'[`r', 2] = e(converged)
            matrix `SENS'[`r', 3] = e(iter)
            matrix `SENS'[`r', 4] = e(ll)
            matrix `SENS'[`r', 5] = e(rcond)
            local line ""
            forvalues i = 1/`M' {
                matrix `SENS'[`r', 5 + `i'] = `EX'[1, `i']
                local line "`line'`: di %8.3f `EX'[1, `i']'"
            }
            di as res %10.4f `a' %6.0f e(converged) %6.0f e(iter) %14.2f e(ll) %11.2e e(rcond) as res "`line'"
        }
        di as txt "{hline 78}"
        di as txt "A scaled rcond below about 1e-5 marks a nearly singular information matrix:"
        di as txt "the sandwich standard errors of the coefficients are then unreliable (the"
        di as txt "bootstrap of Poi's data at alpha_0 = 10 gives 2 to 4 times larger ones)."
        return matrix sensitivity = `SENS'
    }

    * ---- D7: stability to the demographics (estimates) ----
    * Each demographic left out in turn, from the full model and on the same
    * sample: the change dE of the aggregate expenditure and own-price
    * elasticities, its standard error from the difference of the two
    * estimates' influence functions (robust), and z = dE / s.e.  The 2M
    * changes of a demographic are judged at the Bonferroni level (5% two-
    * sided over 2M tests): leaving z out then moves an elasticity beyond
    * sampling noise -- the demographic belongs in the model, often through
    * its correlation with prices.  (A first-order approximation without
    * re-estimating was tried and dropped: it fails exactly when a variable
    * matters, the change being then far from small.)
    tempname STAB
    local K : word count `demographics'
    if "`stability'" != "" & `K' == 0 di _n as txt "D7. Stability: no demographic variable to leave out"
    if "`stability'" != "" & `K' > 0 {
        local zc = invnormal(1 - .05 / (2 * 2 * `M'))
        di _n as txt "D7. Stability: each demographic left out in turn (the model estimated again)"
        di as txt "{hline 78}"
        capture drop _d7f*
        capture drop _d7r*
        quietly equaids `shares' `wopt' if `touse', `popt' `xopt' demographics(`demographics') ///
            anot(`a0') `quadratic' notable nolog saveif(_d7f)
        tempname EF UF SX SU EX UX
        matrix `EF' = e(elas_x)
        matrix `UF' = e(elas_u)
        matrix `STAB' = J(`K' * `M', 6, .)
        matrix colnames `STAB' = dE_x se_x z_x dE_ii se_ii z_ii
        local rn ""
        local k 0
        local zmax 0
        local nsig 0
        foreach z of local demographics {
            local ++k
            local others : list demographics - z
            local dopt = cond("`others'" != "", "demographics(`others')", "")
            capture quietly equaids `shares' `wopt' if `touse', `popt' `xopt' `dopt' ///
                anot(`a0') `quadratic' notable nolog saveif(_d7r)
            if _rc {
                di as err "   `z': the model without it could not be estimated (r(" _rc "))"
                capture drop _d7r*
                continue
            }
            matrix `EX' = e(elas_x)
            matrix `UX' = e(elas_u)
            mata: _eqd_d7se("`touse'", `M', "`SX'", "`SU'")
            di as txt _n "Without " as res "`z'" as txt ":"
            di as txt %-16s "good" %12s "dE_x" %9s "s.e." %8s "z" %14s "dE_ii" %9s "s.e." %8s "z"
            local i 0
            foreach s of local shares {
                local ++i
                local r = (`k' - 1) * `M' + `i'
                local dxe = `EX'[1, `i'] - `EF'[1, `i']
                local due = `UX'[`i', `i'] - `UF'[`i', `i']
                local zx = `dxe' / `SX'[1, `i']
                local zu = `due' / `SU'[`i', `i']
                matrix `STAB'[`r', 1] = (`dxe', `SX'[1, `i'], `zx', `due', `SU'[`i', `i'], `zu')
                local zmax = max(`zmax', abs(`zx'), abs(`zu'))
                local nsig = `nsig' + (abs(`zx') > `zc') + (abs(`zu') > `zc')
                local rn `rn' `z':`s'
                local fx = cond(abs(`zx') > `zc', "*", " ")
                local fu = cond(abs(`zu') > `zc', "*", " ")
                di as txt %-16s abbrev("`s'", 16) as res %12.4f `dxe' %9.4f `SX'[1, `i'] ///
                    %8.2f `zx' as txt "`fx'" as res %13.4f `due' %9.4f `SU'[`i', `i'] %8.2f `zu' as txt "`fu'"
            }
            capture drop _d7r*
        }
        capture drop _d7f*
        matrix rownames `STAB' = `rn'
        di as txt "{hline 78}"
        di as txt "dE: change of the aggregate elasticity (expenditure; own price, ii) when the"
        di as txt "demographic is left out, on the same sample; s.e. robust, from the two sets"
        di as txt "of influence functions.  * : |z| > " %4.2f `zc' " (5%, Bonferroni over the " 2 * `M' " changes"
        di as txt "of a demographic)."
        if `nsig' > 0 {
            di as txt "note: leaving a demographic out moves " as res `nsig' as txt " elasticities beyond sampling noise:"
            di as txt "      the elasticities depend on it.  Most often it belongs in the model"
            di as txt "      (it moves with expenditure or prices, D3); but check also its rho: a"
            di as txt "      fit where 1 + rho'z nears zero for a few households moves them as well"
        }
        return scalar stab_zmax = `zmax'
        return scalar stab_nsig = `nsig'
        return scalar stab_zcrit = `zc'
        return matrix stability = `STAB'
    }

    di _n as txt "{hline 78}"
    if `nwarn' == 0 di as txt "No warning."
    else di as txt "Warnings: " as res `nwarn'
    return scalar N        = `N'
    return scalar N_lost   = `N0' - `N'
    return scalar N_warn   = `nwarn'
    return scalar anot     = `a0'
    return scalar cond_max = `condmax'
    return scalar cond_quad = `condquad'
    return scalar n_l0neg  = `nl0neg'
    return matrix shares   = `SH'
    return matrix bkw      = `BKW'
    if "`demographics'" != "" return matrix demo = `DS'
end

mata:
// ---------------------------------------------------------------------------
// D0-D5.  Prints the tables and warnings; locals nwarn, condmax, condquad,
// nl0neg; matrices shn (shares), bkwn (condition indexes and
// variance-decomposition proportions), dsn (demographics).
// ---------------------------------------------------------------------------
void _eqd_main(string scalar wv, string scalar lpv, string scalar lxv,
               string scalar zv, string scalar omv, string scalar touse,
               real scalar a0, real scalar qd, string scalar shn,
               string scalar bkwn, string scalar dsn)
{
    real matrix    W, LP, Z, R, SH, X, B, DS, C, Cr, Rp, Cz
    real rowvector cz
    real colvector om, w, lx, x, l0, med, mad, z, eta, etaall, etaq
    real rowvector wbar, S, sd, vif
    real scalar    N, M, K, j, k, nw, nout, nd, pr, lo, cmax, cq, nl, ss
    string rowvector wn, zn, xn
    string scalar  lst

    W  = st_data(., tokens(wv), touse)
    LP = st_data(., tokens(lpv), touse)
    lx = st_data(., lxv, touse)
    om = st_data(., omv, touse)
    wn = tokens(wv)
    N = rows(W) ; M = cols(W)
    if (zv != "") {
        Z = st_data(., tokens(zv), touse)
        zn = tokens(zv)
    }
    else {
        Z = J(N, 0, .)
        zn = J(1, 0, "")
    }
    K = cols(Z)
    w = om :/ sum(om)
    nw = 0

    // ---- D0: shares ----
    wbar = colsum(w :* W)
    S    = colsum(w :* exp(lx) :* W) :/ sum(w :* exp(lx))
    SH = J(M, 5, .)
    printf("\n{txt}D0. Budget shares\n{hline 78}\n")
    printf("{txt}%-16s %11s %11s %9s %9s %9s\n", "good", "mean share", "aggregate", "% zero", "% < 0", "% > 1")
    for (j = 1; j <= M; j++) {
        SH[j, .] = (wbar[j], S[j], 100 * sum(w :* (W[., j] :== 0)),
                    100 * sum(w :* (W[., j] :< 0)), 100 * sum(w :* (W[., j] :> 1)))
        printf("{txt}%-16s {res}%11.4f %11.4f %9.1f %9.1f %9.1f\n", abbrev(wn[j], 16),
               SH[j, 1], SH[j, 2], SH[j, 3], SH[j, 4], SH[j, 5])
    }
    printf("{txt}(aggregate = share of the good in total expenditure, sum x w / sum x)\n")
    st_matrix(shn, SH)
    st_matrixrowstripe(shn, (J(M, 1, ""), wn'))
    st_matrixcolstripe(shn, (J(5, 1, ""), ("mean" \ "aggregate" \ "pct_zero" \ "pct_neg" \ "pct_gt1")))
    nd = sum(abs(rowsum(W) :- 1) :> 1e-4)
    if (nd > 0) {
        printf("{err}warning: the shares do not sum to one for %g households (equaids refuses them)\n", nd)
        nw++
    }
    for (j = 1; j <= M; j++) {
        if (SH[j, 3] > 5) {
            printf("{err}warning: %s: %4.1f%% of the shares are zero.  FGNLS does not model censoring,\n", wn[j], SH[j, 3])
            printf("{err}         and the prices of the non-buyers are not observed: they have been\n")
            printf("{err}         filled in (imputed, cluster or regional unit values); check how.\n")
            nw++
        }
        if (SH[j, 4] + SH[j, 5] > 0) {
            printf("{err}warning: %s: shares outside [0, 1]\n", wn[j])
            nw++
        }
    }

    // ---- D1: prices ----
    printf("\n{txt}D1. Relative prices ln(p_k/p_M)\n{hline 78}\n")
    X  = LP[|1, 1 \ N, M - 1|] :- LP[., M]
    C  = _eqd_wcov(X, w)
    sd = sqrt(diagonal(C))'
    printf("{txt}%-16s %11s\n", "good k", "s.d.")
    for (j = 1; j < M; j++) printf("{txt}%-16s {res}%11.4f\n", abbrev(wn[j], 16), sd[j])
    lst = ""
    for (j = 1; j < M; j++) if (sd[j] < 0.05) lst = lst + " " + wn[j]
    if (lst != "") {
        printf("{err}warning: weak variation of relative prices (s.d. below 0.05):%s;\n", lst)
        printf("{err}         the price parameters of these goods are weakly identified\n")
        nw++
    }
    if (M > 2 & min(sd) > 0) {
        Cr = C :/ (sd' * sd)
        Cr = Cr - I(M - 1)
        printf("{txt}largest correlation between two relative prices: {res}%6.3f\n", max(abs(Cr)))
        if (max(abs(Cr)) > 0.95) {
            printf("{err}warning: two relative prices move almost together (|corr| > 0.95)\n")
            nw++
        }
    }
    nout = 0
    for (j = 1; j <= M; j++) {
        x = LP[., j]
        med = _eqd_median(x)
        mad = 1.4826 * _eqd_median(abs(x :- med))
        if (mad > 0) nout = nout + sum(abs(x :- med) :> 5 * mad)
    }
    printf("{txt}log prices more than 5 robust s.d. from their median: {res}%g\n", nout)
    if (nout > 0) nw++

    // ---- D2: expenditure and alpha_0 ----
    // l at the starting point: ln x - alpha_0 - sum_k wbar_k ln p_k
    l0 = lx :- a0 :- LP * wbar'
    nl = sum(l0 :<= 0)
    printf("\n{txt}D2. Expenditure and alpha_0\n{hline 78}\n")
    printf("{txt}ln x: min {res}%8.3f{txt}, median {res}%8.3f{txt}, max {res}%8.3f{txt};  alpha_0 {res}%8.4f\n",
           min(lx), _eqd_median(lx), max(lx), a0)
    printf("{txt}deflated expenditure at the start, l = ln x - alpha_0 - sum wbar_k ln p_k:\n")
    printf("{txt}   min {res}%8.3f{txt}, max {res}%8.3f{txt}; households with l <= 0: {res}%g{txt} (%4.1f%%)\n",
           min(l0), max(l0), nl, 100 * nl / N)
    if (nl > N / 2) {
        printf("{err}warning: alpha_0 is above the deflated log expenditure of most households;\n")
        printf("{err}         a(p) then exceeds their expenditure (Banks, Blundell and Lewbel set\n")
        printf("{err}         alpha_0 just below the smallest ln x; equaids's default does)\n")
        nw++
    }
    st_local("nl0neg", strofreal(nl))

    // ---- D3: demographics ----
    DS = J(K, 9, .)
    if (K > 0) {
        printf("\n{txt}D3. Demographic variables (Ray scaling m0 = 1 + rho'z)\n{hline 78}\n")
        printf("{txt}%-16s %8s %9s %9s %9s %9s %8s %8s\n", "variable", "distinct", "min", "max",
               "mean", "s.d.", "p rare", "Np(1-p)")
        for (k = 1; k <= K; k++) {
            z  = Z[., k]
            nd = rows(uniqrows(z))
            DS[k, 1..5] = (nd, min(z), max(z), sum(w :* z), sqrt(_eqd_wcov(z, w)))
            if (nd == 2) {
                lo = min(z)
                pr = sum(w :* (z :== lo))
                pr = min((pr, 1 - pr))
                DS[k, 6..7] = (pr, N * pr * (1 - pr))
            }
            printf("{txt}%-16s {res}%8.0f %9.3g %9.3g %9.3g %9.3g %8.3f %8.1f\n", abbrev(zn[k], 16),
                   DS[k, 1], DS[k, 2], DS[k, 3], DS[k, 4], DS[k, 5], DS[k, 6], DS[k, 7])
            if (nd == 2 & DS[k, 7] < 30) {
                printf("{err}warning: %s: rare modality, effective size N p(1-p) = %4.1f\n", zn[k], DS[k, 7])
                nw++
            }
            if (nd > 2 & nd <= 12 & all(z :== round(z)) & min(z) != 0) {
                printf("{txt}note: %s has %g integer values: if categorical, use indicators\n", zn[k], nd)
            }
            if (min(z) < 0) {
                printf("{err}warning: %s takes negative values: m0 = 1 + rho'z can approach 0 or turn\n", zn[k])
                printf("{err}         negative; counts and indicators (z >= 0) keep m0 away from 0\n")
                nw++
            }
        }
        // correlation of each demographic with the relative log prices: a
        // demographic that moves with prices (rural households larger and
        // facing lower prices) biases the price elasticities when omitted,
        // and when included lets gamma and eta trade off -- the instability
        // that D7 (stability) measures.  Threshold 0.3 provisional (to be
        // calibrated on constructed and real cases)
        Rp = LP[|1, 1 \ N, M - 1|] :- LP[., M]
        for (k = 1; k <= K; k++) {
            cz = J(1, M - 1, .)
            for (j = 1; j < M; j++) {
                Cz = _eqd_wcov((Z[., k], Rp[., j]), w)
                cz[j] = (Cz[1, 1] > 0 & Cz[2, 2] > 0 ? Cz[1, 2] / sqrt(Cz[1, 1] * Cz[2, 2]) : 0)
            }
            DS[k, 8] = max(abs(cz))
            // and with ln x: omitted, a demographic correlated with
            // expenditure biases the expenditure elasticities (household
            // size on the data of Lecocq and Robin: 0.30, D7 z up to 21)
            Cz = _eqd_wcov((Z[., k], lx), w)
            DS[k, 9] = (Cz[1, 1] > 0 & Cz[2, 2] > 0 ? Cz[1, 2] / sqrt(Cz[1, 1] * Cz[2, 2]) : 0)
            printf("{txt}%-16s largest |corr.| with a relative log price {res}%6.3f{txt}; corr. with ln x {res}%6.3f\n",
                   abbrev(zn[k], 16), DS[k, 8], DS[k, 9])
            if (DS[k, 8] > .3) {
                printf("{txt}note: %s moves with prices: omitting it biases the price elasticities;\n", zn[k])
                printf("{txt}      including it can let gamma and eta trade off (option -stability-)\n")
            }
            if (abs(DS[k, 9]) > .3) {
                printf("{txt}note: %s moves with expenditure: omitting it biases the expenditure\n", zn[k])
                printf("{txt}      elasticities (option -stability-)\n")
            }
        }
        st_matrix(dsn, DS)
        st_matrixrowstripe(dsn, (J(K, 1, ""), zn'))
        st_matrixcolstripe(dsn, (J(9, 1, ""), ("distinct" \ "min" \ "max" \ "mean" \ "sd" \ "p_rare" \ "Neff" \ "corr_p" \ "corr_x")))
    }

    // ---- D4: conditioning of the design at the starting point ----
    // columns of the Jacobian at the start (b = c = 1): constant, relative
    // log prices, l, l^2 (QUAIDS), z, z*l (eta); Belsley, Kuh and Welsch:
    // uncentered columns with the constant, scaled to unit length
    X  = J(N, 1, 1), (LP[|1, 1 \ N, M - 1|] :- LP[., M]), l0
    xn = "constant", ("ln(p" :+ strofreal(1..M - 1) :+ "/p" :+ strofreal(M) :+ ")"), "l"
    if (qd) {
        X  = X, l0 :^ 2
        xn = xn, "l^2"
    }
    if (K > 0) {
        X  = X, Z, Z :* l0
        xn = xn, zn, (zn :+ "*l")
    }
    printf("\n{txt}D4. Conditioning of the regressors at the starting point (Belsley, Kuh, Welsch)\n{hline 78}\n")
    B = _eqd_bkw(X :* sqrt(w), eta)
    etaall = eta
    cmax = max(eta)
    // two levels (Belsley, Kuh and Welsch: 30 to 100 moderate to strong,
    // above 100 strong): calibrated on three data sets (audit/diag/), where
    // the default alpha_0 gives 23 to 36 with standard errors validated by
    // the bootstrap, and the ill-conditioned estimates are all above 100
    printf("{txt}largest condition index: {res}%9.1f{txt}   (30-100: moderate, above 100: strong)\n", cmax)
    for (k = 1; k <= rows(eta); k++) {
        if (eta[k] > 30) {
            lst = ""
            for (j = 1; j <= cols(X); j++) if (B[k, j] > 0.5) lst = lst + " " + xn[j]
            if (lst != "" & eta[k] > 100) {
                printf("{err}warning: strong near dependency (index %6.1f) among:%s\n", eta[k], lst)
                nw++
            }
            else if (lst != "") printf("{txt}note: moderate near dependency (index %6.1f) among:%s\n", eta[k], lst)
        }
    }
    cq = .
    if (qd) {
        // the quadratic block alone: when alpha_0 lies far from the data, l
        // varies little relative to its level and l^2 is almost linear in
        // l over the sample: lambda is weakly identified
        (void) _eqd_bkw((J(N, 1, 1), l0, l0 :^ 2) :* sqrt(w), etaq)
        cq = max(etaq)
        printf("{txt}condition index of (1, l, l^2) alone: {res}%9.1f\n", cq)
        if (cq > 100) {
            printf("{err}warning: l^2 is almost a linear function of l over the sample: the\n")
            printf("{err}         quadratic coefficients lambda are weakly identified and their\n")
            printf("{err}         sandwich standard errors unreliable.  alpha_0 far from the\n")
            printf("{err}         log expenditures causes it (see D2, and D6 with -sensitivity-)\n")
            nw++
        }
        else if (cq > 30) {
            printf("{txt}note: (1, l, l^2) moderately conditioned; the default alpha_0 gives 23 to 36\n")
            printf("{txt}      on the three data sets used to validate equaids, where the standard\n")
            printf("{txt}      errors hold\n")
        }
    }
    // r(bkw): condition index, then the variance-decomposition proportions
    st_matrix(bkwn, (etaall, B))
    st_matrixcolstripe(bkwn, (J(cols(X) + 1, 1, ""), ("cond_index", xn)'))
    st_local("condmax", strofreal(cmax, "%21.0g"))
    st_local("condquad", strofreal(cq, "%21.0g"))

    // ---- D5: small goods ----
    printf("\n{txt}D5. Small goods\n{hline 78}\n")
    lst = ""
    for (j = 1; j <= M; j++) if (S[j] < 0.01) lst = lst + sprintf(" %s (%5.3f%%)", wn[j], 100 * S[j])
    if (lst != "") {
        printf("{err}warning: goods with less than 1%% of total expenditure:%s.\n", lst)
        printf("{err}         QUAIDS can predict their shares near zero or negative, and the mean\n")
        printf("{err}         of the household elasticities divides by these shares; the aggregate\n")
        printf("{err}         elasticities (equaids's default) do not.  Grouping goods also helps.\n")
        nw++
    }
    else printf("{txt}every good has at least 1%% of total expenditure\n")
    st_local("nwarn", strofreal(nw))
}

// condition indexes (column eta) and variance-decomposition proportions
// (rows = components, columns = regressors) of X, columns scaled to unit length
real matrix _eqd_bkw(real matrix X, real colvector eta)
{
    real matrix    Xs, U, Vt, V, Phi
    real colvector s
    Xs = X :/ sqrt(colsum(X :^ 2))
    svd(Xs, U, s, Vt)
    V   = Vt'
    eta = max(s) :/ s
    Phi = (V :^ 2) :/ (s' :^ 2)
    return((Phi :/ rowsum(Phi))')
}

// D7, standard error of the exact change: the difference of the influence
// functions of the full (_d7f*) and restricted (_d7r*) estimates, household
// by household on the same sample, summed robustly
void _eqd_d7se(string scalar tu, real scalar M, string scalar sx, string scalar su)
{
    real matrix    D
    real scalar    N
    string rowvector nx, nu
    nx = "_d7fx" :+ strofreal(1..M)
    nu = "_d7fu" :+ strofreal(1..M * M)
    D  = st_data(., nx, tu) - st_data(., subinstr(nx, "_d7f", "_d7r"), tu)
    N  = rows(D)
    st_matrix(sx, sqrt(N / (N - 1) :* colsum(D :^ 2)))
    D  = st_data(., nu, tu) - st_data(., subinstr(nu, "_d7f", "_d7r"), tu)
    st_matrix(su, rowshape(sqrt(N / (N - 1) :* colsum(D :^ 2)), M))
}

// weighted covariance matrix, w summing to one (Mata's variance(X, w) reads
// w as frequency weights and divides by sum(w) - 1, zero here)
real matrix _eqd_wcov(real matrix X, real colvector w)
{
    real matrix Xc
    Xc = X :- colsum(w :* X)
    return(cross(Xc, w, Xc))
}

real scalar _eqd_median(real colvector x)
{
    real colvector s
    real scalar    n
    s = sort(x, 1)
    n = rows(s)
    return(mod(n, 2) ? s[(n + 1) / 2] : (s[n / 2] + s[n / 2 + 1]) / 2)
}
end

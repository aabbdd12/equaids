*! equaids 1.1.0  2026-09-27  Abdelkrim Araar
*! AIDS and QUAIDS demand systems: iterated FGNLS by Gauss-Newton with an
*! analytic Jacobian, in Mata.
*!
*! Estimation (FGNLS), demographics by Ray scaling, weights; aggregate
*! elasticities (default), at the means and household mean; analytic
*! variances from influence functions: robust, cluster and survey design
*! (svyset), or the bootstrap of the whole procedure; missing prices filled
*! by group means (pimpute); correction for the non-buyers (Shonkwiler and
*! Yen 1999); diagnostics; estat diagnostics and estat engel.
*!
*! Model (Banks, Blundell and Lewbel 1997; Poi 2012 parameterization):
*!   w_i = alpha_i + sum_j gamma_ij ln p_j + (beta_i + eta_i'z) l
*!         + lambda_i l^2 / (b(p) c(p,z))
*!   l   = ln x - ln m0(z) - ln a(p),  m0 = 1 + rho'z
*!   ln a = alpha_0 + sum alpha_k ln p_k + 1/2 sum sum gamma_kl ln p_k ln p_l
*!   b = prod p_k^beta_k,  c = prod p_k^(eta_k'z)
*! Adding-up, homogeneity and symmetry hold by construction (Poi's free
*! parameters); alpha_0 is fixed: by default just below the smallest ln x
*! (min ln x - 0.1), or imposed by anot(#).
program define equaids, eclass
    version 14.2
    if replay() {
        if "`e(cmd)'" != "equaids" error 301
        * internal entry for estat diagnostics: the subprograms of this file
        * are not visible from equaids_estat.ado
        if trim(`"`0'"') == ", _warnings" {
            _equaids_warnings
            exit
        }
        * internal entry for estat engel: the Mata of this file, loaded
        * with it (also after estimates use)
        if substr(trim(`"`0'"'), 1, 8) == ", _engel" {
            _equaids_engel , `=substr(trim(`"`0'"'), 9, .)'
            exit
        }
        _equaids_display `0'
        exit
    }
    * vce(bootstrap): the whole procedure (missing prices, probits, system)
    * estimated again on each replication
    if regexm(`"`0'"', "vce\( *(boot|bootstrap|bs) *[,)]") {
        _equaids_boot `0'
        exit
    }
    _equaids_estimate `0'
end

program define _equaids_estimate, eclass sortpreserve
    version 14.2
    syntax varlist(min=3 numeric) [if] [in] [aweight fweight pweight iweight], ///
        [ ANOT(string) PRices(varlist numeric) LNPRices(varlist numeric)                    ///
          EXPenditure(varname numeric) LNEXPenditure(varname numeric)          ///
          DEMOgraphics(varlist numeric) noQUadratic                            ///
          METhod(string) VCE(string) Level(cilevel)                            ///
          FROM(name) TOLerance(real 1e-6) ITERate(integer 300) noLOg          ///
          SNames(string) DEC(integer 4) DISLAS(integer 1) DREGRES(integer 0)   ///
          ELASticities(string) noELASTse COMPENSated CHECKS DETail STars       ///
          SAVEres(string) noTABle SAVEIF(name) HHSize(varname numeric) PIMPute(varlist)             ///
          SELection SELGoods(namelist) SELVars(string) ]

    timer clear 100
    timer on 100
    local shares `varlist'
    local M : word count `shares'
    local doe = ("`elastse'" == "")
    _equaids_convparse `elasticities'
    local elasticities `s(conv)'
    * elasticities of the individual: each household counts for its weight
    * times its size, hhsize(), in the whole estimation
    local indiv = ("`elasticities'" == "individuals")
    if `indiv' & "`hhsize'" == "" {
        di as err "elasticities(individuals) requires hhsize(varname), the size of the household"
        exit 198
    }
    * ---- selection of the buyers (Shonkwiler and Yen 1999): the goods
    * corrected (selgoods(), default: those with zero shares, the last good
    * excepted) and the variables of the probits only (selvars(): segments
    * separated by ";", "vars" for every corrected good or "good: vars") ----
    if "`selgoods'" != "" | `"`selvars'"' != "" local selection "selection"
    local sel = ("`selection'" != "")
    local zall ""
    local lastgood : word `M' of `shares'
    if `sel' {
        foreach g of local selgoods {
            if !`: list g in shares' {
                di as err "selgoods(): `g' is not one of the shares (`shares')"
                exit 198
            }
            if "`g'" == "`lastgood'" {
                di as err "selgoods(): `g' is the last good, which closes the system and is not corrected"
                exit 198
            }
        }
        local selgoods : list uniq selgoods
        local zcom ""
        foreach g of local shares {
            local z_`g' ""
        }
        local rest `"`selvars'"'
        while `"`rest'"' != "" {
            gettoken seg rest : rest, parse(";")
            local seg = trim(`"`seg'"')
            if `"`seg'"' == ";" | `"`seg'"' == "" continue
            if strpos(`"`seg'"', ":") {
                local gg = trim(substr(`"`seg'"', 1, strpos(`"`seg'"', ":") - 1))
                local vv = trim(substr(`"`seg'"', strpos(`"`seg'"', ":") + 1, .))
                if `"`gg'"' == "" {
                    di as err `"selvars(): a good name is missing before ":" in "`seg'""'
                    exit 198
                }
                if `: word count `gg'' > 1 {
                    di as err `"selvars(): one good per segment, "`gg'" is not one good; write "good1: vars ; good2: vars""'
                    exit 198
                }
                if strpos(`"`vv'"', ":") {
                    di as err `"selvars(): one colon per segment in "`seg'"; separate the goods by ";""'
                    exit 198
                }
                if `"`vv'"' == "" {
                    di as err "selvars(): no variable after `gg':"
                    exit 198
                }
                if !`: list gg in shares' {
                    di as err "selvars(): `gg' is not one of the shares (`shares')"
                    exit 198
                }
                if "`selgoods'" != "" & !`: list gg in selgoods' {
                    di as err "selvars(): `gg' is not corrected; add it to selgoods() (`selgoods')"
                    exit 198
                }
                capture unab vv : `vv'
                if _rc {
                    di as err "selvars(), `gg': variable not found in `vv'"
                    exit 111
                }
                confirm numeric variable `vv'
                local z_`gg' `z_`gg'' `vv'
            }
            else {
                capture unab vv : `seg'
                if _rc {
                    di as err `"selvars(): variable not found in "`seg'""'
                    exit 111
                }
                confirm numeric variable `vv'
                local zcom `zcom' `vv'
            }
        }
        * a probit-only variable is an exclusion restriction: not in the model
        local xin `shares' `prices' `lnprices' `expenditure' `lnexpenditure' `demographics'
        foreach g of local shares {
            local z_`g' : list zcom | z_`g'
            local zall : list zall | z_`g'
        }
        local bad : list zall & xin
        if "`bad'" != "" {
            di as err "selvars(): `bad' already in the model (share, price, expenditure or demographic);"
            di as err "a variable of the probit only must be excluded from the share equations"
            exit 198
        }
        if "`vce'" == "conventional" {
            di as err "selection: use vce(robust), vce(cluster) or vce(svy); the conventional variance ignores the probits"
            exit 198
        }
    }
    * -detail- is the reporting options at once
    if "`detail'" != "" {
        local compensated compensated
        local checks      checks
    }
    if "`snames'" != "" {
        if `: word count `snames'' != `M' {
            di as err "snames(): `: word count `snames'' name(s) for `M' goods"
            exit 198
        }
    }
    if `dec' < 0 | `dec' > 8 {
        di as err "dec() must be between 0 and 8"
        exit 198
    }

    * ---- options ----
    if "`method'" == "" local method fgnls
    if "`method'" != "fgnls" {
        di as err "method(`method') not available yet; only method(fgnls)"
        exit 198
    }
    * variance: robust by default (as easi); pweights imply robust.
    * cluster and svy aggregate the same influence functions by cluster, or
    * by PSU within strata with the finite-population correction (step 3)
    gettoken vt vrest : vce
    local vrest = trim(`"`vrest'"')
    local clustvar ""
    if inlist("`vt'", "", "robust", "r") & `"`vrest'"' == "" local vce robust
    else if inlist("`vt'", "conventional", "conv") & `"`vrest'"' == "" local vce conventional
    else if inlist("`vt'", "cluster", "cl") {
        local vce cluster
        local clustvar `vrest'
        if `: word count `clustvar'' != 1 {
            di as err "vce(cluster clustvar): one cluster variable"
            exit 198
        }
        confirm variable `clustvar'
    }
    else if "`vt'" == "svy" & `"`vrest'"' == "" {
        local vce svy
        if "`weight'" != "" {
            di as err "vce(svy): the sampling weight comes from svyset; do not specify weights"
            exit 198
        }
        quietly svyset
        if "`r(settings)'" == ", clear" | "`r(settings)'" == "" {
            di as err "vce(svy): the data are not svyset"
            exit 119
        }
        local svywv  `r(wvar)'
        local svypsu `r(su1)'
        if "`svypsu'" == "_n" local svypsu ""
        local svystr `r(strata1)'
        local svyfpc `r(fpc1)'
        local single `r(singleunit)'
        if "`single'" == "" local single missing
        if "`svywv'" != "" {
            if "`r(wtype)'" != "pweight" {
                di as err "vce(svy): the svyset weight must be a pweight"
                exit 198
            }
            local weight pweight
            local exp "= `svywv'"
        }
    }
    else {
        di as err "vce(`vce') not allowed; use robust, cluster clustvar, svy or conventional"
        exit 198
    }
    if "`weight'" == "pweight" & "`vce'" == "conventional" {
        di as err "pweights require vce(robust)"
        exit 198
    }
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

    * ---- sample ----
    marksample touse
    * with pimpute(), a missing price does not drop the household: it is filled
    * below from the households of the same group
    if "`pimpute'" == "" markout `touse' `pv'
    markout `touse' `expenditure' `lnexpenditure' `demographics'
    if `sel' & "`zall'" != "" {
        quietly count if `touse'
        local n0 = r(N)
        markout `touse' `zall'
        quietly count if `touse'
        if r(N) < `n0' di as txt "(selvars(): `=string(`n0' - r(N), "%12.0fc")' households dropped for missing values of `zall')"
    }
    if `indiv' {
        markout `touse' `hhsize'
        quietly count if `touse' & `hhsize' <= 0
        if r(N) {
            di as err "hhsize(): `hhsize' has nonpositive values"
            exit 411
        }
    }
    if "`clustvar'" != "" markout `touse' `clustvar', strok
    if "`vce'" == "svy" markout `touse' `svypsu' `svystr' `svyfpc', strok
    if "`weight'" != "" {
        tempvar wt
        quietly gen double `wt' `exp' if `touse'
        quietly replace `touse' = 0 if `wt' <= 0 | missing(`wt')
    }
    quietly count if `touse'
    if r(N) == 0 error 2000
    local N = r(N)

    * ---- logs of prices and expenditure ----
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
        else {
            * a copy when the prices are to be filled
            if "`pimpute'" != "" {
                tempvar l`v'
                quietly gen double `l`v'' = `v' if `touse'
                local lnp `lnp' `l`v''
            }
            else local lnp `lnp' `v'
        }
    }
    * ---- missing prices: the weighted mean of the log prices of the same
    * group (pimpute(); see _equaids_pimpute.ado); the weight is the weight of
    * the estimation (times hhsize() for individuals), as estat rebuilds it
    local pinotes ""
    if "`pimpute'" != "" {
        tempvar pwv
        if "`weight'" != "" quietly gen double `pwv' = `wt' if `touse'
        else quietly gen double `pwv' = 1 if `touse'
        if `indiv' quietly replace `pwv' = `pwv' * `hhsize' if `touse'
        _equaids_pimpute `lnp', touse(`touse') wt(`pwv') groups(`pimpute') names(`shares')
        local pinotes `"`r(notes)'"'
        foreach l of local pinotes {
            di as txt `"`l'"'
        }
        quietly count if `touse'
        if r(N) == 0 error 2000
        local N = r(N)
    }
    * ---- selection: the goods corrected, and good by good the probit-only
    * variables (flags for the engine) ----
    local selgstr ""
    local zmstr ""
    if `sel' {
        quietly count if `touse' & `lastgood' <= 0
        if r(N) {
            di as err "selection: the last good (`lastgood') has zero shares; place last a good that every household buys"
            exit 198
        }
        if "`selgoods'" == "" {
            forvalues k = 1/`=`M'-1' {
                local g : word `k' of `shares'
                quietly count if `touse' & `g' <= 0
                if r(N) local selgoods `selgoods' `g'
            }
            if "`selgoods'" == "" {
                di as err "selection: no good has zero shares; nothing to correct"
                exit 459
            }
        }
        else {
            local sg ""
            foreach g of local selgoods {
                quietly count if `touse' & `g' <= 0
                if r(N) == 0 di as txt "(selection: every household buys `g'; it is left uncorrected)"
                else local sg `sg' `g'
            }
            if "`sg'" == "" {
                di as err "selection: every household buys each of the goods to correct; nothing to correct"
                exit 459
            }
            local selgoods `sg'
        }
        * in the order of the shares, as the parameters delta
        local sgo ""
        forvalues k = 1/`=`M'-1' {
            local g : word `k' of `shares'
            local on : list g in selgoods
            local selgstr `selgstr' `on'
            if `on' local sgo `sgo' `g'
            foreach v of local zall {
                local inz : list v in z_`g'
                local zmstr `zmstr' `=`on' & `inz''
            }
        }
        local selgoods `sgo'
    }
    if "`expenditure'" != "" {
        quietly count if `touse' & `expenditure' <= 0
        if r(N) {
            di as err "`expenditure' has nonpositive values"
            exit 411
        }
        tempvar lnx
        quietly gen double `lnx' = ln(`expenditure') if `touse'
    }
    else local lnx `lnexpenditure'
    quietly summarize `lnx' if `touse', meanonly
    local lnxmin = r(min)
    * alpha_0: by default just below the smallest log expenditure of the
    * estimation sample (Banks, Blundell and Lewbel 1997), so that
    * ln x - ln a(p) stays positive; anot(#) imposes a value (as Poi)
    local a0margin 0.1
    if `"`anot'"' == "" {
        local anot = `lnxmin' - `a0margin'
        local a0rule "min ln x - `a0margin'"
    }
    else {
        capture confirm number `anot'
        if _rc {
            di as err "anot() must be a number"
            exit 198
        }
        local a0rule "user"
    }

    * ---- shares sum to one ----
    tempvar ssum
    quietly egen double `ssum' = rowtotal(`shares') if `touse'
    quietly count if `touse' & abs(`ssum' - 1) > 1e-4
    if r(N) {
        di as err "the shares do not sum to one in " r(N) " observations"
        exit 459
    }

    * ---- weights: aw and pw normalized to sum N, fw and iw as given ----
    * (individuals: the weight times hhsize(), normalized as an aweight when
    * no weight is given)
    tempvar om
    if "`weight'" == "" quietly gen double `om' = 1 if `touse'
    else quietly gen double `om' = `wt' if `touse'
    if `indiv' quietly replace `om' = `om' * `hhsize' if `touse'
    if inlist("`weight'", "aweight", "pweight") | ("`weight'" == "" & `indiv') {
        quietly summarize `om' if `touse', meanonly
        quietly replace `om' = `om' * `N' / r(sum) if `touse'
    }

    * ---- design for the variance: stratum, PSU (or cluster), fpc ----
    tempvar dst dsu dfp
    local vmode = cond("`vce'" == "cluster", 1, cond("`vce'" == "svy", 2, 0))
    if `vmode' {
        if "`svystr'" != "" quietly egen long `dst' = group(`svystr') if `touse'
        else quietly gen long `dst' = 1 if `touse'
        if "`clustvar'`svypsu'" != "" quietly egen long `dsu' = group(`clustvar'`svypsu') if `touse'
        else quietly gen long `dsu' = _n if `touse'
        if "`svyfpc'" != "" quietly gen double `dfp' = `svyfpc' if `touse'
        else quietly gen double `dfp' = . if `touse'
        mata: _eq_desinfo("`dst' `dsu' `dfp'", "`touse'")
        * one PSU (or cluster) in all: the variance is not defined
        if `N_psu' < 2 {
            local u = cond(`vmode' == 1, "cluster", "PSU")
            di as err "vce(`vce'): a single `u' in the estimation sample; the variance is not defined"
            exit 460
        }
        * strata with a single PSU under singleunit(missing): as Stata, the
        * model is estimated and the standard errors are missing (see the note)
    }
    local desv = cond(`vmode', "`dst' `dsu' `dfp'", "")

    * ---- demographics: variability, modalities, proportions ----
    local K : word count `demographics'
    tempname DS VIF
    if `K' > 0 {
        _equaids_demostats `demographics' if `touse', om(`om') mat(`DS')
        local demonotes `"`r(notes)'"'
    }
    * ---- regressors: near collinearity; zero shares; extreme prices ----
    local rnames ""
    forvalues k = 1/`=`M'-1' {
        local rnames `rnames' ln(p`k'/p`M')
    }
    local rnames `rnames' ln_x `demographics'
    mata: _eq_datadiag("`shares'", "`lnp'", "`lnx'", "`demographics'", "`om'", ///
        "`touse'", "`VIF'", "`selgstr'")
    matrix colnames `VIF' = `rnames'
    local demonotes `"`demonotes'`datanotes'"'

    * ---- estimate ----
    tempname bf Vf b V Sig EL
    mata: _eq_fit("`shares'", "`lnp'", "`lnx'", "`demographics'", "`om'", ///
        "`touse'", `anot', `quad', `tolerance', `iterate', "`log'" == "",   ///
        "`from'", "`bf'", "`Vf'", "`b'", "`V'", "`Sig'", "`EL'", `doe', ///
        "`desv'", `vmode', "`single'", `sel', "`zall'", "`selgstr'", "`zmstr'", ///
        "`EL'sdg", "`EL'sal")
    mata: _eq_checks("`EL'")

    * ---- names: Poi's stripe on the full vector ----
    local stripe ""
    forvalues i = 1/`M' {
        local stripe `stripe' alpha:alpha_`i'
    }
    forvalues i = 1/`M' {
        local stripe `stripe' beta:beta_`i'
    }
    forvalues j = 1/`M' {
        forvalues i = `j'/`M' {
            local stripe `stripe' gamma:gamma_`i'_`j'
        }
    }
    if `quad' {
        forvalues i = 1/`M' {
            local stripe `stripe' lambda:lambda_`i'
        }
    }
    foreach z of local demographics {
        forvalues i = 1/`M' {
            local stripe `stripe' eta:eta_`z'_`i'
        }
    }
    foreach z of local demographics {
        local stripe `stripe' rho:rho_`z'
    }
    if `sel' {
        foreach g of local selgoods {
            local stripe `stripe' delta:delta_`g'
        }
    }
    * robust, cluster or design variance of the coefficients from the
    * influence functions
    * the model-based variance A^-1 (inverse information) of the free
    * parameters, kept whatever the vce: it is the curvature of the FGNLS
    * objective, used by equaidsdiag's first-order change of the estimates
    tempname VM
    matrix `VM' = `Vf'
    if "`vce'" != "conventional" {
        matrix `Vf' = `EL'Vt
        mata: st_matrix("`V'", st_matrix("`EL'D") * st_matrix("`Vf'") * st_matrix("`EL'D")')
    }
    matrix colnames `b' = `stripe'
    matrix colnames `V' = `stripe'
    matrix rownames `V' = `stripe'

    * missing variance (single-PSU strata, singleunit(missing)): posted
    * without V, as the standard errors are not defined
    mata: st_local("vmiss", strofreal(hasmissing(st_matrix("`V'"))))
    local pV = cond(`vmiss', "", "`V'")
    * the design df is computed only under svy: a `=...' macro is expanded
    * before the if is tested, and N_psu is empty otherwise
    if `vmode' == 2 {
        local ddf = `N_psu' - `N_strata'
        ereturn post `b' `pV', esample(`touse') obs(`N') dof(`ddf')
    }
    else ereturn post `b' `pV', esample(`touse') obs(`N')
    ereturn scalar V_missing = `vmiss'
    if `vmode' == 1 {
        ereturn scalar N_clust = `N_psu'
        ereturn local  clustvar "`clustvar'"
    }
    if `vmode' == 2 {
        ereturn scalar N_strata = `N_strata'
        ereturn scalar N_psu    = `N_psu'
        ereturn scalar N_single = `N_single'
        ereturn scalar df_r     = `N_psu' - `N_strata'
        ereturn local  singleunit "`single'"
        ereturn local  su1     "`svypsu'"
        ereturn local  strata1 "`svystr'"
        ereturn local  fpc1    "`svyfpc'"
    }
    ereturn matrix b_free = `bf'
    ereturn matrix V_free = `Vf'
    ereturn matrix V_model = `VM'
    ereturn matrix Sigma  = `Sig'
    ereturn scalar ll        = `ll'
    ereturn scalar N_w       = `Nw'
    ereturn scalar anot      = `anot'
    ereturn local  anot_rule "`a0rule'"
    ereturn scalar ngoods    = `M'
    ereturn scalar ndemos    = `K'
    ereturn scalar iter      = `iter'
    ereturn scalar converged = `conv'
    ereturn scalar nrgrad    = `nd'
    ereturn scalar sigdif    = `sd'
    ereturn scalar stalled   = `stall'
    ereturn scalar n_boundary = `nhb'
    ereturn scalar rcond     = `rc'
    ereturn scalar m0_min    = `m0min'
    ereturn scalar n_m0low   = `nm0low'
    ereturn scalar lnx_min   = `lnxmin'
    ereturn scalar n_lneg    = `nlneg'
    ereturn scalar n_shout   = `nneg'
    ereturn scalar cond_x    = `condx'
    ereturn matrix vif       = `VIF'
    if `K' > 0 ereturn matrix demo_stats = `DS'
    ereturn local  data_notes `"`demonotes'"'
    * elasticities: aggregate (default), at the means, household mean
    foreach s in x u c xm um cm xh uh ch S ns {
        matrix colnames `EL'`s' = `shares'
        if inlist("`s'", "u", "c", "um", "cm", "uh", "ch") matrix rownames `EL'`s' = `shares'
    }
    ereturn matrix elas_x   = `EL'x
    ereturn matrix elas_u   = `EL'u
    ereturn matrix elas_c   = `EL'c
    ereturn matrix elas_xm  = `EL'xm
    ereturn matrix elas_um  = `EL'um
    ereturn matrix elas_cm  = `EL'cm
    ereturn matrix elas_xh  = `EL'xh
    ereturn matrix elas_uh  = `EL'uh
    ereturn matrix elas_ch  = `EL'ch
    ereturn matrix aggshare = `EL'S
    matrix colnames `EL'Sm = `shares'
    matrix colnames `EL'Sh = `shares'
    ereturn matrix shares_m = `EL'Sm
    ereturn matrix shares_h = `EL'Sh
    ereturn matrix n_fsmall = `EL'ns
    * variances of the aggregate elasticities (robust, analytic)
    if `doe' {
        local un ""
        foreach i of local shares {
            foreach j of local shares {
                local un `un' `i':`j'
            }
        }
        matrix rownames `EL'Vx = `shares'
        matrix colnames `EL'Vx = `shares'
        matrix rownames `EL'Vu = `un'
        matrix colnames `EL'Vu = `un'
        matrix rownames `EL'Vc = `un'
        matrix colnames `EL'Vc = `un'
        ereturn matrix V_elas_x = `EL'Vx
        ereturn matrix V_elas_u = `EL'Vu
        ereturn matrix V_elas_c = `EL'Vc
        ereturn matrix J_elas_x = `EL'Jx
        ereturn matrix J_elas_u = `EL'Ju
        ereturn matrix J_elas_c = `EL'Jc
        * at the means (households, individuals) and the household mean
        foreach f in m h {
            matrix rownames `EL'Vx`f' = `shares'
            matrix colnames `EL'Vx`f' = `shares'
            matrix rownames `EL'Vu`f' = `un'
            matrix colnames `EL'Vu`f' = `un'
            matrix rownames `EL'Vc`f' = `un'
            matrix colnames `EL'Vc`f' = `un'
            ereturn matrix V_elas_x`f' = `EL'Vx`f'
            ereturn matrix V_elas_u`f' = `EL'Vu`f'
            ereturn matrix V_elas_c`f' = `EL'Vc`f'
        }
    }
    ereturn scalar elastse     = `doe'
    ereturn scalar chk_engel   = `chk_engel'
    ereturn scalar chk_cournot = `chk_cournot'
    ereturn scalar chk_homog   = `chk_homog'
    ereturn scalar chk_symm    = `chk_symm'
    * display settings, used again on replay
    ereturn scalar dec     = `dec'
    ereturn scalar dislas  = `dislas'
    ereturn scalar dregres = `dregres'
    ereturn local  elasticities "`elasticities'"
    ereturn local  report = trim("`compensated' `checks'")
    if "`snames'" != "" ereturn local snames "`snames'"
    ereturn local  lhs          "`shares'"
    if "`prices'" != "" ereturn local prices "`prices'"
    else                ereturn local lnprices "`lnprices'"
    if "`expenditure'" != "" ereturn local expenditure "`expenditure'"
    else                     ereturn local lnexpenditure "`lnexpenditure'"
    ereturn local  demographics "`demographics'"
    ereturn local  model  = cond(`quad', "QUAIDS", "AIDS")
    ereturn local  method "fgnls"
    ereturn local  vce    "`vce'"
    * the weight used, so that estat rebuilds the same one
    if `indiv' {
        if "`weight'" != "" {
            ereturn local wtype "`weight'"
            ereturn local wexp  "= (`=trim(subinstr("`exp'", "=", "", 1))') * (`hhsize')"
        }
        else {
            ereturn local wtype "aweight"
            ereturn local wexp  "= `hhsize'"
        }
    }
    else if "`weight'" != "" {
        ereturn local wtype "`weight'"
        ereturn local wexp  "`exp'"
    }
    if "`hhsize'" != "" ereturn local hhsize "`hhsize'"
    if "`pimpute'" != "" ereturn local pimpute "`pimpute'"
    if `sel' {
        ereturn local selection "shonkwiler-yen"
        ereturn local selgoods "`selgoods'"
        ereturn local selvars `"`selvars'"'
        foreach g of local selgoods {
            ereturn local sel_z_`g' "`z_`g''"
        }
        local g1 ""
        forvalues k = 1/`=`M'-1' {
            local g1 `g1' `: word `k' of `shares''
        }
        matrix colnames `EL'dl  = `g1'
        matrix colnames `EL'sdl = `g1'
        ereturn matrix sel_delta    = `EL'dl
        ereturn matrix se_sel_delta = `EL'sdl
        matrix colnames `EL'sal = `g1'
        matrix rownames `EL'sal = _cons `=cond("`prices'" != "", "`prices'", "`lnprices'")' `demographics' `zall'
        ereturn matrix sel_alpha = `EL'sal
        matrix rownames `EL'sdg = `g1'
        matrix colnames `EL'sdg = Buy_pct Pseudo_R2 Perfect VIF_beta
        ereturn matrix sel_diag = `EL'sdg
        * (theta, alpha of each corrected good) and its variance: estat engel
        ereturn local sel_vars "`zall'"
        ereturn matrix sel_psi   = `EL'psi
        ereturn matrix V_sel_psi = `EL'Vpsi
    }
    * distance of the estimate to the boundary of Ray's scaling, in standard
    * errors: the smallest m0(z) / se(m0(z)) over the households; near the
    * boundary the estimator is not regular and the linearized standard
    * errors are unreliable
    if `K' > 0 & !`vmiss' {
        tempvar esm
        quietly gen byte `esm' = e(sample)
        mata: _eq_m0dist("`demographics'", "`esm'")
        ereturn scalar m0_t = `m0t'
    }
    timer off 100
    quietly timer list 100
    ereturn scalar time = r(t100)
    ereturn local  estat_cmd "equaids_estat"
    ereturn local  cmd "equaids"
    if "`saveres'" != "" local saveres `"saveres(`saveres')"'
    _equaids_display, level(`level') `stars' `saveres' `table'
end

* vce(bootstrap [, reps(#) seed(#) strata(varname) psu(varname) svy]): the
* point estimates of the full sample, and variances from the replications of
* the whole procedure (missing prices filled, probits and system estimated
* again on each resample of households, or of PSUs within strata); the
* households keep their weights. Replaces the variances of the coefficients,
* of delta, of the four types of elasticities, V_free and V_sel_psi.
program define _equaids_boot, eclass
    version 14.2
    syntax varlist(min=3 numeric) [if] [in] [aweight fweight pweight iweight], ///
        VCE(string) [ Level(passthru) STars SAVEres(string) noTABle noLOG     ///
        noELASTse * ]
    * the command's elements, before syntax parses the options of the bootstrap
    local vl `varlist'
    local ifin `if' `in'
    local wt0 `weight'
    local ex0 `"`exp'"'
    local opts `"`options'"'
    gettoken vt vo : vce, parse(",")
    gettoken comma vo : vo, parse(",")
    local 0 `", `vo'"'
    syntax [, Reps(integer 200) SEED(string) STRata(varname) PSU(varname) SVY]
    local varlist `vl'
    local weight `wt0'
    local exp `"`ex0'"'
    local options `"`opts'"'
    if `reps' < 2 {
        di as err "vce(bootstrap): reps() must be at least 2"
        exit 198
    }
    local wgt = cond("`weight'" != "", "[`weight'`exp']", "")
    if "`svy'" != "" {
        quietly svyset
        if "`r(settings)'" == ", clear" | "`r(settings)'" == "" {
            di as err "vce(bootstrap, svy): the data are not svyset"
            exit 119
        }
        if "`psu'" == "" & !inlist("`r(su1)'", "", "_n") local psu `r(su1)'
        if "`strata'" == "" local strata `r(strata1)'
        if "`r(wvar)'" != "" {
            if "`weight'" != "" {
                di as err "vce(bootstrap, svy): the weight comes from svyset; do not specify weights"
                exit 198
            }
            local wgt "[pweight = `r(wvar)']"
        }
    }
    local bso ""
    if "`strata'" != "" local bso `bso' strata(`strata')
    if "`psu'" != "" local bso `bso' cluster(`psu')
    local design = trim(cond("`psu'" != "", "PSUs `psu'", "households") + ///
        cond("`strata'" != "", " within strata `strata'", ""))
    if "`seed'" != "" set seed `seed'
    * the replications iterate up to 1000 times unless iterate() is given: a
    * draw stopped at the limit just above the tolerance would otherwise be
    * dropped (the tail of Gauss-Newton is linear)
    local itopt ""
    if !regexm(`" `options'"', " iter(a|at|ate)?\(") local itopt "iterate(1000)"

    * the full sample: the point estimates
    di as txt "(estimating the full sample)"
    quietly _equaids_estimate `varlist' `ifin' `wgt', `options' `level' vce(robust) notable nolog
    tempvar es
    quietly gen byte `es' = e(sample)
    local sel = ("`e(selection)'" != "")
    local conv = e(converged)
    tempname full
    _estimates hold `full', copy restore
    if !`conv' {
        di as err "the full sample did not converge; see equaids without vce(bootstrap)"
        exit 430
    }

    * the replications
    di as txt "Bootstrap replications (" as res `reps' as txt ", resampling `design')"
    mata: __eq_bs = J(0, 0, .)
    local ok 0
    forvalues r = 1/`reps' {
        preserve
        quietly keep if `es'
        quietly bsample, `bso'
        capture quietly _equaids_estimate `varlist' `wgt', `options' `itopt' vce(robust) noelastse notable nolog
        if !_rc & e(converged) == 1 {
            mata: __eq_bs = (rows(__eq_bs) ? __eq_bs \ _eq_bvec(`sel') : _eq_bvec(`sel'))
            local ++ok
            di as txt "." _c
        }
        else di as err "x" _c
        restore
        if mod(`r', 50) == 0 di as txt " " `r'
    }
    if mod(`reps', 50) di ""
    _estimates unhold `full'
    if `ok' < 2 {
        di as err "vce(bootstrap): fewer than 2 successful replications"
        exit 430
    }

    * the variances from the replications, under the names of the analytic ones
    tempname EB T
    mata: _eq_bpost(__eq_bs, `sel', "`EB'")
    mata: mata drop __eq_bs
    matrix `T' = e(V)
    mata: st_replacematrix("`T'", st_matrix("`EB'Vb"))
    ereturn repost V = `T'
    matrix `T' = e(V_free)
    mata: st_replacematrix("`T'", st_matrix("`EB'Vf"))
    ereturn matrix V_free = `T'
    foreach m in x u c {
        foreach f in "" m h {
            capture confirm matrix e(V_elas_`m'`f')
            if !_rc {
                matrix `T' = e(V_elas_`m'`f')
                mata: st_replacematrix("`T'", st_matrix("`EB'V`m'`f'"))
                ereturn matrix V_elas_`m'`f' = `T'
            }
        }
    }
    if `sel' {
        matrix `T' = e(se_sel_delta)
        mata: st_replacematrix("`T'", st_matrix("`EB'sdl"))
        ereturn matrix se_sel_delta = `T'
        matrix `T' = e(V_sel_psi)
        mata: st_replacematrix("`T'", st_matrix("`EB'Vpsi"))
        ereturn matrix V_sel_psi = `T'
    }
    ereturn scalar V_missing = 0
    ereturn scalar elastse = 1
    ereturn scalar N_reps = `reps'
    ereturn scalar N_reps_ok = `ok'
    ereturn local vce "bootstrap"
    ereturn local vcetype "Bootstrap"
    ereturn local boot_design "`design'"
    if "`seed'" != "" ereturn local boot_seed "`seed'"
    if "`saveres'" != "" local saveres `"saveres(`saveres')"'
    _equaids_display, `level' `stars' `saveres' `table'
end

* elasticities(): which elasticities are reported (names shared with duvm)
*   market       of total demand: at the aggregate shares (the default)
*   households   of the household: at the means (each household its weight)
*   individuals  of the individual: the same, each household counting for its
*                weight times its size, hhsize(), in the whole estimation
*   hhmean       the mean of the household elasticities
* the names of 1.0.0 are accepted: aggregate, means, household
program define _equaids_convparse, sclass
    local c = lower(trim(`"`0'"'))
    if "`c'" == "" | inlist("`c'", "market", "mkt", "aggregate", "agg") local c market
    else if inlist("`c'", "households", "means", "mean", "atmeans") local c households
    else if inlist("`c'", "individuals", "individual", "ind", "persons") local c individuals
    else if inlist("`c'", "hhmean", "household") local c hhmean
    else {
        di as err "elasticities(`0') not allowed; use market (the default), households, individuals or hhmean"
        exit 198
    }
    sreturn local conv "`c'"
end

program define _equaids_display
    syntax [, Level(cilevel) DEC(integer -1) DISLAS(integer -1) DREGRES(integer -1) ///
        ELASticities(string) COMPENSated CHECKS DETail STars SAVEres(string) noTABle]
    if `dec'     == -1 local dec     = e(dec)
    if `dislas'  == -1 local dislas  = e(dislas)
    if `dregres' == -1 local dregres = e(dregres)
    if "`elasticities'" == "" local elasticities `e(elasticities)'
    _equaids_convparse `elasticities'
    local conv `s(conv)'
    * individuals is an estimation (weights times hhsize()), not a display
    _equaids_convparse `e(elasticities)'
    local est `s(conv)'
    if "`est'" == "individuals" & "`conv'" != "individuals" {
        di as err "the estimation is for individuals (weights x `e(hhsize)'): market, households"
        di as err "and hhmean need the estimation without elasticities(individuals)"
        exit 198
    }
    if "`est'" != "individuals" & "`conv'" == "individuals" {
        di as err "elasticities(individuals) needs the estimation with hhsize() and elasticities(individuals)"
        exit 198
    }
    * what to show: what was asked for at estimation, plus anything asked now
    local rep "`e(report)'"
    if "`detail'" != "" local rep "compensated checks"
    foreach o in compensated checks {
        if "``o''" != "" local rep "`rep' `o'"
    }
    local shcq : list posof "compensated" in rep
    local shck : list posof "checks"      in rep

    di ""
    di as txt "{hline 78}"
    di as txt "`e(model)' demand system, iterated FGNLS" ///
        _col(50) "Number of obs  = " as res %10.0fc e(N)
    di as txt "Goods          = " as res %-10.0g e(ngoods) ///
        _col(50) as txt "Demographics   = " as res %10.0g e(ndemos)
    local a0r = cond("`e(anot_rule)'" == "user", "", " (`e(anot_rule)')")
    di as txt "alpha_0        = " as res %-8.4g e(anot) as txt "`a0r'" ///
        _col(50) as txt "Log likelihood = " as res %10.2f e(ll)
    if "`e(wtype)'" != "" di as txt "Weights        = " as res "`e(wtype)' `e(wexp)'"
    di as txt "Variance       = " as res %-10s "`e(vce)'" ///
        _col(50) as txt "Elasticities   = " as res %10s "`conv'"
    if "`e(vce)'" == "bootstrap" di as txt "Replications   = " as res %-10.0fc e(N_reps_ok) ///
        as txt " of " as res e(N_reps) as txt cond("`e(boot_design)'" != "", " (`e(boot_design)')", "")
    if "`e(vce)'" == "cluster" di as txt "Clusters       = " as res %-10.0fc e(N_clust) ///
        as txt " (`e(clustvar)')"
    if "`e(vce)'" == "svy" {
        di as txt "Strata         = " as res %-10.0fc e(N_strata) ///
            _col(50) as txt "PSUs           = " as res %10.0fc e(N_psu)
        di as txt "Design df      = " as res %-10.0fc e(df_r) ///
            _col(50) as txt "Single-PSU str.= " as res %10.0fc e(N_single)
    }
    if e(converged) di as txt "Converged in " as res e(iter) as txt " iteration(s)" ///
        _col(50) as txt "Execution time = " as res %8.2f e(time) as txt " s"
    else di as err "Not converged after " e(iter) " iteration(s) (see below)"
    di as txt "{hline 78}"
    if "`e(selection)'" != "" & "`table'" == "" _equaids_seltable

    if `dregres' {
        di ""
        di as txt "Table 00: Estimated coefficients (Poi's parameterization)"
        ereturn display, level(`level')
    }
    if !inlist(`dislas', 0, 1) {
        di as err "dislas() must be 0 or 1"
        exit 198
    }
    if "`table'" == "" | "`saveres'" != "" {
        local fq ""
        if "`saveres'" != "" local fq `""`saveres'""'
        _equaids_tables `fq', dec(`dec') dislas(`dislas') conv(`conv') ///
            shcq(`shcq') `stars' screen(`="`table'" == ""')
    }
    if "`table'" == "" {
        _equaids_elasnotes, conv(`conv') dislas(`dislas') shcq(`shcq') shck(`shck')
        if "`saveres'" != "" di as txt "Tables written to " as res `"`saveres'"'
    }
    _equaids_warnings
end

* The selection of the buyers: one probit per corrected good, and the
* coefficient delta of its density term in the expected share
program define _equaids_seltable
    tempname D G
    matrix `D' = e(sel_delta)
    matrix `G' = e(sel_diag)
    tempname S
    matrix `S' = e(se_sel_delta)
    local nm `e(snames)'
    if "`nm'" == "" local nm `e(lhs)'
    di ""
    di as txt "Selection of the buyers (Shonkwiler and Yen 1999)"
    di as txt "  E[w_i] = Phi_i f_i + delta_i phi_i, one probit per corrected good"
    di as txt "{hline 13}{c TT}{hline 64}"
    di as txt _col(14) "{c |}" _col(19) "Buyers %" _col(29) "Pseudo R2" _col(40) "Perfect" ///
        _col(51) "VIF(beta)" _col(64) "delta" _col(72) "z"
    di as txt "{hline 13}{c +}{hline 64}"
    local wvif ""
    local wper 0
    forvalues i = 1/`=colsof(`D')' {
        if missing(`D'[1, `i']) continue
        local g : word `i' of `nm'
        local z = `D'[1, `i'] / `S'[1, `i']
        di as txt %12s abbrev("`g'", 12) " {c |}" as res _col(17) %9.1f `G'[`i', 1] ///
            _col(29) %9.3f `G'[`i', 2] _col(40) %7.0f `G'[`i', 3] _col(50) %10.1f `G'[`i', 4] ///
            _col(61) %8.4f `D'[1, `i'] _col(70) %6.2f `z'
        if `G'[`i', 4] > 10 & !missing(`G'[`i', 4]) local wvif "`wvif' `g'"
        if `G'[`i', 3] > 0 & !missing(`G'[`i', 3]) local wper 1
    }
    di as txt "{hline 13}{c BT}{hline 64}"
    local zz ""
    foreach g of local nm {
        if "`e(sel_z_`g')'" != "" local zz "`zz' `g': `e(sel_z_`g')';"
    }
    if "`zz'" != "" di as txt "Probit-only variables:" as res "`zz'"
    else {
        di as txt "Probits on ln(p/x) and the demographics, with no variable of the"
        di as txt "probit only (no exclusion restriction)"
    }
    di as txt "VIF(beta): variance inflation of the expenditure coefficient due to the"
    di as txt "selection term, 1/(1-rho^2); above 10, phi is almost collinear with ln x"
    if "`wvif'" != "" {
        di as err "warning:`wvif': the correction multiplies the std. err. of the expenditure"
        di as err "  coefficient by more than 3 and rests on the curvature of the probit; leave"
        di as err "  the good uncorrected with selgoods(), or add a variable of the probit only"
        di as err "  in selvars()"
    }
    if `wper' di as err "warning: some probits predict purchase perfectly for some households"
end

* The elasticity tables through the private copy of tabstars: on the screen,
* in a file (saveres(): .tex, .csv, .md, .xlsx, .docx), or both.
program define _equaids_tables
    syntax [anything(name=fname)] , DEC(integer) DISLAS(integer) CONV(string) ///
        SHCQ(integer) [STars SCREEN(integer 1)]
    local file `fname'
    local M = e(ngoods)
    local n = cond(`dislas', `M', `M' - 1)
    local nm `e(snames)'
    if "`nm'" == "" local nm `e(lhs)'
    local nm : list retokenize nm
    local sfx = cond("`conv'" == "market", "", cond("`conv'" == "hhmean", "h", "m"))
    local hasse = e(elastse)
    if `hasse' {
        capture confirm matrix e(V_elas_x`sfx')
        local hasse = (_rc == 0)
    }
    local show = cond(`screen', "", "nodisplay")
    local st   = cond("`stars'" != "", "", "nostars")
    tempname T S
    if "`file'" != "" _equaids_tabstars export begin using "`file'", replace
    capture noisily {
        * Table 1: the budget shares at which the elasticities are evaluated
        if "`sfx'" == "" {
            matrix `T' = e(aggshare) * 100
            local tt "Aggregate budget shares (predicted, in % of total expenditure)"
        }
        else if "`sfx'" == "m" {
            matrix `T' = e(shares_m) * 100
            local tt "Predicted budget shares at the means (in %)"
        }
        else {
            matrix `T' = e(shares_h) * 100
            local tt "Mean of the predicted budget shares (in %)"
        }
        matrix `T' = `T'[1, 1..`n']
        matrix rownames `T' = "Share"
        _equaids_setnames `T', names(`nm') n(`n') cols
        _equaids_tab, num(1) title("`tt'") est(`T') dec(`dec') `show'
        * Table 2: expenditure elasticities
        matrix `T' = e(elas_x`sfx')
        matrix `T' = `T'[1, 1..`n']
        matrix rownames `T' = "Elasticity"
        _equaids_setnames `T', names(`nm') n(`n') cols
        local tt "Expenditure elasticities, d ln q / d ln x"
        if `hasse' {
            mata: _eq_semat("e(V_elas_x`sfx')", 1, "`S'")
            matrix `S' = `S'[1, 1..`n']
            _equaids_setnames `S', names(`nm') n(`n') cols
            _equaids_tab, num(2) title("`tt'") est(`T') se(`S') dec(`dec') `st' `show'
        }
        else _equaids_tab, num(2) title("`tt'") est(`T') dec(`dec') `show'
        * Table 3 (and 4): price elasticities, rows = goods, columns = prices
        local k 2
        foreach w in u c {
            if "`w'" == "c" & !`shcq' continue
            local ++k
            if "`w'" == "u" local tt "Uncompensated price elasticities, d ln q_i / d ln p_j (rows i: goods, columns j: prices)"
            else            local tt "Compensated price elasticities, at constant utility (rows i: goods, columns j: prices)"
            matrix `T' = e(elas_`w'`sfx')
            matrix `T' = `T'[1..`n', 1..`n']
            _equaids_setnames `T', names(`nm') n(`n') rows cols
            if `hasse' {
                mata: _eq_semat("e(V_elas_`w'`sfx')", `M', "`S'")
                matrix `S' = `S'[1..`n', 1..`n']
                _equaids_setnames `S', names(`nm') n(`n') rows cols
                _equaids_tab, num(`k') title("`tt'") est(`T') se(`S') dec(`dec') rowtitle("Good") `st' `show'
            }
            else _equaids_tab, num(`k') title("`tt'") est(`T') dec(`dec') rowtitle("Good") `show'
        }
    }
    if _rc {
        local rc = _rc
        capture _equaids_tabstars export clear
        if "`file'" != "" di as err "saveres(): the tables could not be written to `file'"
        exit `rc'
    }
    if "`file'" != "" _equaids_tabstars export end
end

* the first n names of the goods as row and/or column names
program define _equaids_setnames
    syntax name(name=mat) , names(string) n(integer) [rows cols]
    local l ""
    forvalues i = 1/`n' {
        local l `l' `: word `i' of `names''
    }
    if "`rows'" != "" matrix rownames `mat' = `l'
    if "`cols'" != "" matrix colnames `mat' = `l'
end

* one table (and its standard errors) through the private tabstars
program define _equaids_tab
    syntax, num(string) title(string) est(name) [se(name) dec(integer 4) noSTars rowtitle(string) NODISplay]
    if "`nodisplay'" == "" di ""
    * with a survey design the tests use t with the design degrees of freedom
    local df = cond("`e(vce)'" == "svy", e(df_r), 0)
    if "`se'" != "" _equaids_tabstars `est', se(`se') dec(`dec') title(`"`title'"') tabnumber(`num') rowtitle(`"`rowtitle'"') df(`df') `nodisplay' `stars'
    else _equaids_tabstars `est', dec(`dec') title(`"`title'"') tabnumber(`num') rowtitle(`"`rowtitle'"') `nodisplay'
end

* notes under the elasticity tables
program define _equaids_elasnotes
    syntax, conv(string) dislas(integer) shcq(integer) shck(integer)
    * aggregation identities: shown on request, and whenever one fails (it
    * would mean a wrong formula, which the user must know)
    local tol 1e-6
    * under selection, symmetry holds for the latent shares, not for the
    * expected ones: it is not an identity of the reported elasticities
    local sel = ("`e(selection)'" != "")
    local bad = (e(chk_engel) > `tol') | (e(chk_cournot) > `tol') | ///
        (e(chk_homog) > `tol') | (!`sel' & e(chk_symm) > `tol')
    if `shck' | `bad' {
        di ""
        di as txt "Aggregation identities of the aggregate elasticities (exact by construction)"
        di as txt "  Engel        |sum_i S_i E_i - 1|"            _col(52) as res %11.3e e(chk_engel)
        di as txt "  Cournot      max_j |sum_i S_i E_ij + S_j|"   _col(52) as res %11.3e e(chk_cournot)
        di as txt "  Homogeneity  max_i |sum_j E_ij + E_i|"       _col(52) as res %11.3e e(chk_homog)
        if !`sel' di as txt "  Symmetry     max |S_i E*_ij - S_j E*_ji|"    _col(52) as res %11.3e e(chk_symm)
        else      di as txt "  Symmetry     holds for the latent shares only (selection)"
        if `bad' {
            di as err "  warning: an aggregation identity fails by more than `tol';"
            di as err "  the reported elasticities should not be trusted; please report this."
        }
    }
    di ""
    if "`conv'" == "market" {
        di as txt "Market elasticities: the elasticities of total demand; each household"
        di as txt "counts in proportion to its expenditure on the good."
    }
    else if "`conv'" == "households" di as txt "Elasticities of the household, at the means of ln p, ln x and z."
    else if "`conv'" == "individuals" {
        di as txt "Elasticities of the individual: at the means of ln p, ln x and z, each household"
        di as txt "counting for its weight times its size (`e(hhsize)') in the whole estimation."
    }
    else di as txt "Mean of the household elasticities (each household counts by its weight only)."
    if e(elastse) {
        if "`e(vce)'" == "svy" local vd "by the survey design (linearized: strata, PSUs, fpc)"
        else if "`e(vce)'" == "cluster" local vd "by cluster (`e(clustvar)')"
        else local vd "robust"
        if "`e(vce)'" == "bootstrap" {
            di as txt "Standard errors: bootstrap, " e(N_reps_ok) " replications of the whole procedure"
            di as txt "(" cond("`e(pimpute)'" != "", "prices filled, ", "") ///
                cond("`e(selection)'" != "", "probits, ", "") "system estimated again on each)."
        }
        else {
            di as txt "Standard errors: analytic, `vd'; they combine the sampling of the"
            di as txt "households and the estimation of the coefficients (influence functions)."
        }
    }
    else di as txt "{bf:noelastse}: standard errors not computed."
    if "`conv'" == "hhmean" {
        * the household elasticities divide by the predicted shares
        tempname NS
        matrix `NS' = e(n_fsmall)
        local nsm 0
        forvalues i = 1/`=colsof(`NS')' {
            local nsm = max(`nsm', `NS'[1, `i'])
        }
        if `nsm' > 0 | e(n_shout) > 0 {
            di as err "warning: the household elasticities divide by the predicted shares;"
            if `nsm' > 0 di as err "  up to `nsm' households have a predicted share below 0.001 for a good,"
            if e(n_shout) > 0 di as err "  " e(n_shout) " households have predicted shares outside [0, 1];"
            di as err "  their elasticities are extreme and drive the mean.  The aggregate"
            di as err "  elasticities (the default) are not affected."
        }
    }
    if `dislas' == 0 di as txt "(last good omitted: {bf:dislas(1)} to show it)"
    if !`shcq' di as txt "{bf:compensated} adds the Hicksian table, {bf:checks} the aggregation identities."
end

* Diagnostics of the estimate, each with its reason.
program define _equaids_warnings
    if (!e(converged)) {
        di as err "convergence not achieved" _c
        if (e(stalled)) di as err ": no progress over 20 iterations (stall: neither the log likelihood nor the scaled gradient moved)"
        else            di as err ": iteration limit reached"
        if (e(n_lneg) >= e(N) / 2) di as err "  - alpha_0 = " e(anot) " is above the deflated log expenditure of " ///
            e(n_lneg) " of " e(N) " households;" _n ///
            "    choose alpha_0 below the smallest log expenditure (" %6.3f e(lnx_min) ")"
        if (e(n_boundary) > 0) di as err "  - " e(n_boundary) ///
            " steps were cut at the boundary m0(z) = 1 + rho'z > 0 (Ray scaling)"
        if (e(rcond) < 1e-10) di as err "  - the information matrix is nearly singular (scaled rcond " ///
            %8.2e e(rcond) "): some parameters are weakly identified"
        di as err "  - see the demographics diagnostics and try from() with other starting values"
    }
    else {
        if (e(rcond) < 1e-10) di as txt "note: information matrix nearly singular (scaled rcond " ///
            %8.2e e(rcond) "); some parameters are weakly identified"
        if (e(n_boundary) > 0) di as txt "note: " e(n_boundary) ///
            " steps were cut at the boundary m0(z) = 1 + rho'z > 0"
    }
    if (e(ndemos) > 0) {
        if (e(m0_min) > 10) {
            di as txt "note: m0(z) = 1 + rho'z exceeds 10 for every household (smallest " %9.3g e(m0_min) "):" _n ///
                "      the 1 of Ray's scaling is negligible, so the level of rho acts as a free" _n ///
                "      constant in ln x - ln m0 - ln a(p) and competes with alpha_0.  alpha_0 is" _n ///
                "      likely too high: choose it below the smallest log expenditure."
        }
        if (e(m0_t) < 3 & "`e(vce)'" != "bootstrap") {
            di as err "warning: the estimate lies within " %4.2f e(m0_t) " standard errors of the boundary of" _n ///
                "  Ray's scaling, m0(z) = 1 + rho'z > 0: the estimator is not regular there and" _n ///
                "  the linearized standard errors understate the uncertainty; use vce(bootstrap)," _n ///
                "  or demographics that keep m0(z) away from 0"
        }
        if (e(m0_min) < 0.05) {
            di as txt "note: the optimum lies near the boundary of the Ray scaling: m0(z) = 1 + rho'z" _n ///
                "      falls to " %8.5f e(m0_min) " (below 0.05 for " e(n_m0low) " households);" _n ///
                "      the estimates are driven by these households.  Demographics that keep" _n ///
                "      m0(z) away from 0 (nonnegative counts), or demographic translation," _n ///
                "      avoid this boundary."
        }
    }
    if (e(n_lneg) > 0) di as txt "note: ln(x/(m0 a(p))) <= 0 for " e(n_lneg) ///
        " households: their expenditure is below the estimated cost a(p) (alpha_0 = " e(anot) ")" _n ///
        "      a(p) is read as a subsistence cost (Deaton and Muellbauer 1980)"
    if (e(n_shout) > 0) di as txt "note: predicted shares outside [0, 1] for " e(n_shout) " households"
    if ("`e(vce)'" == "svy" & e(V_missing) == 1) {
        di as txt "note: missing standard errors because of " e(N_single) " stratum(s) with a single sampling" _n ///
            "      unit (svyset's singleunit(missing), Stata's default); svyset ..., singleunit(certainty)," _n ///
            "      singleunit(scaled) or singleunit(centered) chooses a treatment of these strata"
    }
    if (`"`e(data_notes)'"' != "") {
        local notes `"`e(data_notes)'"'
        while (`"`notes'"' != "") {
            gettoken n notes : notes, parse("|")
            if (`"`n'"' != "|") di as txt "note: `n'"
        }
    }
end

* Descriptive diagnostics of the demographics (weighted by om): distinct
* values, min, max, mean, sd; for a binary variable, the weighted share of
* the rarer modality and the effective size N p (1 - p), to which the
* variance of its coefficients is inversely proportional.
program define _equaids_demostats, rclass
    syntax varlist(numeric) if, om(varname) mat(name)
    marksample touse
    local K : word count `varlist'
    matrix `mat' = J(`K', 7, .)
    matrix colnames `mat' = distinct min max mean sd p_rare Neff
    matrix rownames `mat' = `varlist'
    quietly count if `touse'
    local N = r(N)
    local notes ""
    local k 0
    foreach z of local varlist {
        local ++k
        mata: st_local("nd", strofreal(rows(uniqrows(st_data(., "`z'", "`touse'")))))
        quietly summarize `z' [aw=`om'] if `touse'
        matrix `mat'[`k', 1] = `nd'
        matrix `mat'[`k', 2] = r(min)
        matrix `mat'[`k', 3] = r(max)
        matrix `mat'[`k', 4] = r(mean)
        matrix `mat'[`k', 5] = r(sd)
        if (`nd' == 2) {
            local lo = r(min)
            quietly summarize `om' if `touse' & `z' == `lo', meanonly
            local p1 = r(sum)
            quietly summarize `om' if `touse', meanonly
            local p = `p1' / r(sum)
            local prare = min(`p', 1 - `p')
            matrix `mat'[`k', 6] = `prare'
            matrix `mat'[`k', 7] = `N' * `prare' * (1 - `prare')
            local neff = `N' * `prare' * (1 - `prare')
            if (`neff' < 30) {
                local sp : di %5.3f `prare'
                local sn : di %5.1f `neff'
                local notes `"`notes'`z' is binary with a rare modality (share `sp', effective size N p(1-p) = `sn'): its coefficients are imprecise|"'
            }
        }
        else if (`nd' <= 12) {
            quietly count if `touse' & `z' != round(`z')
            if (r(N) == 0) {
                local notes `"`notes'`z' has `nd' integer values: if it is categorical, use indicator variables (Ray scaling is linear in z)|"'
            }
        }
    }
    return local notes `"`notes'"'
end

version 14.2
* ---------------------------------------------------------------------------
* estat engel (internal): the fitted shares
*   grid  on the ln x of the grid rows, prices and demographics at the vectors
*         lpm and zm; standard errors from the analytic Jacobian and
*         e(V_free); turning points of the quadratic Engel curves
*   obs   for every household, as observed
* writes <out>w1..wM (and <out>se1..seM on the grid)
program define _equaids_engel, rclass
    syntax , MODE(string) LX(varname) TOUSE(varname) OUT(string) ///
        [LP(varlist) Z(varlist) Q(varlist) LPM(name) ZM(name) QM(name) NOVAR]
    local doV = ("`novar'" == "")
    if "`mode'" == "grid" {
        mata: _eq_engel_grid("`lx'", "`touse'", "`lpm'", "`zm'", "`out'", `doV', "`qm'")
        return matrix lnx_turn = __eq_turn
        return scalar ok = `eqok'
    }
    else mata: _eq_engel_obs("`lp'", "`lx'", "`z'", "`touse'", "`out'", "`q'")
end

mata:
mata set matastrict on

// the free -> full map of _eq_fit (full(th) = c0 + th * D')
real matrix _eq_dmat(real scalar P, real scalar M, real scalar K, real scalar qd)
{
    real rowvector c0, e
    real matrix    D
    real scalar    r
    c0 = _eq_full(J(1, P, 0), M, K, qd)
    D  = J(cols(c0), P, 0)
    for (r = 1; r <= P; r++) {
        e = J(1, P, 0) ; e[r] = 1
        D[., r] = (_eq_full(e, M, K, qd) - c0)'
    }
    return(D)
}

// fitted shares (N x M) and, when wanted, their Jacobians (G, M-1 pointers)
real scalar _eq_fitted(real rowvector th, real matrix LP, real colvector LX,
                       real matrix Z, real scalar a0, real scalar qd,
                       real scalar M, real matrix D, real matrix F,
                       pointer(real matrix) rowvector G, real scalar wantG)
{
    real matrix U
    if (!_eq_model(th, J(rows(LX), M, 0), LP, LX, Z, a0, qd, D, U, G, wantG)) return(0)
    F = -U
    F = F, 1 :- rowsum(F)
    return(1)
}

// selection: the corrected goods (cix), the variables of each probit
// (zmask over e(sel_vars)) and the number of free parameters of the model
void _eq_selctx(real scalar M, real rowvector cix, real matrix zmask,
                real scalar P)
{
    string rowvector sh, sg, zv, zi
    real scalar      i, C
    sh = tokens(st_global("e(lhs)"))
    sg = tokens(st_global("e(selgoods)"))
    zv = tokens(st_global("e(sel_vars)"))
    cix = J(1, M - 1, 0) ; zmask = J(M - 1, cols(zv), 0)
    C = 0
    for (i = 1; i < M; i++) {
        if (!anyof(sg, sh[i])) continue
        cix[i] = ++C
        zi = tokens(st_global("e(sel_z_" + sh[i] + ")"))
        if (cols(zv) & cols(zi)) zmask[i, .] = colsum(J(cols(zi), 1, zv) :== J(1, cols(zv), zi'))
    }
    P = cols(st_matrix("e(b_free)")) - C
}

// expected shares E[w] (N x M) under selection, and when wanted their
// standard errors from e(V_sel_psi) (derivatives by central differences)
void _eq_selengel(real matrix LP, real colvector LX, real matrix Z,
                  real matrix Q, real scalar doV, real matrix EW,
                  real matrix SE)
{
    real rowvector psi, pp, pm, cix
    real matrix    zmask, V, DX, DP, E1, E2, GG, Gi
    real scalar    M, qd, a0, P, k, h, i, n
    M   = st_numscalar("e(ngoods)")
    qd  = (st_global("e(model)") == "QUAIDS")
    a0  = st_numscalar("e(anot)")
    psi = st_matrix("e(sel_psi)")
    _eq_selctx(M, cix, zmask, P)
    _eq_selparts(psi, LP, LX, Z, Q, a0, qd, M, P, cix, zmask, EW, DX, DP)
    SE = J(0, 0, .)
    if (!doV) return
    n  = rows(LX)
    V  = st_matrix("e(V_sel_psi)")
    GG = J(n * M, cols(psi), .)
    for (k = 1; k <= cols(psi); k++) {
        h = 1e-5 * max((1, abs(psi[k])))
        pp = psi ; pp[k] = psi[k] + h
        pm = psi ; pm[k] = psi[k] - h
        _eq_selparts(pp, LP, LX, Z, Q, a0, qd, M, P, cix, zmask, E1, DX, DP)
        _eq_selparts(pm, LP, LX, Z, Q, a0, qd, M, P, cix, zmask, E2, DX, DP)
        GG[., k] = vec((E1 - E2) :/ (2 * h))
    }
    SE = J(n, M, .)
    for (i = 1; i <= M; i++) {
        Gi = GG[|(i - 1) * n + 1, 1 \ i * n, cols(psi)|]
        SE[., i] = sqrt(rowsum((Gi * V) :* Gi))
    }
}

void _eq_engel_grid(string scalar lxv, string scalar touse, string scalar lpm,
                    string scalar zm, string scalar out, real scalar doV,
                    | string scalar qm)
{
    real matrix    EW, SE, Q
    real rowvector th, mlp, mz, tp, al, be, la, rh, f
    real matrix    V, D, F, LP, Z, g, Ga, et
    real colvector LX, se
    real scalar    M, K, qd, a0, n, i, j, ok, lna, lnm0, bc, Bi
    pointer(real matrix) rowvector G

    th  = st_matrix("e(b_free)")
    M   = st_numscalar("e(ngoods)")
    K   = st_numscalar("e(ndemos)")
    qd  = (st_global("e(model)") == "QUAIDS")
    a0  = st_numscalar("e(anot)")
    LX  = st_data(., lxv, touse)
    n   = rows(LX)
    mlp = st_matrix(lpm)
    mz  = (K > 0 ? st_matrix(zm) : J(1, 0, .))
    LP  = J(n, 1, 1) * mlp
    Z   = J(n, 1, 1) * mz
    // selection: the expected shares, at the means of the probit-only
    // variables too; no turning point (E[w] is not quadratic in ln x)
    if (st_global("e(selection)") != "") {
        if (args() > 6 & qm != "") Q = J(n, 1, 1) * st_matrix(qm)
        else                       Q = J(n, 0, .)
        _eq_selengel(LP, LX, Z, Q, doV, EW, SE)
        st_local("eqok", strofreal(!hasmissing(EW)))
        st_matrix("__eq_turn", J(1, M, .))
        for (i = 1; i <= M; i++) {
            (void) st_addvar("double", out + "w" + strofreal(i))
            st_store(., out + "w" + strofreal(i), touse, EW[., i])
            if (doV) {
                (void) st_addvar("double", out + "se" + strofreal(i))
                st_store(., out + "se" + strofreal(i), touse, SE[., i])
            }
        }
        return
    }
    D   = _eq_dmat(cols(th), M, K, qd)
    ok  = _eq_fitted(th, LP, LX, Z, a0, qd, M, D, F, G, doV)
    st_local("eqok", strofreal(ok))
    st_matrix("__eq_turn", J(1, M, .))
    if (!ok) return
    for (i = 1; i <= M; i++) {
        (void) st_addvar("double", out + "w" + strofreal(i))
        st_store(., out + "w" + strofreal(i), touse, F[., i])
    }
    if (doV) {
        V = st_matrix("e(V_free)")
        for (i = 1; i <= M; i++) {
            if (i < M) g = *G[i]
            else {
                g = J(n, cols(th), 0)
                for (j = 1; j < M; j++) g = g - *G[j]
            }
            se = sqrt(rowsum((g * V) :* g))
            (void) st_addvar("double", out + "se" + strofreal(i))
            st_store(., out + "se" + strofreal(i), touse, se)
        }
    }
    // turning points: d f_i / d l = B_i + 2 lambda_i l / (bc) = 0, so
    // l* = -B_i bc / (2 lambda_i), ln x* = l* + ln m0 + ln a (at the means)
    if (qd) {
        f = _eq_full(th, M, K, qd)
        _eq_unpack(f, M, K, qd, al, be, Ga, la, et, rh)
        lna  = a0 + mlp * al' + 0.5 * mlp * Ga * mlp'
        lnm0 = (K > 0 ? ln(1 + mz * rh') : 0)
        bc   = exp(mlp * be' + (K > 0 ? (mz * et) * mlp' : 0))
        tp   = J(1, M, .)
        for (i = 1; i <= M; i++) {
            Bi = be[i] + (K > 0 ? mz * et[., i] : 0)
            if (la[i] != 0) tp[i] = -Bi * bc / (2 * la[i]) + lnm0 + lna
        }
        st_matrix("__eq_turn", tp)
    }
}

void _eq_engel_obs(string scalar lpv, string scalar lxv, string scalar zv,
                   string scalar touse, string scalar out, | string scalar qv)
{
    real matrix    EW, SE, Q
    real rowvector th
    real matrix    LP, Z, D, F
    real colvector LX
    real scalar    M, K, qd, i
    pointer(real matrix) rowvector G

    th = st_matrix("e(b_free)")
    M  = st_numscalar("e(ngoods)")
    K  = st_numscalar("e(ndemos)")
    qd = (st_global("e(model)") == "QUAIDS")
    LP = st_data(., tokens(lpv), touse)
    LX = st_data(., lxv, touse)
    Z  = (K > 0 ? st_data(., tokens(zv), touse) : J(rows(LX), 0, .))
    if (st_global("e(selection)") != "") {
        if (args() > 5 & qv != "") Q = st_data(., tokens(qv), touse)
        else                       Q = J(rows(LX), 0, .)
        _eq_selengel(LP, LX, Z, Q, 0, EW, SE)
        for (i = 1; i <= M; i++) {
            (void) st_addvar("double", out + "w" + strofreal(i))
            st_store(., out + "w" + strofreal(i), touse, EW[., i])
        }
        return
    }
    D  = _eq_dmat(cols(th), M, K, qd)
    if (!_eq_fitted(th, LP, LX, Z, st_numscalar("e(anot)"), qd, M, D, F, G, 0)) {
        errprintf("equaids: m0(z) <= 0 for some household at the estimate\n")
        exit(459)
    }
    for (i = 1; i <= M; i++) {
        (void) st_addvar("double", out + "w" + strofreal(i))
        st_store(., out + "w" + strofreal(i), touse, F[., i])
    }
}

// ---------------------------------------------------------------------------
// free parameters -> full vector (Poi's parameterization)
//   free: alpha(M-1) beta(M-1) vech(Gamma[1..M-1,1..M-1]) lambda(M-1)
//         eta (K x (M-1), by demographic) rho(K)
//   full: alpha(M) beta(M) vech(Gamma) lambda(M) vec(eta')' rho(K)
// ---------------------------------------------------------------------------
real rowvector _eq_full(real rowvector th, real scalar M, real scalar K,
                        real scalar qd)
{
    real rowvector al, be, la, rh, out
    real matrix    Ga, et
    real scalar    c, i, j, d

    c  = 1
    al = J(1, M, 0) ; al[M] = 1
    for (i = 1; i < M; i++) {
        al[i] = th[c++]
        al[M] = al[M] - al[i]
    }
    be = J(1, M, 0)
    for (i = 1; i < M; i++) {
        be[i] = th[c++]
        be[M] = be[M] - be[i]
    }
    Ga = J(M, M, 0)
    for (j = 1; j < M; j++) {
        for (i = j; i < M; i++) {
            Ga[i, j] = th[c]
            Ga[j, i] = th[c]
            c++
        }
    }
    for (i = 1; i < M; i++) {
        for (j = 1; j < M; j++) Ga[i, M] = Ga[i, M] - Ga[i, j]
        Ga[M, i] = Ga[i, M]
    }
    for (i = 1; i < M; i++) Ga[M, M] = Ga[M, M] - Ga[i, M]
    out = al, be, vech(Ga)'
    if (qd) {
        la = J(1, M, 0)
        for (i = 1; i < M; i++) {
            la[i] = th[c++]
            la[M] = la[M] - la[i]
        }
        out = out, la
    }
    if (K > 0) {
        et = J(K, M, 0)
        for (d = 1; d <= K; d++) {
            for (j = 1; j < M; j++) {
                et[d, j] = th[c++]
                et[d, M] = et[d, M] - et[d, j]
            }
        }
        rh = th[|c \ c + K - 1|]
        out = out, vec(et')', rh
    }
    return(out)
}

void _eq_unpack(real rowvector f, real scalar M, real scalar K,
                real scalar qd, real rowvector al, real rowvector be,
                real matrix Ga, real rowvector la, real matrix et,
                real rowvector rh)
{
    real scalar c, nv
    nv = M * (M + 1) / 2
    al = f[|1 \ M|]
    be = f[|M + 1 \ 2 * M|]
    Ga = invvech(f[|2 * M + 1 \ 2 * M + nv|]')
    c  = 2 * M + nv + 1
    if (qd) {
        la = f[|c \ c + M - 1|]
        c = c + M
    }
    else la = J(1, M, 0)
    if (K > 0) {
        et = colshape(f[|c \ c + K * M - 1|], M)
        c = c + K * M
        rh = f[|c \ c + K - 1|]
    }
    else {
        et = J(0, M, 0)
        rh = J(1, 0, 0)
    }
}

// ---------------------------------------------------------------------------
// model: residuals of the first M-1 equations and, when wanted, their
// analytic Jacobians with respect to the free parameters (G[i] = N x P)
// ---------------------------------------------------------------------------
real scalar _eq_model(real rowvector th, real matrix W, real matrix LP,
                      real colvector LX, real matrix Z, real scalar a0,
                      real scalar qd, real matrix D, real matrix U,
                      pointer(real matrix) rowvector G, real scalar wantG,
                      | real matrix PH, real matrix PHI, real rowvector cix)
{
    real scalar    M, K, N, i, r, c, d, k, col, Pf, sel, C, dl, P0
    real rowvector al, be, la, rh, f
    real matrix    Ga, et, ZE, Gf, Gs
    real colvector fi
    real colvector lna, lnb, lnc, m0, bc, ell, q, B, mu

    N = rows(W) ; M = cols(W) ; K = cols(Z)
    // selection (Shonkwiler and Yen): E[w_i] = Phi_i f_i + delta_i phi_i
    sel = (args() > 11)
    P0 = cols(D)
    C  = cols(th) - P0
    f = _eq_full(th, M, K, qd)
    _eq_unpack(f, M, K, qd, al, be, Ga, la, et, rh)
    lna = a0 :+ LP * al' :+ 0.5 :* rowsum((LP * Ga) :* LP)
    lnb = LP * be'
    if (K > 0) {
        ZE  = Z * et
        lnc = rowsum(ZE :* LP)
        m0  = 1 :+ Z * rh'
        if (min(m0) <= 0) return(0)
    }
    else {
        ZE  = J(N, M, 0)
        lnc = J(N, 1, 0)
        m0  = J(N, 1, 1)
    }
    bc  = exp(lnb + lnc)
    ell = LX - ln(m0) - lna
    q   = (ell :^ 2) :/ bc
    U   = J(N, M - 1, .)
    if (wantG) {
        Pf = cols(f)
        G  = J(1, M - 1, NULL)
    }
    for (i = 1; i < M; i++) {
        B  = be[i] :+ ZE[., i]
        mu = B
        U[., i] = W[., i] - (al[i] :+ LP * Ga[., i] + B :* ell)
        if (qd) {
            U[., i] = U[., i] - la[i] :* q
            mu = mu + 2 * la[i] :* ell :/ bc
        }
        if (sel) {
            fi = W[., i] - U[., i]
            dl = (cix[i] ? th[P0 + cix[i]] : 0)
            U[., i] = W[., i] - (PH[., i] :* fi + dl :* PHI[., i])
        }
        if (!wantG) continue
        Gf  = J(N, Pf, 0)
        // alpha
        Gf[|1, 1 \ N, M|] = -mu :* LP
        Gf[., i] = Gf[., i] :+ 1
        // beta
        if (qd) Gf[|1, M + 1 \ N, 2 * M|] = -la[i] :* q :* LP
        Gf[., M + i] = Gf[., M + i] + ell
        // gamma (vech order: column c, rows r >= c)
        col = 2 * M
        for (c = 1; c <= M; c++) {
            for (r = c; r <= M; r++) {
                col++
                if (r == c) Gf[., col] = -mu :* (0.5 :* LP[., r] :^ 2)
                else        Gf[., col] = -mu :* (LP[., r] :* LP[., c])
                if (r == c) {
                    if (r == i) Gf[., col] = Gf[., col] + LP[., r]
                }
                else {
                    if (r == i) Gf[., col] = Gf[., col] + LP[., c]
                    if (c == i) Gf[., col] = Gf[., col] + LP[., r]
                }
            }
        }
        // lambda
        if (qd) {
            Gf[., col + i] = q
            col = col + M
        }
        // eta, then rho
        if (K > 0) {
            for (d = 1; d <= K; d++) {
                for (k = 1; k <= M; k++) {
                    col++
                    if (qd) Gf[., col] = -la[i] :* q :* Z[., d] :* LP[., k]
                    if (k == i) Gf[., col] = Gf[., col] + Z[., d] :* ell
                }
            }
            for (d = 1; d <= K; d++) {
                col++
                Gf[., col] = -mu :* Z[., d] :/ m0
            }
        }
        if (sel) {
            Gs = PH[., i] :* (Gf * D)
            if (C) {
                Gs = Gs, J(N, C, 0)
                if (cix[i]) Gs[., P0 + cix[i]] = PHI[., i]
            }
            G[i] = &(Gs :* 1)
        }
        else G[i] = &(Gf * D)
    }
    return(1)
}

// ---------------------------------------------------------------------------
// elasticities: the pieces, for all M goods, at given data rows
//   F   N x M      predicted shares f_i
//   MU  N x M      mu_i  = d f_i / d ln x = B_i + 2 lambda_i l/(bc)
//   MUJ N x M^2    mu_ij = d f_i / d ln p_j (column (i-1)M + j)
//                        = gamma_ij - mu_i pi_j - lambda_i B_j l^2/(bc)
//   with B_i = beta_i + eta_i'z and pi_j = alpha_j + sum_k gamma_jk ln p_k
// household elasticities, predicted shares:
//   e_i = 1 + mu_i/f_i ; e_ij = -delta_ij + mu_ij/f_i ; e*_ij = e_ij + e_i f_j
// ---------------------------------------------------------------------------
void _eq_parts(real rowvector th, real matrix LP, real colvector LX,
               real matrix Z, real scalar a0, real scalar qd, real scalar M,
               real matrix F, real matrix MU, real matrix MUJ)
{
    real scalar    N, K, i, j
    real rowvector al, be, la, rh, f
    real matrix    Ga, et, ZE, B, PI
    real colvector lna, lnc, m0, bc, ell, q

    N = rows(LP) ; K = cols(Z)
    f = _eq_full(th, M, K, qd)
    _eq_unpack(f, M, K, qd, al, be, Ga, la, et, rh)
    lna = a0 :+ LP * al' :+ 0.5 :* rowsum((LP * Ga) :* LP)
    if (K > 0) {
        ZE  = Z * et
        lnc = rowsum(ZE :* LP)
        m0  = 1 :+ Z * rh'
    }
    else {
        ZE  = J(N, M, 0)
        lnc = J(N, 1, 0)
        m0  = J(N, 1, 1)
    }
    bc  = exp(LP * be' + lnc)
    ell = LX - ln(m0) - lna
    q   = (ell :^ 2) :/ bc
    B   = be :+ ZE                         // N x M
    PI  = al :+ LP * Ga                     // N x M (Ga symmetric)
    F   = PI + B :* ell + q * la
    MU  = B + 2 :* ((ell :/ bc) * la)
    MUJ = J(N, M * M, .)
    for (i = 1; i <= M; i++) {
        for (j = 1; j <= M; j++) {
            // explicit parentheses: in Mata, :- binds less tightly than -
            MUJ[., (i - 1) * M + j] = (Ga[i, j] :- (MU[., i] :* PI[., j])) - (la[i] :* (B[., j] :* q))
        }
    }
}

// the three summaries of the elasticities (point values)
//   aggregate (default): weights omega x f_i, the expenditure on the good:
//     E_i = 1 + sum(w x mu_i)/S_i, E_ij = -delta_ij + sum(w x mu_ij)/S_i,
//     E*_ij = E_ij + sum(w x (f_i + mu_i) f_j)/S_i,  S_i = sum(w x f_i)
//   at the means: the household at the weighted means of ln p, ln x, z
//   household mean: sum(w e_h)/sum(w)
void _eq_elas(real rowvector th, real matrix LP, real colvector LX,
              real matrix Z, real colvector om, real scalar a0,
              real scalar qd, real scalar M, string scalar base)
{
    real matrix    F, MU, MUJ, Fm, MUm, MUJm, Eu, Ec, Eum, Ecm, Euh, Ech, D0
    real colvector wx
    real rowvector S, Ex, Exm, Exh, mlp, mz, nsm
    real scalar    N, i, j, sw, mlx

    N  = rows(LP)
    sw = sum(om)
    _eq_parts(th, LP, LX, Z, a0, qd, M, F, MU, MUJ)
    D0 = I(M)
    // aggregate
    wx = om :* exp(LX)
    S  = colsum(wx :* F)
    Ex = 1 :+ colsum(wx :* MU) :/ S
    Eu = J(M, M, .) ; Ec = J(M, M, .)
    for (i = 1; i <= M; i++) {
        for (j = 1; j <= M; j++) {
            Eu[i, j] = -D0[i, j] + sum(wx :* MUJ[., (i - 1) * M + j]) / S[i]
            Ec[i, j] = Eu[i, j] + sum(wx :* (F[., i] + MU[., i]) :* F[., j]) / S[i]
        }
    }
    // at the means
    mlp = colsum(om :* LP) :/ sw
    mlx = sum(om :* LX) / sw
    if (cols(Z) > 0) mz = colsum(om :* Z) :/ sw
    else             mz = J(1, 0, .)
    _eq_parts(th, mlp, mlx, mz, a0, qd, M, Fm, MUm, MUJm)
    Exm = 1 :+ MUm :/ Fm
    Eum = J(M, M, .) ; Ecm = J(M, M, .)
    for (i = 1; i <= M; i++) {
        for (j = 1; j <= M; j++) {
            Eum[i, j] = -D0[i, j] + MUJm[1, (i - 1) * M + j] / Fm[i]
            Ecm[i, j] = Eum[i, j] + Exm[i] * Fm[j]
        }
    }
    // household mean
    Exh = colsum(om :* (1 :+ MU :/ F)) :/ sw
    Euh = J(M, M, .) ; Ech = J(M, M, .)
    for (i = 1; i <= M; i++) {
        for (j = 1; j <= M; j++) {
            Euh[i, j] = sum(om :* (-D0[i, j] :+ MUJ[., (i - 1) * M + j] :/ F[., i])) / sw
            Ech[i, j] = sum(om :* (-D0[i, j] :+ MUJ[., (i - 1) * M + j] :/ F[., i]
                         + (1 :+ MU[., i] :/ F[., i]) :* F[., j])) / sw
        }
    }
    nsm = colsum(F :< 0.001)
    // the shares at which the elasticities are evaluated: at the means, and
    // the mean of the predicted shares
    st_matrix(base + "Sm", Fm)
    st_matrix(base + "Sh", colsum(om :* F) :/ sw)
    st_matrix(base + "x",  Ex)
    st_matrix(base + "u",  Eu)
    st_matrix(base + "c",  Ec)
    st_matrix(base + "S",  S :/ sum(S))
    st_matrix(base + "xm", Exm)
    st_matrix(base + "um", Eum)
    st_matrix(base + "cm", Ecm)
    st_matrix(base + "xh", Exh)
    st_matrix(base + "uh", Euh)
    st_matrix(base + "ch", Ech)
    st_matrix(base + "ns", nsm)
}

// ---------------------------------------------------------------------------
// standard errors of the aggregate elasticities, analytic
//
// Each aggregate elasticity is a ratio of weighted totals, E = T/S, of
// functions of theta.  Its influence function has two terms:
//   sample term     omega_h x_h (t_h - E s_h) / S      (the households drawn)
//   parameter term  J IF_theta,h,  J = (dT/dtheta - E dS/dtheta) / S
// with IF_theta,h = A^-1 omega_h G_h' Sigma^-1 u_h (FGNLS score).  All
// derivatives are analytic: they reduce to weighted column sums of
//   DL = d l / d theta_full  and  DB = d ln(b c) / d theta_full  (N x Pf):
//   d mu_i   = e_beta_i + z e_eta_.i + (2l/bc) e_lambda_i
//              + (2 lambda_i/bc) DL - (2 lambda_i l/bc) DB
//   d pi_j   = e_alpha_j + sum_k ln p_k e_gamma_jk
//   d B_j    = e_beta_j + z e_eta_.j
//   d q      = (2l/bc) DL - q DB
//   d mu_ij  = e_gamma_ij - pi_j d mu_i - mu_i d pi_j - B_j q e_lambda_i
//              - lambda_i q d B_j - lambda_i B_j d q
// Variances: sum of squared influence functions (robust); the same influence
// functions are aggregated by cluster and by design in step 3.
// ---------------------------------------------------------------------------
real scalar _eq_vidx(real scalar r, real scalar c, real scalar M)
{
    real scalar a, b, k, idx
    a = max((r, c)) ; b = min((r, c))
    idx = 0
    for (k = 1; k < b; k++) idx = idx + (M - k + 1)
    return(idx + a - b + 1)
}

void _eq_elasvar(real rowvector th, real matrix LP, real colvector LX,
                 real matrix Z, real colvector om, real scalar a0,
                 real scalar qd, real scalar M, real matrix D, real matrix U,
                 pointer(real matrix) rowvector G, real matrix S,
                 real matrix A, string scalar base, real scalar doe,
                 real matrix DS, real scalar vmode, string scalar single)
{
    real scalar    N, K, P, Pf, nv, ob, oG, oL, oE, oR, i, j, k, d, c
    real rowvector al, be, la, rh, f, Sg, Ex, sT, sS, sC, sE
    real matrix    Ga, et, ZE, B, PI, F, MU, MUJ, DL, DB, IFt, GG, Ai,
                   IFx, IFu, IFc, Jx, Ju, Jc, US
    real colvector lna, lnc, m0, bc, ell, q, wx, v, Ex_i
    pointer(real matrix) rowvector GF

    N = rows(LP) ; K = cols(Z) ; P = cols(th)
    f = _eq_full(th, M, K, qd)
    Pf = cols(f)
    _eq_unpack(f, M, K, qd, al, be, Ga, la, et, rh)
    nv = M * (M + 1) / 2
    ob = M ; oG = 2 * M ; oL = 2 * M + nv
    oE = oL + qd * M ; oR = oE + K * M
    // pieces
    lna = a0 :+ LP * al' :+ 0.5 :* rowsum((LP * Ga) :* LP)
    if (K > 0) {
        ZE = Z * et ; lnc = rowsum(ZE :* LP) ; m0 = 1 :+ Z * rh'
    }
    else {
        ZE = J(N, M, 0) ; lnc = J(N, 1, 0) ; m0 = J(N, 1, 1)
    }
    bc  = exp(LP * be' + lnc)
    ell = LX - ln(m0) - lna
    q   = (ell :^ 2) :/ bc
    B   = ZE :+ be
    PI  = LP * Ga :+ al
    F   = PI + B :* ell + q * la
    MU  = B + 2 :* ((ell :/ bc) * la)
    // DL and DB
    DL = J(N, Pf, 0) ; DB = J(N, Pf, 0)
    for (k = 1; k <= M; k++) {
        DL[., k] = -LP[., k]
        DB[., ob + k] = LP[., k]
    }
    for (c = 1; c <= M; c++) {
        for (k = c; k <= M; k++) {
            if (k == c) DL[., oG + _eq_vidx(k, c, M)] = -0.5 :* LP[., k] :^ 2
            else        DL[., oG + _eq_vidx(k, c, M)] = -LP[., k] :* LP[., c]
        }
    }
    for (d = 1; d <= K; d++) {
        for (k = 1; k <= M; k++) DB[., oE + (d - 1) * M + k] = Z[., d] :* LP[., k]
        DL[., oR + d] = -Z[., d] :/ m0
    }
    // derivatives of the fitted shares in the free space, all M goods
    GF = J(1, M, NULL)
    for (i = 1; i < M; i++) GF[i] = G[i]
    GG = J(N, P, 0)
    for (i = 1; i < M; i++) GG = GG - *G[i]
    GF[M] = &GG
    // influence functions of theta (N x P)
    Ai = invsym(A)
    US = U * S
    IFt = J(N, P, 0)
    for (i = 1; i < M; i++) IFt = IFt + (*G[i]) :* (om :* US[., i])
    IFt = IFt * Ai
    st_matrix(base + "Vt", _eq_vagg(IFt, DS, vmode, single))
    // saveif(stub), for the tests only: the influence functions as variables
    _eq_saveif(IFt, "t")
    // noelastse: the robust variance of the coefficients only
    if (!doe) return
    // weights and totals
    wx = om :* exp(LX)
    Sg = colsum(wx :* F)
    Ex = 1 :+ colsum(wx :* MU) :/ Sg
    IFx = J(N, M, .) ; Jx = J(M, P, .)
    IFu = J(N, M * M, .) ; Ju = J(M * M, P, .)
    IFc = J(N, M * M, .) ; Jc = J(M * M, P, .)
    for (i = 1; i < M + 1; i++) {
        sS = wx' * (*GF[i])                                   // dS_i (free)
        // d T_i, full space: weighted sums of d mu_i
        sT = _eq_dmu(wx, i, la, ell, bc, Z, DL, DB, M, ob, oL, oE, qd) * D
        Jx[i, .] = (sT - (Ex[i] - 1) :* sS) :/ Sg[i]
        IFx[., i] = wx :* (MU[., i] - (Ex[i] - 1) :* F[., i]) :/ Sg[i] + IFt * Jx[i, .]'
        for (j = 1; j <= M; j++) {
            c = (i - 1) * M + j
            // T_ij = sum wx mu_ij, full-space derivative
            sT = J(1, Pf, 0)
            sT[oG + _eq_vidx(i, j, M)] = sum(wx)
            sT = sT - _eq_dmu(wx :* PI[., j], i, la, ell, bc, Z, DL, DB, M, ob, oL, oE, qd)
            v = wx :* MU[., i]                                // - mu_i d pi_j
            sT[j] = sT[j] - sum(v)
            for (k = 1; k <= M; k++) sT[oG + _eq_vidx(j, k, M)] = sT[oG + _eq_vidx(j, k, M)] - sum(v :* LP[., k])
            if (qd) {
                sT[oL + i] = sT[oL + i] - sum(wx :* B[., j] :* q)
                v = la[i] :* wx :* q                          // - lambda_i q d B_j
                sT[ob + j] = sT[ob + j] - sum(v)
                for (d = 1; d <= K; d++) sT[oE + (d - 1) * M + j] = sT[oE + (d - 1) * M + j] - sum(v :* Z[., d])
                v = la[i] :* wx :* B[., j]                    // - lambda_i B_j d q
                sT = sT - ((v :* (2 :* ell :/ bc))' * DL - (v :* q)' * DB)
            }
            sT = sT * D
            MUJ = (Ga[i, j] :- (MU[., i] :* PI[., j])) - (la[i] :* (B[., j] :* q))
            Ju[c, .] = (sT - ((sum(wx :* MUJ) / Sg[i]) :* sS)) :/ Sg[i]
            IFu[., c] = wx :* (MUJ - (sum(wx :* MUJ) / Sg[i]) :* F[., i]) :/ Sg[i] + IFt * Ju[c, .]'
            // compensated: E*_ij = E_ij + C_ij/S_i, C_ij = sum wx (f_i + mu_i) f_j
            sC = (wx :* F[., j])' * (*GF[i]) +
                 _eq_dmu(wx :* F[., j], i, la, ell, bc, Z, DL, DB, M, ob, oL, oE, qd) * D +
                 (wx :* (F[., i] + MU[., i]))' * (*GF[j])
            v = wx :* (F[., i] + MU[., i]) :* F[., j]
            Jc[c, .] = Ju[c, .] + (sC - ((sum(v) / Sg[i]) :* sS)) :/ Sg[i]
            IFc[., c] = IFu[., c] + (v - (sum(v) / Sg[i]) :* wx :* F[., i]) :/ Sg[i] +
                        IFt * ((sC - ((sum(v) / Sg[i]) :* sS)) :/ Sg[i])'
        }
    }
    _eq_saveif(IFx, "x")
    _eq_saveif(IFu, "u")
    st_matrix(base + "Vx", _eq_vagg(IFx, DS, vmode, single))
    st_matrix(base + "Vu", _eq_vagg(IFu, DS, vmode, single))
    st_matrix(base + "Vc", _eq_vagg(IFc, DS, vmode, single))
    st_matrix(base + "Jx", Jx)
    st_matrix(base + "Ju", Ju)
    st_matrix(base + "Jc", Jc)
    // at the means (households, individuals) and the household mean
    _eq_elasvar_mh(th, LP, LX, Z, om, a0, qd, M, IFt, DS, vmode, single, base)
}

// ---------------------------------------------------------------------------
// the elasticities at the means (households, individuals) and the mean of the
// household elasticities (hhmean), and their variances
//
// at the means: E(theta, m), m = the weighted means of ln p, ln x and z;
//   IF_h = J_theta IF_theta,h + J_m omega_h (m_h - m) / sum omega
// household mean: E = sum omega e_h(theta) / sum omega;
//   IF_h = omega_h (e_h - E) / sum omega + J_theta IF_theta,h
// J_theta and J_m by central differences (step 1e-5 max(1, |v|): error of the
// order of 1e-10, well below what a standard error needs)
// the elasticities are stacked as (E_x (M), E_u (M^2), E_c (M^2)), row-major
// ---------------------------------------------------------------------------
real rowvector _eq_elas_at(real rowvector th, real rowvector lp, real scalar lx,
                           real rowvector z, real scalar a0, real scalar qd,
                           real scalar M)
{
    real matrix    F, MU, MUJ
    real rowvector ex, eu, ec
    real scalar    i, j, c
    _eq_parts(th, lp, lx, z, a0, qd, M, F, MU, MUJ)
    ex = 1 :+ MU :/ F
    eu = J(1, M * M, .) ; ec = J(1, M * M, .)
    for (i = 1; i <= M; i++) {
        for (j = 1; j <= M; j++) {
            c = (i - 1) * M + j
            eu[c] = -(i == j) + MUJ[1, c] / F[i]
            ec[c] = eu[c] + ex[i] * F[j]
        }
    }
    return((ex, eu, ec))
}

// the elasticities of every household, N x (M + 2 M^2)
real matrix _eq_elas_hh(real rowvector th, real matrix LP, real colvector LX,
                        real matrix Z, real scalar a0, real scalar qd,
                        real scalar M)
{
    real matrix    F, MU, MUJ, H
    real scalar    i, j, c
    _eq_parts(th, LP, LX, Z, a0, qd, M, F, MU, MUJ)
    H = J(rows(LP), M + 2 * M * M, .)
    H[., 1..M] = 1 :+ MU :/ F
    for (i = 1; i <= M; i++) {
        for (j = 1; j <= M; j++) {
            c = (i - 1) * M + j
            H[., M + c] = -(i == j) :+ MUJ[., c] :/ F[., i]
            H[., M + M * M + c] = H[., M + c] + H[., i] :* F[., j]
        }
    }
    return(H)
}

void _eq_elasvar_mh(real rowvector th, real matrix LP, real colvector LX,
                    real matrix Z, real colvector om, real scalar a0,
                    real scalar qd, real scalar M, real matrix IFt,
                    real matrix DS, real scalar vmode, string scalar single,
                    string scalar base)
{
    real scalar    N, K, P, q, k, h, sw, nm
    real rowvector m, mp, mm, tp, tm, E0
    real matrix    Jt, Jm, MV, IF, H, V
    N = rows(LP) ; K = cols(Z) ; P = cols(th) ; q = M + 2 * M * M
    sw = sum(om)
    // at the means
    MV = (LP, LX, Z)
    nm = cols(MV)
    m  = colsum(om :* MV) :/ sw
    Jt = J(q, P, .)
    for (k = 1; k <= P; k++) {
        h = 1e-5 * max((1, abs(th[k])))
        tp = th ; tp[k] = th[k] + h
        tm = th ; tm[k] = th[k] - h
        Jt[., k] = ((_eq_elas_at(tp, m[1..M], m[M+1], (K ? m[M+2..nm] : J(1, 0, .)), a0, qd, M) -
                     _eq_elas_at(tm, m[1..M], m[M+1], (K ? m[M+2..nm] : J(1, 0, .)), a0, qd, M)) :/ (2 * h))'
    }
    Jm = J(q, nm, .)
    for (k = 1; k <= nm; k++) {
        h = 1e-5 * max((1, abs(m[k])))
        mp = m ; mp[k] = m[k] + h
        mm = m ; mm[k] = m[k] - h
        Jm[., k] = ((_eq_elas_at(th, mp[1..M], mp[M+1], (K ? mp[M+2..nm] : J(1, 0, .)), a0, qd, M) -
                     _eq_elas_at(th, mm[1..M], mm[M+1], (K ? mm[M+2..nm] : J(1, 0, .)), a0, qd, M)) :/ (2 * h))'
    }
    IF = IFt * Jt' + ((om :* (MV :- m)) :/ sw) * Jm'
    V  = _eq_vagg(IF, DS, vmode, single)
    st_matrix(base + "Vxm", V[1..M, 1..M])
    st_matrix(base + "Vum", V[M+1..M+M*M, M+1..M+M*M])
    st_matrix(base + "Vcm", V[M+M*M+1..q, M+M*M+1..q])
    // household mean
    H  = _eq_elas_hh(th, LP, LX, Z, a0, qd, M)
    E0 = colsum(om :* H) :/ sw
    Jt = J(q, P, .)
    for (k = 1; k <= P; k++) {
        h = 1e-5 * max((1, abs(th[k])))
        tp = th ; tp[k] = th[k] + h
        tm = th ; tm[k] = th[k] - h
        Jt[., k] = ((colsum(om :* _eq_elas_hh(tp, LP, LX, Z, a0, qd, M)) -
                     colsum(om :* _eq_elas_hh(tm, LP, LX, Z, a0, qd, M))) :/ (2 * h * sw))'
    }
    IF = (om :* (H :- E0)) :/ sw + IFt * Jt'
    V  = _eq_vagg(IF, DS, vmode, single)
    st_matrix(base + "Vxh", V[1..M, 1..M])
    st_matrix(base + "Vuh", V[M+1..M+M*M, M+1..M+M*M])
    st_matrix(base + "Vch", V[M+M*M+1..q, M+M*M+1..q])
}

// ---------------------------------------------------------------------------
// variance from influence functions IF (N x q), one row per household
//   vmode 0  robust    N/(N-1) sum IF_h IF_h'
//   vmode 1  cluster   G/(G-1) sum z_g z_g',  z_g = sum of IF over cluster g
//   vmode 2  svy       sum_h (1-f_h) n_h/(n_h-1) sum_i (z_hi - zbar_h)(...)'
//                      z_hi = sum of IF over PSU i of stratum h (first stage,
//                      linearization as Stata's svy); f_h = fpc if <= 1 (a
//                      sampling rate), n_h/fpc if > 1 (a population size)
// The IF already carry the weights (their scale cancels), so the PSU totals
// are the linearized totals of the design.  DS = (stratum, PSU, fpc).
// Strata with one PSU (svyset's singleunit()): certainty = no contribution;
// scaled = the variance of the other strata scaled by H/(H - H1);
// centered = deviation from the grand mean of the PSU totals.
// ---------------------------------------------------------------------------
real matrix _eq_vagg(real matrix IF, real matrix DS, real scalar vmode,
                     string scalar single)
{
    real scalar    N, G, q, h, nh, f, H, H1
    real colvector p, chg, psuid, zst, zfp
    real matrix    D2, info, Z, info2, Zh, V, zall
    N = rows(IF)
    q = cols(IF)
    if (vmode == 0) return(N / (N - 1) :* cross(IF, IF))
    p  = order(DS[., (1, 2)], (1, 2))
    D2 = DS[p, .]
    if (N > 1) chg = 1 \ ((D2[|2, 1 \ N, 1|] :!= D2[|1, 1 \ N - 1, 1|]) :|
                          (D2[|2, 2 \ N, 2|] :!= D2[|1, 2 \ N - 1, 2|]))
    else       chg = 1
    psuid = runningsum(chg)
    info  = panelsetup(psuid, 1)
    Z     = panelsum(IF[p, .], info)
    G     = rows(Z)
    if (vmode == 1) return(G / (G - 1) :* cross(Z, Z))
    zst   = D2[info[., 1], 1]
    zfp   = D2[info[., 1], 3]
    info2 = panelsetup(zst, 1)
    H     = rows(info2)
    H1    = 0
    zall  = mean(Z)
    V     = J(q, q, 0)
    for (h = 1; h <= H; h++) {
        Zh = panelsubmatrix(Z, h, info2)
        nh = rows(Zh)
        f  = zfp[info2[h, 1]]
        if (f >= .)     f = 0
        else if (f > 1) f = nh / f
        if (nh == 1) {
            H1++
            if (single == "centered") V = V + (1 - f) :* cross(Zh :- zall, Zh :- zall)
            else if (single == "missing") return(J(q, q, .))
            continue
        }
        V = V + (1 - f) * nh / (nh - 1) :* cross(Zh :- mean(Zh), Zh :- mean(Zh))
    }
    if (single == "scaled" & H1 > 0 & H > H1) V = V :* (H / (H - H1))
    return(V)
}

// saveif(stub): store the influence functions in stub<k><j> (tests only);
// reads the locals saveif and touse of the calling ado
void _eq_saveif(real matrix IF, string scalar k)
{
    string scalar  stub, tu
    string rowvector nm
    real scalar    j
    stub = st_local("saveif")
    if (stub == "") return
    tu = st_local("touse")
    nm = J(1, cols(IF), "")
    for (j = 1; j <= cols(IF); j++) {
        nm[j] = stub + k + strofreal(j)
        (void) st_addvar("double", nm[j])
    }
    st_store(., nm, tu, IF)
}

// numbers of strata, PSUs and single-PSU strata of the design (locals)
void _eq_desinfo(string scalar desv, string scalar touse)
{
    real matrix    DS, u, info
    real colvector nps
    DS = st_data(., tokens(desv), touse)
    u  = uniqrows(DS[., (1, 2)])
    info = panelsetup(u[., 1], 1)
    nps  = info[., 2] - info[., 1] :+ 1
    st_local("N_strata", strofreal(rows(info)))
    st_local("N_psu",    strofreal(rows(u)))
    st_local("N_single", strofreal(sum(nps :== 1)))
}

// weighted sums of d mu_i / d theta_full: v' d mu_i (1 x Pf)
real rowvector _eq_dmu(real colvector v, real scalar i, real rowvector la,
                       real colvector ell, real colvector bc, real matrix Z,
                       real matrix DL, real matrix DB, real scalar M,
                       real scalar ob, real scalar oL, real scalar oE,
                       real scalar qd)
{
    real rowvector s
    real scalar    d
    s = J(1, cols(DL), 0)
    if (qd) s = (v :* (2 * la[i] :/ bc))' * DL - (v :* (2 * la[i] :* ell :/ bc))' * DB
    s[ob + i] = s[ob + i] + sum(v)
    for (d = 1; d <= cols(Z); d++) s[oE + (d - 1) * M + i] = s[oE + (d - 1) * M + i] + sum(v :* Z[., d])
    if (qd) s[oL + i] = s[oL + i] + sum(v :* (2 :* ell :/ bc))
    return(s)
}

// derivative of residuals is minus the derivative of the fitted shares:
// the G returned above are derivatives of the fitted shares f_i.

real scalar _eq_obj(real matrix U, real colvector om, real matrix S)
{
    return(sum(om :* rowsum((U * S) :* U)))
}

// one Gauss-Newton step with step-halving at a given S = Sigma^-1;
// updates th, U and G, returns the relative change of th, and in nd the
// scaled gradient g'A^-1 g before the step (the Newton decrement: the
// second-order prediction of the decrease of the objective, which is -2
// times the log likelihood up to a constant; it does not depend on units)
real scalar _eq_step(real rowvector th, real matrix S, real matrix W,
                     real matrix LP, real colvector LX, real matrix Z,
                     real scalar a0, real scalar qd, real matrix D,
                     real colvector om, real matrix U,
                     pointer(real matrix) rowvector G, real scalar nd,
                     real scalar tt, real scalar hb,
                     | real matrix PH, real matrix PHI, real rowvector cix)
{
    real matrix    A, Un, H
    real colvector g, step
    real rowvector thn
    real scalar    P, M, i, j, t, Q, Qn, ok, r
    pointer(real matrix) rowvector Gn

    hb = 0          // 1 if a trial step crossed the boundary m0(z) > 0

    P = cols(th) ; M = cols(W)
    // A = sum_ij S_ij G_i'WG_j, g = sum_ij S_ij G_i'W u_j, grouped (S is
    // symmetric): H_i = sum_j S_ij G_j, A = sum_i G_i'W H_i, g = sum_i H_i'W u_i
    // -- M-1 cross products of N x P instead of (M-1)^2
    A = J(P, P, 0) ; g = J(P, 1, 0)
    for (i = 1; i < M; i++) {
        H = J(rows(W), P, 0)
        for (j = 1; j < M; j++) H = H + S[i, j] :* *G[j]
        A = A + cross(*G[i], om, H)
        g = g + cross(H, om, U[., i])
    }
    step = cholsolve(A, g)
    if (hasmissing(step)) step = invsym(A) * g
    nd = g' * step
    Q = _eq_obj(U, om, S)
    t = 1
    do {
        thn = th + t :* step'
        if (args() > 15) ok = _eq_model(thn, W, LP, LX, Z, a0, qd, D, Un, Gn, 0, PH, PHI, cix)
        else             ok = _eq_model(thn, W, LP, LX, Z, a0, qd, D, Un, Gn, 0)
        if (ok) Qn = _eq_obj(Un, om, S)
        else {
            Qn = .
            hb = 1
        }
        if (ok & Qn <= Q) break
        t = t / 2
    } while (t > 1e-12)
    tt = t
    if (t <= 1e-12) return(0)
    r  = mreldif(thn, th)
    th = thn
    if (args() > 15) (void) _eq_model(th, W, LP, LX, Z, a0, qd, D, U, G, 1, PH, PHI, cix)
    else             (void) _eq_model(th, W, LP, LX, Z, a0, qd, D, U, G, 1)
    return(r)
}

// Data diagnostics before estimation (notes returned in the local datanotes,
// separated by |; the VIF in the matrix vifn; the condition number in condx):
//  - near collinearity of the regressors ln(p_k/p_M), ln x and z: condition
//    number of their weighted correlation matrix (Belsley, Kuh and Welsch
//    1980: above 30, harmful) and variance inflation factors
//  - zero shares by good: FGNLS does not model censoring
//  - extreme log prices: more than 5 robust standard deviations (1.4826 MAD)
//    from the median
void _eq_datadiag(string scalar wv, string scalar lpv, string scalar lxv,
                  string scalar zv, string scalar omv, string scalar touse,
                  string scalar vifn, | string scalar selgs)
{
    real rowvector sgf
    real matrix    W, LP, Z, X, Xc, C, R
    real colvector om, w, x, ev, vif, med, mad
    real rowvector mu, sd
    real scalar    N, M, K, j, condx, pz, nout
    string scalar  notes, lst
    string rowvector wn, zn

    W  = st_data(., tokens(wv), touse)
    LP = st_data(., tokens(lpv), touse)
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
    X = (LP[|1, 1 \ N, M - 1|] :- LP[., M]), st_data(., lxv, touse), Z
    w = om :/ sum(om)
    mu = colsum(w :* X)
    Xc = X :- mu
    C  = cross(Xc, w, Xc)
    sd = sqrt(diagonal(C))'
    notes = ""
    if (min(sd) > 0) {
        R = C :/ (sd' * sd)
        ev = symeigenvalues(R)'
        condx = sqrt(max(ev) / max((min(ev), 1e-300)))
        vif = diagonal(invsym(R))
        st_matrix(vifn, vif')
        if (condx > 30) {
            lst = ""
            for (j = 1; j <= cols(X); j++) {
                if (vif[j] > 10) {
                    if (j < M)            lst = lst + sprintf(" ln(p%g/p%g) %4.0f", j, M, vif[j])
                    else if (j == M)      lst = lst + sprintf(" ln_x %4.0f", vif[j])
                    else                  lst = lst + sprintf(" %s %4.0f", zn[j - M], vif[j])
                }
            }
            notes = notes + sprintf("near collinearity among the regressors: condition number %4.0f > 30 (Belsley, Kuh and Welsch 1980); VIF above 10:%s|", condx, lst)
        }
    }
    else {
        condx = .
        st_matrix(vifn, J(1, cols(X), .))
    }
    // weak price variation: the gamma of good k are identified by the
    // variation of ln(p_k/p_M); their variance is inversely proportional to
    // N var(ln(p_k/p_M)).  Standard deviation below 0.05 (prices varying by
    // less than about 5 percent across households)
    lst = ""
    for (j = 1; j < M; j++) {
        if (sd[j] < 0.05) lst = lst + sprintf(" ln(p%g/p%g) %6.4f", j, M, sd[j])
    }
    if (lst != "") notes = notes + sprintf("weak variation of relative prices (standard deviation below 0.05):%s; the price parameters of these goods are weakly identified|", lst)
    // zero shares by good
    // (with selection, the goods corrected are not flagged)
    if (args() > 7 & selgs != "") sgf = strtoreal(tokens(selgs)), 0
    else                          sgf = J(1, M, 0)
    for (j = 1; j <= M; j++) {
        pz = sum(w :* (W[., j] :<= 0))
        if (pz > 0.05 & !sgf[j]) notes = notes + sprintf("%s: %4.1f%% of the shares are zero; FGNLS does not model censored demand (see selection)|", wn[j], 100 * pz)
    }
    // extreme log prices
    nout = 0
    for (j = 1; j <= M; j++) {
        x = LP[., j]
        med = mm_median_eq(x)
        mad = 1.4826 * mm_median_eq(abs(x :- med))
        if (mad > 0) nout = nout + sum(abs(x :- med) :> 5 * mad)
    }
    if (nout > 0) notes = notes + sprintf("%g log prices lie more than 5 robust standard deviations from their median|", nout)
    st_local("datanotes", notes)
    st_local("condx", strofreal(condx, "%21.0g"))
}

real scalar mm_median_eq(real colvector x)
{
    real colvector s
    real scalar    n
    s = sort(x, 1)
    n = rows(s)
    if (mod(n, 2)) return(s[(n + 1) / 2])
    return((s[n / 2] + s[n / 2 + 1]) / 2)
}


// ===========================================================================
// Selection of the buyers (option selection): Shonkwiler and Yen (1999)
//
// For each corrected good i, a probit of buying on s = (1, ln(p_1/x), ...,
// ln(p_M/x), z, the probit-only variables of the good): homogeneous of degree
// zero in prices and expenditure. Its fitted Phi_i, phi_i enter the system
// as  E[w_i] = Phi_i f_i + delta_i phi_i  (the last good, never zero, closes
// the system: E[w_M] = 1 - sum of the others). The delta_i are appended to
// the free parameters.
// ===========================================================================

// weighted probit of dd (0/1) on S by Fisher scoring (concave log
// likelihood); Phi and phi at the estimate, IF the household influence
// functions of the coefficients (w r s') I^-1. ok = 0 when everybody or
// nobody buys.
real colvector _eq_probit(real colvector dd, real matrix S, real colvector w,
                          real colvector Ph, real colvector ph, real matrix IF,
                          real scalar ok)
{
    real colvector g, xb, P, f, r, a, step
    real matrix    Ii
    real scalar    it, n, n1
    n = rows(dd)
    g = J(cols(S), 1, 0)
    n1 = sum(dd)
    ok = (n1 > 0 & n1 < n)
    if (!ok) {
        Ph = J(n, 1, 1) ; ph = J(n, 1, 0) ; IF = J(n, cols(S), 0)
        return(g)
    }
    for (it = 1; it <= 200; it++) {
        xb = S * g
        P  = rowmin((rowmax((normal(xb), J(n, 1, 1e-15))), J(n, 1, 1 - 1e-15)))
        f  = normalden(xb)
        r  = (dd - P) :* f :/ (P :* (1 :- P))
        a  = (f :^ 2) :/ (P :* (1 :- P))
        Ii = invsym(quadcross(S, w :* a, S))
        step = Ii * quadcross(S, w :* r)
        g = g + step
        if (max(abs(step)) < 1e-11) break
    }
    xb = S * g
    P  = rowmin((rowmax((normal(xb), J(n, 1, 1e-15))), J(n, 1, 1 - 1e-15)))
    f  = normalden(xb)
    r  = (dd - P) :* f :/ (P :* (1 :- P))
    a  = (f :^ 2) :/ (P :* (1 :- P))
    Ii = invsym(quadcross(S, w :* a, S))
    Ph = normal(xb)
    ph = f
    IF = ((w :* r) :* S) * Ii
    return(g)
}

// the probit design of good i at data rows: (1, ln p - ln x, z, q[mask_i])
real matrix _eq_seldesign(real matrix LP, real colvector LX, real matrix Z,
                          real matrix Q, real rowvector mask)
{
    real matrix S
    S = (J(rows(LP), 1, 1), LP :- LX, Z)
    if (cols(Q)) {
        if (any(mask)) S = S, Q[., selectindex(mask)]
    }
    return(S)
}

// expected shares and their derivatives at data rows, all M goods:
//   EW N x M, DX N x M (d/d ln x), DP N x M^2 (d E[w_i] / d ln p_j at (i-1)M+j)
// psi = (theta free, delta, alpha of each corrected good in order)
void _eq_selparts(real rowvector psi, real matrix LP, real colvector LX,
                  real matrix Z, real matrix Q, real scalar a0, real scalar qd,
                  real scalar M, real scalar P, real rowvector cix,
                  real matrix zmask, real matrix EW, real matrix DX,
                  real matrix DP)
{
    real matrix    F, MU, MUJ, S
    real colvector xb, Ph, ph, k
    real rowvector th, al
    real scalar    N, C, i, j, c, pos, p, dl, ax
    N = rows(LP)
    C = sum(cix :> 0)
    th = psi[|1 \ P|]
    _eq_parts(th, LP, LX, Z, a0, qd, M, F, MU, MUJ)
    EW = F ; DX = MU ; DP = MUJ
    pos = P + C
    for (i = 1; i < M; i++) {
        if (!cix[i]) continue
        S  = _eq_seldesign(LP, LX, Z, Q, zmask[i, .])
        p  = cols(S)
        al = psi[|pos + 1 \ pos + p|]
        pos = pos + p
        dl = psi[P + cix[i]]
        xb = S * al'
        Ph = normal(xb)
        ph = normalden(xb)
        // d s'alpha / d ln p_j = alpha_pj ; d / d ln x = - sum_j alpha_pj
        k  = ph :* (F[., i] - dl :* xb)
        ax = -sum(al[|2 \ M + 1|])
        EW[., i] = Ph :* F[., i] + dl :* ph
        DX[., i] = Ph :* MU[., i] + ax :* k
        for (j = 1; j <= M; j++) {
            c = (i - 1) * M + j
            DP[., c] = Ph :* MUJ[., c] + al[1 + j] :* k
        }
    }
    // the last good closes the system
    EW[., M] = 1 :- rowsum(EW[|1, 1 \ N, M - 1|])
    DX[., M] = -rowsum(DX[|1, 1 \ N, M - 1|])
    for (j = 1; j <= M; j++) {
        c = (M - 1) * M + j
        DP[., c] = J(N, 1, 0)
        for (i = 1; i < M; i++) DP[., c] = DP[., c] - DP[., (i - 1) * M + j]
    }
}

// the elasticities stacked (E_x (M), E_u (M^2), E_c (M^2)) from the pieces:
//   market: weights omega x on the expected expenditure (T/S)
//   one household (a row): its elasticities
real rowvector _eq_selstack_mkt(real matrix EW, real matrix DX, real matrix DP,
                                real colvector wx, real scalar M)
{
    real rowvector S, ex, eu, ec
    real scalar    i, j, c
    S  = colsum(wx :* EW)
    ex = 1 :+ colsum(wx :* DX) :/ S
    eu = J(1, M * M, .) ; ec = J(1, M * M, .)
    for (i = 1; i <= M; i++) {
        for (j = 1; j <= M; j++) {
            c = (i - 1) * M + j
            eu[c] = -(i == j) + sum(wx :* DP[., c]) / S[i]
            ec[c] = eu[c] + sum(wx :* (EW[., i] + DX[., i]) :* EW[., j]) / S[i]
        }
    }
    return((ex, eu, ec))
}

real matrix _eq_selstack_hh(real matrix EW, real matrix DX, real matrix DP,
                            real scalar M)
{
    real matrix H
    real scalar i, j, c
    H = J(rows(EW), M + 2 * M * M, .)
    H[., 1..M] = 1 :+ DX :/ EW
    for (i = 1; i <= M; i++) {
        for (j = 1; j <= M; j++) {
            c = (i - 1) * M + j
            H[., M + c] = -(i == j) :+ DP[., c] :/ EW[., i]
            H[., M + M * M + c] = H[., M + c] + H[., i] :* EW[., j]
        }
    }
    return(H)
}

// the three summaries at psi: market (row 1), at the means (row 2) and the
// household mean (row 3), each 1 x (M + 2 M^2)
real matrix _eq_selsum(real rowvector psi, real matrix LP, real colvector LX,
                       real matrix Z, real matrix Q, real colvector om,
                       real scalar a0, real scalar qd, real scalar M,
                       real scalar P, real rowvector cix, real matrix zmask,
                       real rowvector m)
{
    real matrix    EW, DX, DP, EWm, DXm, DPm
    real scalar    K, nq, sw
    real rowvector out1, out2, out3
    K = cols(Z) ; nq = cols(Q) ; sw = sum(om)
    _eq_selparts(psi, LP, LX, Z, Q, a0, qd, M, P, cix, zmask, EW, DX, DP)
    out1 = _eq_selstack_mkt(EW, DX, DP, om :* exp(LX), M)
    out3 = colsum(om :* _eq_selstack_hh(EW, DX, DP, M)) :/ sw
    _eq_selparts(psi, m[|1 \ M|], m[M + 1],
                 (K ? m[|M + 2 \ M + 1 + K|] : J(1, 0, .)),
                 (nq ? m[|M + 2 + K \ M + 1 + K + nq|] : J(1, 0, .)),
                 a0, qd, M, P, cix, zmask, EWm, DXm, DPm)
    out2 = _eq_selstack_hh(EWm, DXm, DPm, M)
    return((out1 \ out2 \ out3))
}

// estimation of the probits, before the system
void _eq_selprobits(real matrix W, real matrix LP, real colvector LX,
                    real matrix Z, real matrix Q, real colvector om,
                    real rowvector sg, real matrix zmask, real matrix PH,
                    real matrix PHI, real rowvector cix, real rowvector alpha,
                    real matrix IFA, string scalar diagn, string scalar alphan)
{
    real matrix    S, IFg, AL, DG
    real colvector dd, Ph, ph, g, xb, pr
    real scalar    N, M, i, C, okp, pb, pz, ll, ll0, pbar
    N = rows(W) ; M = cols(W)
    PH = J(N, M - 1, 1) ; PHI = J(N, M - 1, 0) ; cix = J(1, M - 1, 0)
    alpha = J(1, 0, .) ; IFA = J(N, 0, .)
    pb = 1 + M + cols(Z) ; pz = pb + cols(Q)
    AL = J(pz, M - 1, .)
    DG = J(M - 1, 4, .)
    C = 0
    for (i = 1; i < M; i++) {
        if (!sg[i]) continue
        C++
        cix[i] = C
        S  = _eq_seldesign(LP, LX, Z, Q, zmask[i, .])
        dd = W[., i] :> 0
        g  = _eq_probit(dd, S, om, Ph, ph, IFg, okp)
        if (!okp) {
            errprintf("equaids: selection: every household (or none) buys good %g\n", i)
            exit(459)
        }
        PH[., i] = Ph ; PHI[., i] = ph
        alpha = alpha, g'
        IFA = IFA, IFg
        if (cols(Q)) {
            if (any(zmask[i, .])) AL[(1::pb) \ (pb :+ selectindex(zmask[i, .])'), i] = g
            else AL[1::pb, i] = g
        }
        else AL[1::pb, i] = g
        // diagnostics: buyers, McFadden pseudo-R2, households predicted with
        // probability 0 or 1 (separation)
        pbar = quadsum(om :* dd) / quadsum(om)
        xb = S * g
        ll  = quadsum(om :* (dd :* lnnormal(xb) + (1 :- dd) :* lnnormal(-xb)))
        ll0 = quadsum(om) * (pbar * ln(pbar) + (1 - pbar) * ln(1 - pbar))
        pr = normal(xb)
        DG[i, 1] = 100 * pbar
        DG[i, 2] = 1 - ll / ll0
        DG[i, 3] = sum((pr :< 1e-10) :| (pr :> 1 - 1e-10))
    }
    st_matrix(diagn, DG)
    st_matrix(alphan, AL)
}

// after the system: the influence functions of theta with the estimation of
// the probits, the variances of delta, the elasticities of the four types on
// the expected shares and their variances (derivatives with respect to all
// the parameters by central differences), the identification of delta
void _eq_selpost(real rowvector th, real matrix W, real matrix LP,
                 real colvector LX, real matrix Z, real matrix Q,
                 real colvector om, real scalar a0, real scalar qd,
                 real scalar M, real scalar P, real matrix D, real matrix U,
                 pointer(real matrix) rowvector G, real matrix S,
                 real matrix A, real matrix PH, real matrix PHI,
                 real rowvector cix, real matrix zmask, real rowvector alpha,
                 real matrix IFA, string scalar base, real scalar doe,
                 real matrix DS, real scalar vmode, string scalar single,
                 string scalar diagn)
{
    real matrix    Ai, US, IFt, H, Sd, F, MU, MUJ, Cj, IFp, JJ, Jm, IFs, MV,
                   EW, DX, DP, V, Hh, Gi, DG, Gx, EWm, DXm, DPm
    real rowvector psi, pp, pm, m, mp, mm, E0, dl, sdl, S0
    real colvector wx, fsel, xb, rr
    real scalar    N, C, Pt, q, i, j, k, h, pos, p, nm, sw, R2, K, nq, c, r
    N = rows(W) ; C = sum(cix :> 0) ; Pt = cols(th) ; q = M + 2 * M * M
    K = cols(Z) ; nq = cols(Q) ; sw = sum(om)
    // influence functions of theta, with the probits: IF = A^-1 (s - C IF_alpha)
    Ai  = invsym(A)
    US  = U * S
    IFt = J(N, Pt, 0)
    for (i = 1; i < M; i++) IFt = IFt + (*G[i]) :* (om :* US[., i])
    _eq_parts(th[|1 \ P|], LP, LX, Z, a0, qd, M, F, MU, MUJ)
    pos = 0
    for (j = 1; j < M; j++) {
        if (!cix[j]) continue
        Sd = _eq_seldesign(LP, LX, Z, Q, zmask[j, .])
        p  = cols(Sd)
        H  = J(N, Pt, 0)
        for (i = 1; i < M; i++) H = H + S[i, j] :* *G[i]
        xb = Sd * alpha[|pos + 1 \ pos + p|]'
        fsel = PHI[., j] :* (F[., j] - th[P + cix[j]] :* xb)
        Cj = cross(H, om, fsel :* Sd)
        IFt = IFt - IFA[|1, pos + 1 \ N, pos + p|] * Cj'
        pos = pos + p
    }
    IFt = IFt * Ai
    st_matrix(base + "Vt", _eq_vagg(IFt, DS, vmode, single))
    _eq_saveif(IFt, "t")
    // delta and its standard error
    V = _eq_vagg(IFt, DS, vmode, single)
    dl = J(1, M - 1, .) ; sdl = J(1, M - 1, .)
    for (i = 1; i < M; i++) {
        if (!cix[i]) continue
        dl[i]  = th[P + cix[i]]
        sdl[i] = sqrt(V[P + cix[i], P + cix[i]])
    }
    st_matrix(base + "dl", dl)
    st_matrix(base + "sdl", sdl)
    // what the correction costs the expenditure coefficient beta_i: the
    // variance inflation 1/(1 - rho^2), rho the correlation of phi_i and of
    // the gradient column of beta_i given the other columns of the gradient
    // of equation i (as the VIF of the quality elasticity in duvm)
    DG = st_matrix(diagn)
    for (i = 1; i < M; i++) {
        if (!cix[i]) continue
        Gi = *G[i]
        Gx = Gi
        Gx[., P + cix[i]] = J(N, 1, 0)
        Gx[., M - 1 + i]  = J(N, 1, 0)
        Gx = select(Gx, colsum(abs(Gx)) :> 0)
        Gx = J(N, 1, 1), Gx
        Hh = invsym(cross(Gx, om, Gx))
        rr = PHI[., i] - Gx * (Hh * cross(Gx, om, PHI[., i]))
        wx = Gi[., M - 1 + i] - Gx * (Hh * cross(Gx, om, Gi[., M - 1 + i]))
        R2 = cross(rr, om, wx)^2 / (cross(rr, om, rr) * cross(wx, om, wx))
        DG[i, 4] = 1 / (1 - R2)
    }
    st_matrix(diagn, DG)
    // the elasticities and the shares at which they are evaluated
    psi = th, alpha
    MV  = (LP, LX, Z, Q)
    nm  = cols(MV)
    m   = colsum(om :* MV) :/ sw
    E0  = _eq_selsum(psi, LP, LX, Z, Q, om, a0, qd, M, P, cix, zmask, m)
    _eq_selparts(psi, LP, LX, Z, Q, a0, qd, M, P, cix, zmask, EW, DX, DP)
    wx = om :* exp(LX)
    S0 = colsum(wx :* EW)
    st_matrix(base + "S",  S0 :/ sum(S0))
    st_matrix(base + "Sh", colsum(om :* EW) :/ sw)
    st_matrix(base + "ns", colsum(EW :< 0.001))
    st_matrix(base + "x",  E0[1, 1..M])
    st_matrix(base + "u",  rowshape(E0[1, M+1..M+M*M], M))
    st_matrix(base + "c",  rowshape(E0[1, M+M*M+1..q], M))
    st_matrix(base + "xm", E0[2, 1..M])
    st_matrix(base + "um", rowshape(E0[2, M+1..M+M*M], M))
    st_matrix(base + "cm", rowshape(E0[2, M+M*M+1..q], M))
    st_matrix(base + "xh", E0[3, 1..M])
    st_matrix(base + "uh", rowshape(E0[3, M+1..M+M*M], M))
    st_matrix(base + "ch", rowshape(E0[3, M+M*M+1..q], M))
    _eq_selparts(psi, m[|1 \ M|], m[M + 1], (K ? m[|M + 2 \ M + 1 + K|] : J(1, 0, .)),
                 (nq ? m[|M + 2 + K \ nm|] : J(1, 0, .)), a0, qd, M, P, cix, zmask, EWm, DXm, DPm)
    st_matrix(base + "Sm", EWm)
    // psi = (theta, alpha) and its variance: the Engel curves of estat
    IFp = IFt, IFA
    st_matrix(base + "psi", psi)
    st_matrix(base + "Vpsi", _eq_vagg(IFp, DS, vmode, single))
    if (!doe) return
    // derivatives of the three summaries with respect to psi and to the means
    JJ = J(3 * q, cols(psi), .)
    for (k = 1; k <= cols(psi); k++) {
        h = 1e-5 * max((1, abs(psi[k])))
        pp = psi ; pp[k] = psi[k] + h
        pm = psi ; pm[k] = psi[k] - h
        JJ[., k] = vec((_eq_selsum(pp, LP, LX, Z, Q, om, a0, qd, M, P, cix, zmask, m) -
                       _eq_selsum(pm, LP, LX, Z, Q, om, a0, qd, M, P, cix, zmask, m))')  :/ (2 * h)
    }
    Jm = J(q, nm, .)
    for (k = 1; k <= nm; k++) {
        h = 1e-5 * max((1, abs(m[k])))
        mp = m ; mp[k] = m[k] + h
        mm = m ; mm[k] = m[k] - h
        Jm[., k] = ((_eq_selsum(psi, LP, LX, Z, Q, om, a0, qd, M, P, cix, zmask, mp)[2, .] -
                     _eq_selsum(psi, LP, LX, Z, Q, om, a0, qd, M, P, cix, zmask, mm)[2, .]) :/ (2 * h))'
    }
    // market: sampling term of the ratios of totals
    IFs = J(N, q, .)
    for (i = 1; i <= M; i++) {
        IFs[., i] = wx :* (DX[., i] - (E0[1, i] - 1) :* EW[., i]) :/ S0[i]
        for (j = 1; j <= M; j++) {
            c = (i - 1) * M + j
            IFs[., M + c] = wx :* (DP[., c] - (E0[1, M + c] + (i == j)) :* EW[., i]) :/ S0[i]
            r = sum(wx :* (EW[., i] + DX[., i]) :* EW[., j]) / S0[i]
            IFs[., M + M * M + c] = IFs[., M + c] +
                wx :* ((EW[., i] + DX[., i]) :* EW[., j] - r :* EW[., i]) :/ S0[i]
        }
    }
    IFs = IFs + IFp * JJ[1..q, .]'
    V = _eq_vagg(IFs, DS, vmode, single)
    _eq_saveif(IFs[., 1..M], "x")
    _eq_saveif(IFs[., M+1..M+M*M], "u")
    st_matrix(base + "Vx", V[1..M, 1..M])
    st_matrix(base + "Vu", V[M+1..M+M*M, M+1..M+M*M])
    st_matrix(base + "Vc", V[M+M*M+1..q, M+M*M+1..q])
    st_matrix(base + "Jx", JJ[1..M, 1..P])
    st_matrix(base + "Ju", JJ[M+1..M+M*M, 1..P])
    st_matrix(base + "Jc", JJ[M+M*M+1..q, 1..P])
    // at the means
    V = _eq_vagg(IFp * JJ[q+1..2*q, .]' + ((om :* (MV :- m)) :/ sw) * Jm', DS, vmode, single)
    st_matrix(base + "Vxm", V[1..M, 1..M])
    st_matrix(base + "Vum", V[M+1..M+M*M, M+1..M+M*M])
    st_matrix(base + "Vcm", V[M+M*M+1..q, M+M*M+1..q])
    // household mean
    Hh = _eq_selstack_hh(EW, DX, DP, M)
    V = _eq_vagg((om :* (Hh :- E0[3, .])) :/ sw + IFp * JJ[2*q+1..3*q, .]', DS, vmode, single)
    st_matrix(base + "Vxh", V[1..M, 1..M])
    st_matrix(base + "Vuh", V[M+1..M+M*M, M+1..M+M*M])
    st_matrix(base + "Vch", V[M+M*M+1..q, M+M*M+1..q])
}

// the smallest m0(z) / se(m0(z)) over the sample: the distance of the estimate
// to the boundary of Ray's scaling in standard errors (local m0t)
void _eq_m0dist(string scalar zv, string scalar touse)
{
    string matrix  st
    real colvector k, m0, se
    real rowvector b, rh
    real matrix    Z, V
    st = st_matrixcolstripe("e(b)")
    k  = selectindex(st[., 1] :== "rho")
    b  = st_matrix("e(b)")
    V  = st_matrix("e(V)")
    rh = b[k']
    Z  = st_data(., tokens(zv), touse)
    m0 = 1 :+ Z * rh'
    se = sqrt(rowsum((Z * V[k, k]) :* Z))
    st_local("m0t", strofreal(min(m0 :/ se), "%21.0g"))
}

// one bootstrap replication: e(b), e(b_free), the nine elasticity matrices
// (row by row), and under selection delta (0 where not corrected) and psi
real rowvector _eq_bvec(real scalar sel)
{
    real rowvector v
    string rowvector nm
    real scalar i
    nm = ("x", "u", "c", "xm", "um", "cm", "xh", "uh", "ch")
    v = st_matrix("e(b)"), st_matrix("e(b_free)")
    for (i = 1; i <= 9; i++) v = v, vec(st_matrix("e(elas_" + nm[i] + ")")')'
    if (sel) v = v, editmissing(st_matrix("e(sel_delta)"), 0), st_matrix("e(sel_psi)")
    return(v)
}

// the variances of the replications, in the layout of _eq_bvec
void _eq_bpost(real matrix B, real scalar sel, string scalar base)
{
    real matrix    C
    real rowvector sd, dl
    string rowvector nm
    real scalar    nb, nf, M, pos, i, k, np
    C  = variance(B)
    nb = cols(st_matrix("e(b)"))
    nf = cols(st_matrix("e(b_free)"))
    M  = st_numscalar("e(ngoods)")
    nm = ("x", "u", "c", "xm", "um", "cm", "xh", "uh", "ch")
    st_matrix(base + "Vb", C[|1, 1 \ nb, nb|])
    pos = nb
    st_matrix(base + "Vf", C[|pos + 1, pos + 1 \ pos + nf, pos + nf|])
    pos = pos + nf
    for (i = 1; i <= 9; i++) {
        k = (mod(i - 1, 3) == 0 ? M : M * M)
        st_matrix(base + "V" + nm[i], C[|pos + 1, pos + 1 \ pos + k, pos + k|])
        pos = pos + k
    }
    if (sel) {
        dl = st_matrix("e(sel_delta)")
        sd = sqrt(diagonal(C[|pos + 1, pos + 1 \ pos + M - 1, pos + M - 1|]))'
        for (i = 1; i < M; i++) if (missing(dl[i])) sd[i] = .
        st_matrix(base + "sdl", sd)
        pos = pos + M - 1
        np = cols(st_matrix("e(sel_psi)"))
        st_matrix(base + "Vpsi", C[|pos + 1, pos + 1 \ pos + np, pos + np|])
    }
}

void _eq_fit(string scalar wv, string scalar lpv, string scalar lxv,
             string scalar zv, string scalar omv, string scalar touse,
             real scalar a0, real scalar qd, real scalar tol,
             real scalar maxit, real scalar showlog, string scalar fromn,
             string scalar bfn, string scalar Vfn, string scalar bn,
             string scalar Vn, string scalar Sgn, string scalar elbase,
             real scalar doe, string scalar desv, real scalar vmode,
             string scalar single, real scalar sel, string scalar qv,
             string scalar selgs, string scalar zms, string scalar seldiag,
             string scalar selalpha)
{
    real matrix    Q, zmask, PH, PHI, IFA, Dx
    real rowvector sg, cix, alpha
    real scalar    C, Pt
    real matrix    W, LP, Z, D, U, S, Sig, Sigo, A, V, Vf, As, X, DS, H
    real colvector LX, om, ndh, ngh, dA, m0, ell
    real rowvector th, c0, bfull, e, al, be, la, rh
    real matrix    Ga, et, F
    real scalar    N, M, K, P, Pf, i, j, r, ok, outer, inner,
                   it, conv, Nw, ll, sw, nd, sd, tt, hb, nhb, stall,
                   rc, m0min, nlneg, nneg
    pointer(real matrix) rowvector G

    W  = st_data(., tokens(wv), touse)
    LP = st_data(., tokens(lpv), touse)
    LX = st_data(., lxv, touse)
    om = st_data(., omv, touse)
    if (zv != "") Z = st_data(., tokens(zv), touse)
    else          Z = J(rows(W), 0, .)
    N = rows(W) ; M = cols(W) ; K = cols(Z)
    P = 2 * (M - 1) + (M - 1) * M / 2 + qd * (M - 1) + K * (M - 1) + K

    // ---- identification, before estimating ----
    // gamma: under homogeneity only relative prices matter, so the log
    // prices relative to good M must vary independently of each other and
    // of the constant
    X = J(N, 1, 1), LP[|1, 1 \ N, M - 1|] :- LP[., M]
    if (rank(X) < M) {
        errprintf("equaids: the relative prices ln(p_k/p_M) do not vary independently;\n")
        errprintf("         the price parameters gamma are not identified\n")
        exit(459)
    }
    // beta (and lambda): log expenditure must vary
    if (rank((J(N, 1, 1), LX)) < 2) {
        errprintf("equaids: log expenditure does not vary\n")
        exit(459)
    }
    // eta, rho: the demographics must not be collinear with the constant
    if (K > 0) {
        if (rank((J(N, 1, 1), Z)) < K + 1) {
            errprintf("equaids: the demographics are collinear (with each other or with a constant)\n")
            exit(459)
        }
    }

    // affine map free -> full: full(th) = c0 + th * D'
    c0 = _eq_full(J(1, P, 0), M, K, qd)
    Pf = cols(c0)
    D  = J(Pf, P, 0)
    for (r = 1; r <= P; r++) {
        e = J(1, P, 0) ; e[r] = 1
        D[., r] = (_eq_full(e, M, K, qd) - c0)'
    }

    // selection: the probits of the goods corrected, before the system
    C = 0
    if (sel) {
        if (qv != "") Q = st_data(., tokens(qv), touse)
        else          Q = J(N, 0, .)
        sg = strtoreal(tokens(selgs))
        if (cols(Q)) zmask = rowshape(strtoreal(tokens(zms)), M - 1)
        else         zmask = J(M - 1, 0, 0)
        _eq_selprobits(W, LP, LX, Z, Q, om, sg, zmask, PH, PHI, cix, alpha, IFA, seldiag, selalpha)
        C = sum(cix :> 0)
    }

    // start: from() if given, else alpha at the weighted mean shares and
    // the rest at zero
    sw = sum(om)
    if (fromn != "") {
        th = st_matrix(fromn)
        if (rows(th) != 1 | (cols(th) != P & cols(th) != P + C)) {
            errprintf("equaids: from() must be a 1 x %g row vector of free parameters (e(b_free))\n", P + C)
            exit(503)
        }
        if (cols(th) == P & C) th = th, J(1, C, 0)
    }
    else {
        th = J(1, P, 0)
        th[|1 \ M - 1|] = (colsum(om :* W[|1, 1 \ N, M - 1|]) :/ sw)
        if (C) th = th, J(1, C, 0)
    }
    // phase 1: NLS (Sigma = I) to a loose tolerance, a start for phase 2
    it = 0
    conv = 0
    if (sel) ok = _eq_model(th, W, LP, LX, Z, a0, qd, D, U, G, 1, PH, PHI, cix)
    else     ok = _eq_model(th, W, LP, LX, Z, a0, qd, D, U, G, 1)
    if (!ok) {
        errprintf("equaids: m0(z) = 1 + rho'z is not positive at the start\n")
        exit(459)
    }
    S = I(M - 1)
    nhb = 0
    for (inner = 1; inner <= maxit; inner++) {
        it++
        if (sel) r = _eq_step(th, S, W, LP, LX, Z, a0, qd, D, om, U, G, nd, tt, hb, PH, PHI, cix)
        else     r = _eq_step(th, S, W, LP, LX, Z, a0, qd, D, om, U, G, nd, tt, hb)
        nhb = nhb + hb
        if (nd < 1e-4) break
    }
    // phase 2: iterated FGNLS, Sigma updated at every Gauss-Newton step;
    // the fixed point is the Gaussian maximum-likelihood estimate.
    // Convergence: the Newton decrement g'A^-1 g below tol (the score is
    // zero up to tol, whatever the units of the parameters) and Sigma
    // stable at the fixed point (relative change below tol)
    // The tolerance has a statistical meaning: g'A^-1 g is the squared
    // distance to the optimum in the metric of the information, so
    // g'A^-1 g < 1e-6 (the default) puts the estimate within 0.001 standard
    // errors of the optimum, in all directions jointly (Stata's ml uses
    // nrtolerance(1e-5)); Sigma stable to the same relative order.
    // Stall: over the last 20 iterations NEITHER measure of progress moved --
    // the log likelihood rose by less than 1e-7 per observation AND the
    // scaled gradient did not fall at all (steps cut to nothing at a
    // boundary, a singular information matrix).  Either kind of progress
    // is not a stall: Gauss-Newton converges only linearly when the residuals
    // are large (it leaves out the residual times second-derivative term of
    // the Hessian), so its predicted gain g'A^-1 g / 2 overstates the realized
    // gain even when it converges, and near the optimum the gradient may fall
    // by only a few percent per iteration (both seen on bootstrap
    // replications); along a curved ridge the likelihood keeps rising while
    // the gradient does not fall (the iteration limit applies).
    Sig = cross(U, om, U) :/ sw
    ndh = J(maxit, 1, .)
    ngh = J(maxit, 1, .)
    stall = 0
    for (outer = 1; outer <= maxit; outer++) {
        it++
        S = invsym(Sig)
        if (sel) r = _eq_step(th, S, W, LP, LX, Z, a0, qd, D, om, U, G, nd, tt, hb, PH, PHI, cix)
        else     r = _eq_step(th, S, W, LP, LX, Z, a0, qd, D, om, U, G, nd, tt, hb)
        nhb = nhb + hb
        Sigo = Sig
        Sig  = cross(U, om, U) :/ sw
        sd   = mreldif(Sig, Sigo)
        ndh[outer] = -sw / 2 * ln(det(Sig))      // log likelihood up to a constant
        if (showlog) printf("{txt}FGNLS iteration %g: scaled gradient %9.2e, Sigma change %9.2e, step %g%s\n",
            outer, nd, sd, tt, (hb ? " (boundary m0 > 0)" : ""))
        if (nd < tol & sd < tol) {
            conv = 1
            break
        }
        ngh[outer] = nd
        if (outer > 20) {
            if (ndh[outer] - ndh[outer - 20] < 1e-7 * sw & nd >= ngh[outer - 20]) {
                stall = 1
                break
            }
        }
    }

    // conventional variance at the final Sigma
    S = invsym(Sig)
    if (sel) (void) _eq_model(th, W, LP, LX, Z, a0, qd, D, U, G, 1, PH, PHI, cix)
    else     (void) _eq_model(th, W, LP, LX, Z, a0, qd, D, U, G, 1)
    // grouped as in _eq_step
    Pt = cols(th)
    A = J(Pt, Pt, 0)
    for (i = 1; i < M; i++) {
        H = J(N, Pt, 0)
        for (j = 1; j < M; j++) H = H + S[i, j] :* *G[j]
        A = A + cross(*G[i], om, H)
    }
    Vf = invsym(A)
    // with selection, the deltas follow the free parameters of the model
    if (C) {
        Dx = (D, J(rows(D), C, 0)) \ (J(C, cols(D), 0), I(C))
        V  = Dx * Vf * Dx'
        bfull = (c0 + th[|1 \ P|] * D'), th[|P + 1 \ Pt|]
    }
    else {
        Dx = D
        V  = D * Vf * D'
        bfull = c0 + th * D'
    }
    Nw = sw
    ll = -Nw / 2 * ((M - 1) * (1 + ln(2 * pi())) + ln(det(Sig)))

    // ---- diagnostics at the estimate ----
    // conditioning of the information matrix, scaled to unit diagonal so
    // that it does not depend on the units of the parameters
    dA = sqrt(diagonal(A))
    As = A :/ (dA * dA')
    rc = 1 / cond(As)
    // deflated expenditure ln(x / (m0 a(p))) and the scale m0(z)
    _eq_unpack(bfull, M, K, qd, al, be, Ga, la, et, rh)
    if (K > 0) m0 = 1 :+ Z * rh'
    else       m0 = J(N, 1, 1)
    m0min = min(m0)
    st_local("nm0low", strofreal(sum(m0 :< 0.05)))
    ell = LX - ln(m0) - (a0 :+ LP * al' :+ 0.5 :* rowsum((LP * Ga) :* LP))
    nlneg = sum(ell :<= 0)
    // predicted shares outside [0, 1]
    F = W[|1, 1 \ N, M - 1|] - U
    F = F, 1 :- rowsum(F)
    nneg = sum(rowmin(F) :< 0 :| rowmax(F) :> 1)

    st_matrix(bfn, th)
    st_matrix(Vfn, Vf)
    st_matrix(bn, bfull)
    st_matrix(Vn, V)
    st_matrix(Sgn, Sig)
    st_local("ll",   strofreal(ll, "%21.0g"))
    st_local("iter", strofreal(it))
    st_local("conv", strofreal(conv))
    st_local("Nw",   strofreal(Nw, "%21.0g"))
    st_local("nd",   strofreal(nd, "%21.0g"))
    st_local("sd",   strofreal(sd, "%21.0g"))
    st_local("stall", strofreal(stall))
    st_local("nhb",   strofreal(nhb))
    st_local("rc",    strofreal(rc, "%21.0g"))
    st_local("m0min", strofreal(m0min, "%21.0g"))
    st_local("nlneg", strofreal(nlneg))
    st_local("nneg",  strofreal(nneg))
    if (vmode > 0) DS = st_data(., tokens(desv), touse)
    else           DS = J(0, 3, .)
    if (!sel) {
        _eq_elas(th, LP, LX, Z, om, a0, qd, M, elbase)
        _eq_elasvar(th, LP, LX, Z, om, a0, qd, M, D, U, G, S, A, elbase, doe, DS, vmode, single)
    }
    else _eq_selpost(th, W, LP, LX, Z, Q, om, a0, qd, M, P, D, U, G, S, A, PH, PHI,
                     cix, zmask, alpha, IFA, elbase, doe, DS, vmode, single, seldiag)
    st_matrix(elbase + "D", Dx)
}

// aggregation identities of the aggregate elasticities (exact by
// construction; a failure would mean a wrong formula): Engel, Cournot,
// homogeneity, symmetry of the compensated (Slutsky) matrix
void _eq_checks(string scalar base)
{
    real rowvector S, Ex
    real matrix    Eu, Ec
    S  = st_matrix(base + "S")
    Ex = st_matrix(base + "x")
    Eu = st_matrix(base + "u")
    Ec = st_matrix(base + "c")
    st_local("chk_engel",   strofreal(abs(S * Ex' - 1), "%21.0g"))
    st_local("chk_cournot", strofreal(max(abs(S * Eu + S)), "%21.0g"))
    st_local("chk_homog",   strofreal(max(abs(rowsum(Eu) + Ex')), "%21.0g"))
    st_local("chk_symm",    strofreal(max(abs(S' :* Ec - (S' :* Ec)')), "%21.0g"))
}

// standard errors from a variance matrix, shaped r x c (row-major)
void _eq_semat(string scalar vn, real scalar r, string scalar out)
{
    st_matrix(out, rowshape(sqrt(diagonal(st_matrix(vn)))', r))
}

end

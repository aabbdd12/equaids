*! equaids_estat 1.0.0  2026-09-25  Abdelkrim Araar
*! estat after equaids:
*!   estat diagnostics   the diagnostics of the estimate and of the data
*!   estat engel         the Engel curves (layout of easi and duvm)
program define equaids_estat, rclass
    version 14.2
    if "`e(cmd)'" != "equaids" error 301
    gettoken sub rest : 0, parse(" ,")
    local lsub = length("`sub'")
    if `lsub' >= 4 & "`sub'" == substr("diagnostics", 1, `lsub') {
        _equaids_estat_diag `rest'
        return add
    }
    else if "`sub'" == "engel" {
        _equaids_estat_engel `rest'
        return add
    }
    else {
        di as err "estat `sub' not allowed after equaids; use estat diagnostics or estat engel"
        exit 321
    }
end

* ---------------------------------------------------------------------------
* estat engel: the Engel curve of every good, one panel per good, against the
* percentiles of total expenditure (the layout of easi's and duvm's estat
* engel).  At the means (default): the model's share over a grid of ln x,
* log prices and demographics at their weighted means, with a delta-method
* band and, for QUAIDS, the turning point.  As observed: the fitted shares of
* the households, smoothed; observed adds the smoothed observed shares.
program define _equaids_estat_engel, rclass
    syntax [if] [in] [, ATMeans ASObserved OBServed N(integer 100) TRIM(real 1) ///
        Level(cilevel) NOCI DATA(string) SAVing(string asis) NODRAW          ///
        BWidth(real 0) LNX NOTURN *]
    if "`atmeans'" != "" & "`asobserved'" != "" {
        di as err "specify atmeans or asobserved, not both"
        exit 198
    }
    local asobs = ("`asobserved'" != "")
    if `bwidth' != 0 & !`asobs' {
        di as err "bwidth() goes with asobserved: the Engel curve at the means is exact, there is nothing to smooth"
        exit 198
    }
    if "`observed'" != "" & !`asobs' {
        di as err "observed goes with asobserved: the observed shares vary with prices and"
        di as err "demographics as the fitted shares as observed do, not as the curve at the means"
        exit 198
    }
    if `n' < 2 {
        di as err "n() must be at least 2"
        exit 198
    }
    if `trim' < 0 | `trim' >= 50 {
        di as err "trim() must be in [0, 50)"
        exit 198
    }
    local shares `e(lhs)'
    local M = e(ngoods)
    local K = e(ndemos)
    local names `shares'
    if "`e(snames)'" != "" local names `e(snames)'
    local quad = ("`e(model)'" == "QUAIDS")

    marksample touse, novarlist
    qui replace `touse' = 0 if !e(sample)
    tempvar wt lnxv
    local wexp = trim(subinstr(`"`e(wexp)'"', "=", "", 1))
    if "`wexp'" == "" qui gen double `wt' = 1
    else              qui gen double `wt' = `wexp'
    if "`e(expenditure)'" != "" qui gen double `lnxv' = ln(`e(expenditure)')
    else                        qui gen double `lnxv' = `e(lnexpenditure)'
    local lnp
    local k 0
    foreach v in `e(prices)'`e(lnprices)' {
        local ++k
        tempvar lp`k'
        if "`e(prices)'" != "" qui gen double `lp`k'' = ln(`v')
        else                   qui gen double `lp`k'' = `v'
        local lnp `lnp' `lp`k''
    }
    local demos `e(demographics)'

    * the band: the variance of the estimate (robust, cluster or design), t
    * with the design degrees of freedom under vce(svy)
    local ci = ("`noci'" == "" & !`asobs')
    if `ci' {
        tempname VF
        matrix `VF' = e(V_free)
        if matmissing(`VF') {
            di as txt "(no confidence band: the variance of the estimate is missing)"
            local ci 0
        }
    }
    if e(df_r) < . local zc = invttail(e(df_r), (100 - `level') / 200)
    else           local zc = invnormal((100 + `level') / 200)

    * weighted means of the log prices and demographics (as elasticities(means))
    tempname LPM ZM
    matrix `LPM' = J(1, `M', .)
    local k 0
    foreach v of local lnp {
        local ++k
        qui summarize `v' [aw=`wt'] if `touse', meanonly
        matrix `LPM'[1, `k'] = r(mean)
    }
    local zopt
    if `K' > 0 {
        matrix `ZM' = J(1, `K', .)
        local k 0
        foreach v of local demos {
            local ++k
            qui summarize `v' [aw=`wt'] if `touse', meanonly
            matrix `ZM'[1, `k'] = r(mean)
        }
        local zopt zm(`ZM')
    }

    preserve
    qui keep if `touse'
    if `asobs' {
        local zv
        if `K' > 0 local zv z(`demos')
        quietly equaids, _engel mode(obs) lp(`lnp') lx(`lnxv') `zv' touse(`touse') out(_eqr)
    }
    * only what the curves need (the data may hold a pctile or an lnexp)
    local need `lnxv' `wt' `touse' `lnp' `demos' `shares'
    if `asobs' local need `need' _eqrw*
    keep `need'
    local nobs = _N
    if `n' > `nobs' local n = `nobs'
    tempvar gx gp gt
    qui gen double `gx' = .
    qui gen double `gp' = .
    forvalues i = 1/`n' {
        local q = `trim' + (100 - 2 * `trim') * (`i' - 0.5) / `n'
        qui _pctile `lnxv' [aw=`wt'], percentiles(`q')
        qui replace `gx' = r(r1) in `i'
        qui replace `gp' = `q' in `i'
    }
    qui gen byte `gt' = (_n <= `n')
    tempname TURN
    matrix `TURN' = J(2, `M', .)
    matrix colnames `TURN' = `names'
    matrix rownames `TURN' = lnx pctile
    if !`asobs' {
        local nv = cond(`ci', "", "novar")
        quietly equaids, _engel mode(grid) lx(`gx') touse(`gt') out(_eqg) ///
            lpm(`LPM') `zopt' `nv'
        if !r(ok) {
            di as err "m0(z) <= 0 at the means of the demographics: no Engel curve"
            exit 459
        }
        tempname LT
        matrix `LT' = r(lnx_turn)
        * the turning point on the percentile scale of ln x
        qui summarize `wt', meanonly
        local wsum = r(sum)
        forvalues j = 1/`M' {
            local t = `LT'[1, `j']
            if `t' < . {
                matrix `TURN'[1, `j'] = `t'
                qui summarize `wt' if `lnxv' <= `t', meanonly
                matrix `TURN'[2, `j'] = 100 * cond(r(N), r(sum), 0) / `wsum'
            }
        }
        forvalues j = 1/`M' {
            qui gen double _w`j' = _eqgw`j'
            if `ci' qui gen double _se`j' = _eqgse`j'
        }
        local what "at the means of log prices and demographics"
    }
    else {
        if `bwidth' == 0 {
            qui lpoly _eqrw1 `lnxv' [aw=`wt'], degree(1) nograph
            local bwidth = r(bwidth)
        }
        forvalues j = 1/`M' {
            qui lpoly _eqrw`j' `lnxv' [aw=`wt'], degree(1) bwidth(`bwidth') at(`gx') nograph generate(_w`j')
            if "`observed'" != "" {
                qui lpoly `: word `j' of `shares'' `lnxv' [aw=`wt'], degree(1) bwidth(`bwidth') ///
                    at(`gx') nograph generate(_o`j')
            }
        }
        local bws = string(`bwidth', "%6.4f")
        local what "fitted shares as observed, local linear, bandwidth `bws'"
    }
    qui keep in 1/`n'
    qui gen double pctile = `gp'
    qui gen double lnexp = `gx'
    forvalues j = 1/`M' {
        if `ci' {
            qui gen double _lo`j' = _w`j' - `zc' * _se`j'
            qui gen double _hi`j' = _w`j' + `zc' * _se`j'
        }
        label variable _w`j' "`: word `j' of `names''"
        if "`observed'" != "" label variable _o`j' "`: word `j' of `names'', observed"
    }
    label variable pctile "Percentiles of total expenditure"
    label variable lnexp "Log of total expenditure"
    local keepv pctile lnexp _w*
    if `ci' local keepv `keepv' _se* _lo* _hi*
    if "`observed'" != "" local keepv `keepv' _o*
    keep `keepv'
    order pctile lnexp

    if `"`data'"' != "" {
        * data(filename) or data(filename, replace): the file is replaced
        gettoken dfile : data, parse(",")
        local dfile = trim(`"`dfile'"')
        qui save `"`dfile'"', replace
        di as txt `"curve data saved to {bf:`dfile'}"'
    }

    if "`nodraw'" == "" {
        if "`lnx'" != "" {
            local xv lnexp
            local xl "xlabel(, labsize(vsmall))"
            local xt "Log of total expenditure"
            local tr 1
        }
        else {
            local xv pctile
            local xl "xlabel(0(20)100, labsize(vsmall))"
            local xt "Percentiles of total expenditure"
            local tr 2
        }
        qui summarize `xv', meanonly
        local xmin = r(min)
        local xmax = r(max)
        local plots
        local anyturn 0
        forvalues j = 1/`M' {
            local band
            if `ci' local band (rarea _lo`j' _hi`j' `xv', color(navy%25) lwidth(none))
            local obs
            if "`observed'" != "" local obs (line _o`j' `xv', lcolor(maroon) lpattern(dash) lwidth(medium))
            local xline
            local t = `TURN'[`tr', `j']
            if "`noturn'" == "" & `t' < . & `t' >= `xmin' & `t' <= `xmax' {
                local xline xline(`t', lcolor(gs8) lpattern(shortdash))
                local anyturn 1
            }
            tempname g`j'
            twoway `band' (line _w`j' `xv', lcolor(navy) lpattern(solid) lwidth(medthick)) `obs', ///
                title("`: word `j' of `names''", size(medsmall)) ///
                ytitle("Budget share", size(vsmall)) xtitle("") ylabel(, labsize(vsmall) angle(0)) ///
                `xl' `xline' legend(off) graphregion(color(white)) ///
                name(`g`j'', replace) nodraw
            local plots `plots' `g`j''
        }
        local nt `""Budget share, `what'""'
        local nt2
        if !`asobs' {
            if `quad' local nt2 "quadratic in ln x"
            else      local nt2 "linear in ln x (AIDS)"
            if `anyturn' local nt2 "`nt2'; vertical line: turning point"
        }
        if "`observed'" != "" local nt2 "dashed: observed shares, same smoother"
        if `trim' > 0 {
            if "`nt2'" != "" local nt2 "`nt2'; "
            local nt2 "`nt2'tails trimmed at `trim'%"
        }
        if "`nt2'" != "" local nt `"`nt' "`nt2'""'
        if `ci' local nt `"`nt' "`level'% confidence band, delta method, vce(`e(vce)')""'
        local ncol = ceil(sqrt(`M'))
        local nrow = ceil(`M' / `ncol')
        local grid
        if !strpos(`"`options'"', "cols(") & !strpos(`"`options'"', "rows(") local grid cols(`ncol')
        local gsize
        if !strpos(`"`options'"', "xsize(") & !strpos(`"`options'"', "ysize(") {
            local gsize xsize(`=min(2.3 * `ncol', 12)') ysize(`=min(1.9 * `nrow' + 1, 12)')
        }
        graph combine `plots', `grid' `gsize' ///
            title("Engel curves, `e(model)' (equaids)") ///
            b1title("`xt'", size(small)) ///
            note(`nt', size(vsmall)) graphregion(color(white)) `options'
        if `"`saving'"' != "" graph save `saving'
        graph drop `plots'
    }
    restore
    return scalar n = `n'
    if `asobs' return scalar bwidth = `bwidth'
    else       return matrix turn = `TURN'
end

program define _equaids_estat_diag, rclass
    syntax
    di _n as txt "Diagnostics of the estimate: " as res "`e(model)'" as txt ", " ///
        as res e(N) as txt " households, alpha_0 = " as res %8.4f e(anot) ///
        as txt " (" as res "`e(anot_rule)'" as txt ")"
    di as txt "{hline 78}"
    di as txt "Convergence" _col(40) as res cond(e(converged), "yes", "no") ///
        as txt ", " as res e(iter) as txt " iterations"
    di as txt "  Newton decrement g'A^-1 g" _col(40) as res %10.2e e(nrgrad)
    di as txt "  relative change of Sigma" _col(40) as res %10.2e e(sigdif)
    di as txt "  stopped for lack of progress" _col(40) as res cond(e(stalled), "yes", "no")
    di as txt "  steps cut at the boundary m0 > 0" _col(40) as res %10.0f e(n_boundary)
    di as txt "Information matrix, scaled rcond" _col(40) as res %10.2e e(rcond) ///
        as txt cond(e(rcond) < 1e-5, "  (nearly singular)", "")
    di as txt "Regressors, condition number" _col(40) as res %10.1f e(cond_x)
    if e(ndemos) > 0 di as txt "Ray scaling, smallest m0(z)" _col(40) as res %10.4f e(m0_min) ///
        as txt "  (below 0.05: " as res e(n_m0low) as txt " households)"
    di as txt "Households with ln x below ln m0 + ln a(p)" _col(40) as res %10.0f e(n_lneg)
    di as txt "Households with predicted shares outside [0,1]" _col(40) as res %10.0f e(n_shout)
    di as txt "{hline 78}"
    * the same warnings and notes as at estimation
    equaids, _warnings
    return scalar converged = e(converged)
    return scalar rcond     = e(rcond)
    return scalar n_lneg    = e(n_lneg)
    return scalar n_shout   = e(n_shout)
end

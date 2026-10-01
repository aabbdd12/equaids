*! equaids_examples 1.2.1  2026-10-01  Abdelkrim Araar
*! The examples of help equaids and help equaidsdiag, run from their links.
*!   equaids_examples #          run example # in the command window
*!   equaids_examples #, db      open the dialog box of equaids filled in for example #
*!   equaids_examples #, do      open example # as a do-file in the Do-file Editor
*! The data in memory are never lost: the run keeps them (preserve) and gives
*! them back at the end, even after an error or a Break; the do-file does the
*! same; the dialog box, which needs the example data in memory, refuses to
*! replace data of the user that have unsaved changes (the example data it
*! loads are marked and can be replaced). Files written by the examples go to
*! Stata's temporary folder, c(tmpdir), never to the working folder.
*! The example data of the package are ancillary files: read from the
*! current folder (where "ssc install equaids, all" or "net get equaids"
*! copies them), else from the SSC archive, else from GitHub
*! (_equaids_exload); nothing is written.
program define equaids_examples
    version 14.2
    syntax anything(name=ex id="example number") [, DB DO NOEDIT]
    capture confirm integer number `ex'
    if _rc | !inrange(`ex', 1, 8) {
        di as err "equaids_examples: the examples are numbered 1 to 8 (see help equaids)"
        exit 198
    }
    if "`db'" != "" & "`do'" != "" {
        di as err "equaids_examples: db and do cannot be combined"
        exit 198
    }
    local T = c(tmpdir)
    local T : subinstr local T "\" "/", all
    if substr("`T'", -1, 1) != "/" local T "`T'/"
    local F "w1-w4, prices(p1-p4) expenditure(expfd)"
    local W "wcorn wwheat wrice wother wcomp"
    local P "pcorn pwheat price pother pcomp"
    local n 0
    * the data of each example: Poi's food data (webuse), the Mexican cereals or
    * the simulated non-buyers, ancillary files of the package (_equaids_exload)
    local data "webuse food, clear"
    local dfile ""
    if `ex' == 1 {
        local title "The elasticities of four food groups (Poi's data)"
        local c1 "equaids `F' snames(meat fruitveg bread dairy)"
        local c2 "equaids, compensated checks stars"
        local n 2
    }
    else if `ex' == 2 {
        local title "The types of elasticities: households (the default), market, reference, hhmean"
        local c1 "equaids `F' snames(meat fruitveg bread dairy) notable"
        local c2 "equaids, elasticities(households)"
        local c3 "equaids, elasticities(market)"
        local c4 "equaids, elasticities(reference)"
        local c5 "equaids, elasticities(hhmean)"
        local n 5
    }
    else if `ex' == 3 {
        local dfile "mexico_2014_cereals"
        local data "use mexico_2014_cereals, clear"
        local title "Survey design (Mexican cereals)"
        local c1 "svyset psu [pweight=sweight], strata(strata) vce(linearized) singleunit(centered)"
        local c2 "equaids `W', prices(`P') expenditure(hh_current_inc) demographics(hhsize isMale) vce(svy)"
        local n 2
    }
    else if `ex' == 4 {
        local dfile "mexico_2014_cereals"
        local data "use mexico_2014_cereals, clear"
        local title "The elasticities of the individual (Mexican cereals)"
        local c1 "equaids `W' [aw=sweight], prices(`P') expenditure(hh_current_inc) hhsize(hhsize)"
        local n 1
    }
    else if `ex' == 5 {
        local title "Engel curves after estimation"
        local c1 "equaids `F' notable"
        local c2 "estat engel"
        local c3 `"estat engel, lnx level(90) data("`T'equaids_curves", replace)"'
        local c4 "estat engel, asobserved observed"
        local n 4
    }
    else if `ex' == 6 {
        local title "Diagnose a specification before estimating it"
        local c1 "equaidsdiag `F' sensitivity"
        local n 1
    }
    else if `ex' == 7 {
        local dfile "mexico_2014_cereals"
        local data "use mexico_2014_cereals, clear"
        local title "The non-buyers (Mexican cereals)"
        local c1 "equaids `W' [pw=sweight], prices(`P') expenditure(hh_current_inc) demographics(hhsize isMale) pimpute(psu rururb) selection selvars(perc_ocupa) vce(bootstrap, reps(50) seed(1))"
        local c2 "estat engel"
        local n 2
    }
    else if `ex' == 8 {
        local dfile "equaids_nonbuyers"
        local data "use equaids_nonbuyers, clear"
        local title "The non-buyers with analytic standard errors (simulated data)"
        local c1 "equaids w1 w2 w3, prices(p1 p2 p3) expenditure(x) noquadratic demographics(hs) anot(0) pimpute(grp) selection selvars(w1: q1 ; w2: q2) vce(cluster grp)"
        local c2 `"display "e(m0_t), distance to the boundary m0(z) > 0: " %4.1f e(m0_t) " standard errors""'
        local n 2
    }

    * ---- as a do-file, in Stata's temporary folder ----
    if "`do'" != "" {
        local fn "`T'equaids_example_`ex'.do"
        tempname fh
        file open `fh' using "`fn'", write text replace
        file write `fh' "* equaids, example `ex': `title'" _n
        file write `fh' "* Written by equaids_examples in Stata's temporary folder; save it elsewhere to keep it." _n
        file write `fh' "* preserve keeps the data in memory and gives them back when this do-file ends;" _n
        file write `fh' "* delete that line to keep working on the example data." _n
        file write `fh' "preserve" _n
        if "`dfile'" == "" file write `fh' "`data'" _n
        else {
            local l = substr("`dfile'", 1, 1)
            file write `fh' "* the example data: the current folder (where ssc install equaids, all copies" _n
            file write `fh' "* them), else the SSC archive, else GitHub" _n
            file write `fh' `"capture use `dfile', clear"' _n
            file write `fh' `"if _rc capture use "http://fmwww.bc.edu/repec/bocode/`l'/`dfile'.dta", clear"' _n
            file write `fh' `"if _rc use "https://raw.githubusercontent.com/aabbdd12/equaids/main/examples/`dfile'.dta", clear"' _n
        }
        forvalues i = 1/`n' {
            file write `fh' `"`c`i''"' _n
        }
        file close `fh'
        if "`noedit'" == "" doedit "`fn'"
        di as txt "(example `ex' written to " as res `"`fn'"' as txt ")"
        exit
    }

    * ---- in the dialog box: needs the example data in memory ----
    if "`db'" != "" {
        if inlist(`ex', 2, 5) {
            di as err "equaids_examples: example `ex' runs commands after the estimation; run it in the command window"
            exit 198
        }
        * the example data loaded for a dialog box carry a mark, here
        * _dta[equaids_example], and _dta[easi_example] or _dta[duvm_example]
        * from the examples of easi and duvm (easi builds prices in levels on
        * them, which sets c(changed)): they are replaced without asking; the
        * user's own data with changes are not
        local isex 0
        local cl : char _dta[]
        foreach c of local cl {
            if substr("`c'", -8, .) == "_example" & `"`: char _dta[`c']'"' == "1" local isex 1
        }
        if c(changed) & !`isex' {
            di as err "equaids_examples, db: the data in memory have changes not saved;"
            di as err "save them (or clear) first: the dialog box needs the example data in memory"
            exit 4
        }
        if "`dfile'" == "" qui `data'
        else _equaids_exload `dfile'
        char _dta[equaids_example] "1"
        if `ex' == 3 qui svyset psu [pweight=sweight], strata(strata) vce(linearized) singleunit(centered)
        di as txt "(example data loaded for the dialog box)"
        db equaids
        * Stata keeps the state of a dialog between two openings: every control
        * an example may set is first put back to its default
        .equaids_dlg.main.cb_act.setvalue "est"
        .equaids_dlg.main.rb_plev.seton
        .equaids_dlg.main.rb_xlev.seton
        .equaids_dlg.main.rb_quaids.seton
        .equaids_dlg.main.ed_snames.setvalue ""
        .equaids_dlg.main.vl_demo.setvalue ""
        .equaids_dlg.main.ck_anot.setoff
        .equaids_dlg.main.ed_anot.setvalue ""
        .equaids_dlg.weights.vl_wgt.setvalue ""
        .equaids_dlg.weights.rb_none.seton
        .equaids_dlg.se.rb_robust.seton
        .equaids_dlg.se.vn_clust.setvalue ""
        .equaids_dlg.se.sp_reps.setvalue 200
        .equaids_dlg.se.ed_seed.setvalue ""
        .equaids_dlg.se.ck_bsvy.setoff
        .equaids_dlg.rpt.cb_el.setvalue "households"
        .equaids_dlg.rpt.vn_hhs.setvalue ""
        .equaids_dlg.rpt.ck_comp.setoff
        .equaids_dlg.rpt.ck_checks.setoff
        .equaids_dlg.rpt.ck_stars.setoff
        .equaids_dlg.rpt.ck_notab.setoff
        .equaids_dlg.rpt.fi_save.setvalue ""
        .equaids_dlg.dg.ck_sens.setoff
        .equaids_dlg.main.vl_pimp.setvalue ""
        .equaids_dlg.sel.ed_selg.setvalue ""
        .equaids_dlg.sel.vl_all.setvalue ""
        forvalues r = 1/10 {
            .equaids_dlg.sel.ed_g`r'.setvalue ""
            .equaids_dlg.sel.vl_v`r'.setvalue ""
        }
        .equaids_dlg.sel.cb_ng.setvalue "0"
        .equaids_dlg.sel.ck_sel.setoff
        * the example
        if `ex' == 1 {
            .equaids_dlg.main.vl_shares.setvalue "w1 w2 w3 w4"
            .equaids_dlg.main.vl_prices.setvalue "p1 p2 p3 p4"
            .equaids_dlg.main.vn_exp.setvalue "expfd"
            .equaids_dlg.main.ed_snames.setvalue "meat fruitveg bread dairy"
        }
        if `ex' == 3 {
            .equaids_dlg.main.vl_shares.setvalue "`W'"
            .equaids_dlg.main.vl_prices.setvalue "`P'"
            .equaids_dlg.main.vn_exp.setvalue "hh_current_inc"
            .equaids_dlg.main.vl_demo.setvalue "hhsize isMale"
            .equaids_dlg.se.rb_svy.seton
        }
        if `ex' == 4 {
            .equaids_dlg.main.vl_shares.setvalue "`W'"
            .equaids_dlg.main.vl_prices.setvalue "`P'"
            .equaids_dlg.main.vn_exp.setvalue "hh_current_inc"
            .equaids_dlg.weights.rb_aw.seton
            .equaids_dlg.weights.vl_wgt.setvalue "sweight"
            .equaids_dlg.rpt.cb_el.setvalue "individuals"
            .equaids_dlg.rpt.vn_hhs.setvalue "hhsize"
        }
        if `ex' == 7 {
            .equaids_dlg.main.vl_shares.setvalue "`W'"
            .equaids_dlg.main.vl_prices.setvalue "`P'"
            .equaids_dlg.main.vn_exp.setvalue "hh_current_inc"
            .equaids_dlg.main.vl_demo.setvalue "hhsize isMale"
            .equaids_dlg.main.vl_pimp.setvalue "psu rururb"
            .equaids_dlg.weights.rb_pw.seton
            .equaids_dlg.weights.vl_wgt.setvalue "sweight"
            .equaids_dlg.sel.ck_sel.seton
            .equaids_dlg.sel.vl_all.setvalue "perc_ocupa"
            .equaids_dlg.se.rb_boot.seton
            .equaids_dlg.se.sp_reps.setvalue 50
            .equaids_dlg.se.ed_seed.setvalue "1"
        }
        if `ex' == 8 {
            .equaids_dlg.main.vl_shares.setvalue "w1 w2 w3"
            .equaids_dlg.main.vl_prices.setvalue "p1 p2 p3"
            .equaids_dlg.main.vn_exp.setvalue "x"
            .equaids_dlg.main.rb_aids.seton
            .equaids_dlg.main.vl_demo.setvalue "hs"
            .equaids_dlg.main.ck_anot.seton
            .equaids_dlg.main.ed_anot.setvalue "0"
            .equaids_dlg.main.vl_pimp.setvalue "grp"
            .equaids_dlg.sel.ck_sel.seton
            .equaids_dlg.sel.cb_ng.setvalue "2"
            .equaids_dlg.sel.ed_g1.setvalue "w1"
            .equaids_dlg.sel.vl_v1.setvalue "q1"
            .equaids_dlg.sel.ed_g2.setvalue "w2"
            .equaids_dlg.sel.vl_v2.setvalue "q2"
            .equaids_dlg.se.rb_cluster.seton
            .equaids_dlg.se.vn_clust.setvalue "grp"
        }
        if `ex' == 6 {
            .equaids_dlg.main.cb_act.setvalue "diag"
            .equaids_dlg.main.vl_shares.setvalue "w1 w2 w3 w4"
            .equaids_dlg.main.vl_prices.setvalue "p1 p2 p3 p4"
            .equaids_dlg.main.vn_exp.setvalue "expfd"
            .equaids_dlg.dg.ck_sens.seton
        }
        exit
    }

    * ---- in the command window: the data in memory are kept ----
    preserve
    if "`dfile'" == "" {
        qui `data'
        local what "`data'"
    }
    else {
        _equaids_exload `dfile'
        local what "`dfile'.dta, read from `r(source)'"
    }
    di as txt _n "{hline 78}" _n "equaids, example `ex': " as res "`title'" _n as txt "{hline 78}"
    di as txt `"(example data: `what'; the data in memory come back at the end)"'
    forvalues i = 1/`n' {
        di as txt _n `". `c`i''"'
        `c`i''
    }
    if `ex' == 5 di as txt _n `"(files written to `T')"'
end

* ============================================================================
* the data of an example, an ancillary file of the package: from the current
* folder (where "ssc install equaids, all" or "net get equaids" copies it),
* else from the SSC archive, else from GitHub; nothing is written to disk
program define _equaids_exload, rclass
    args f
    capture confirm file "`f'.dta"
    if !_rc {
        quietly use "`f'.dta", clear
        return local source "the current folder"
        exit
    }
    local l = substr("`f'", 1, 1)
    capture quietly use "http://fmwww.bc.edu/repec/bocode/`l'/`f'.dta", clear
    if !_rc {
        return local source "the SSC archive"
        exit
    }
    capture quietly use "https://raw.githubusercontent.com/aabbdd12/equaids/main/examples/`f'.dta", clear
    if !_rc {
        return local source "GitHub"
        exit
    }
    di as err "equaids_examples: `f'.dta is not in the current folder (`c(pwd)'),"
    di as err "  and neither the SSC archive nor GitHub could be reached."
    di as txt "  Copy the example data into the current folder with"
    di as txt `"  {stata "ssc install equaids, all replace"} (from SSC), or"'
    di as txt `"  {stata "net get equaids, from(https://raw.githubusercontent.com/aabbdd12/equaids/main)"} (from GitHub)."'
    exit 601
end

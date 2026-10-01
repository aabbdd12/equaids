*! _equaids_pimpute 1.1.0  2026-09-27  Abdelkrim Araar
*! Missing log prices of the households of the sample (option pimpute() of
*! equaids): each is replaced by the weighted mean of the log prices of the
*! households of the same group that have one, the first grouping variable
*! first, the next ones when nobody in the group has a price. Households still
*! without a price leave the sample, and the means are computed again on the
*! sample so reduced, until it no longer changes: the donors are the
*! households of the final sample, so that estat, run on e(sample), finds the
*! same values. Used by equaids and by estat after equaids.
*!
*!   _equaids_pimpute lnpvars, touse(var) wt(var) groups(varlist) [names(list)]
*!
*! lnpvars are modified in place (they must be copies); touse is updated.
*! r(n_left): households that left the sample; r(notes): one note per good.
program define _equaids_pimpute, rclass
    version 14.2
    syntax varlist(numeric), TOUSE(varname) WT(varname) GROUPS(varlist) [NAMES(string)]
    local M : word count `varlist'
    if "`names'" == "" local names `varlist'
    * the log prices as observed, kept for every pass
    local k 0
    foreach v of local varlist {
        local ++k
        tempvar o`k'
        quietly gen double `o`k'' = `v'
    }
    local nleft 0
    local pass 0
    while 1 {
        local ++pass
        local k 0
        foreach v of local varlist {
            local ++k
            quietly replace `v' = `o`k''
            local lev 0
            foreach g of local groups {
                local ++lev
                tempvar num den
                quietly egen double `num' = total(cond(`touse' & !missing(`o`k''), `wt' * `o`k'', .)), by(`g')
                quietly egen double `den' = total(cond(`touse' & !missing(`o`k''), `wt', .)), by(`g')
                quietly count if `touse' & missing(`v') & `den' > 0 & !missing(`den') & !missing(`g')
                local n`k'_`lev' = r(N)
                quietly replace `v' = `num' / `den' if `touse' & missing(`v') & `den' > 0 & !missing(`den') & !missing(`g')
                drop `num' `den'
            }
        }
        * households still without a price for some good leave the sample
        tempvar miss
        quietly egen byte `miss' = rowmiss(`varlist') if `touse'
        quietly count if `touse' & `miss' > 0 & !missing(`miss')
        local drop = r(N)
        quietly replace `touse' = 0 if `touse' & `miss' > 0 & !missing(`miss')
        drop `miss'
        local nleft = `nleft' + `drop'
        if `drop' == 0 | `pass' >= 20 continue, break
    }
    * one note per good: the prices filled at each level
    local notes ""
    local k 0
    foreach v of local varlist {
        local ++k
        local nm : word `k' of `names'
        local t ""
        local lev 0
        local tot 0
        foreach g of local groups {
            local ++lev
            local tot = `tot' + `n`k'_`lev''
            local t "`t'`=cond("`t'" == "", "", ", ")'`=string(`n`k'_`lev'', "%12.0fc")' from `g'"
        }
        if `tot' > 0 local notes `"`notes' "(pimpute(): price of `nm' filled for `t')""'
        return scalar n_filled`k' = `tot'
    }
    if `nleft' > 0 local notes `"`notes' "(pimpute(): `=string(`nleft', "%12.0fc")' households without a price in any group left the sample)""'
    return scalar n_left = `nleft'
    return scalar passes = `pass'
    return local notes `"`notes'"'
end

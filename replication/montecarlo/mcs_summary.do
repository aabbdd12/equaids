* mcs_summary.do -- the summary of the Monte Carlo of mcs_equaids.do: Table
* tab:mc of the note (coverage of the 95% intervals and mean analytic
* standard error over the Monte Carlo standard deviation), and the numbers of
* its text (bias, far draws, iterations).  Reads raw/mcs_<cell>_truth.dta and
* raw/mcs_<cell>_1.dta (the replications, shipped: this summary runs in
* seconds; mcs_equaids.do recomputes them).  Writes out/mcs_summary.csv.
*
*   cd <repository>/replication/montecarlo
*   do mcs_summary.do
*
* For each cell and each family of statistics (the market expenditure,
* uncompensated and compensated elasticities, and the free coefficients
* theta): coverage = share of the replications whose interval
* estimate +/- c se covers the truth (the census), mean over the family;
* SE/SD = mean analytic standard error / standard deviation of the
* estimates, median over the family.  c = 1.96 under simple random sampling,
* the t quantile with the design degrees of freedom (190) under the survey
* design.  Monte Carlo bands (1,000 draws): coverage 0.95 +/- 0.017, SE/SD
* 1 +/- 0.056.
version 14.2
clear all
set more off
capture mkdir out
mata:
real rowvector _q(real matrix X, real scalar p)
{
    // quantile of each column, linear interpolation between order statistics
    real rowvector q
    real colvector v
    real scalar j, n, h, lo
    q = J(1, cols(X), .)
    for (j = 1; j <= cols(X); j++) {
        v = sort(X[., j], 1) ; n = rows(v)
        h = (n - 1) * p + 1 ; lo = floor(h)
        q[j] = (lo >= n ? v[n] : v[lo] + (h - lo) * (v[lo + 1] - v[lo]))
    }
    return(q)
}
real scalar _med(real rowvector x) return(_q(x', .5))
void _cell(string scalar tag, real scalar c, real scalar M, real matrix OUT)
{
    real matrix    T, D, E, S, TH, ST, X, SE
    real rowvector tr, tt, sd, cov, rat, bias, rs, md
    real scalar    K, P, f, lo, hi, R, far, it
    string rowvector v
    stata(`"quietly use "raw/mcs_"' + tag + `"_truth.dta", clear"')
    v = st_varname(1..st_nvar())
    K = sum(regexm(v, "^e[0-9]+$")) ; P = sum(regexm(v, "^t[0-9]+$"))
    tr = st_data(1, "e" :+ strofreal(1..K)) ; tt = st_data(1, "t" :+ strofreal(1..P))
    stata(`"quietly use "raw/mcs_"' + tag + `"_1.dta", clear"')
    stata("quietly keep if conv == 1")
    R  = st_nobs()
    E  = st_data(., "e" :+ strofreal(1..K)) ; S  = st_data(., "s" :+ strofreal(1..K))
    TH = st_data(., "t" :+ strofreal(1..P)) ; ST = st_data(., "st" :+ strofreal(1..P))
    it = _med(st_data(., "iter")')
    // estimates more than 10 robust standard deviations from the median
    rs  = (_q(E, .75) - _q(E, .25)) :/ 1.349
    md  = _q(E, .5)
    far = sum(abs(E :- md) :> 10 :* rs)
    OUT = J(4, 6, .)
    for (f = 1; f <= 4; f++) {
        if (f == 4) {
            X = TH ; SE = ST ; T = tt
        }
        else {
            lo = (f == 1 ? 1 : (f == 2 ? M + 1 : M + M * M + 1))
            hi = (f == 1 ? M : (f == 2 ? M + M * M : K))
            X = E[., lo..hi] ; SE = S[., lo..hi] ; T = tr[lo..hi]
        }
        sd   = sqrt(diagonal(variance(X)))'
        cov  = colsum(abs(X :- T) :<= c :* SE) :/ R
        rat  = colsum(SE) :/ R :/ sd
        bias = abs(colsum(X) :/ R :- T) :/ sd
        OUT[f, .] = (mean(cov'), _med(rat), max(bias), R, it, far)
    }
}
end

tempname fh
file open `fh' using "out/mcs_summary.csv", write text replace
file write `fh' "model,errors,design,family,coverage,se_sd,max_bias_sd,R,median_iter,far_draws" _n
local fam "expenditure uncompensated compensated theta"
mata: OUT = J(0, 0, .)
display as text _n "Table tab:mc: coverage (range over the four families) and SE/SD (range)"
display as text "{hline 72}"
local cmin 1
local cmax 0
local rmin 9
local rmax 0
foreach m in aids quaids {
    foreach e in homo hetero clus {
        foreach d in srs svy {
            local c = cond("`d'" == "svy", invttail(190, .025), 1.96)
            mata: _cell("`m'_`e'_`d'_n5000_psu25", `c', 4, OUT) ; st_matrix("OUT", OUT)
            mata: st_numscalar("c0", min(OUT[., 1])); st_numscalar("c1", max(OUT[., 1])); ///
                st_numscalar("r0", min(OUT[., 2])); st_numscalar("r1", max(OUT[., 2]))
            display as text %-7s "`m'" %-8s "`e'" %-5s "`d'" "  coverage " as result %5.3f c0 "--" %5.3f c1 ///
                as text "   SE/SD " as result %5.3f r0 "--" %5.3f r1
            forvalues f = 1/4 {
                file write `fh' "`m',`e',`d',`: word `f' of `fam''"
                forvalues k = 1/6 {
                    file write `fh' "," %12.6f (OUT[`f', `k'])
                }
                file write `fh' _n
            }
        }
    }
}
file close `fh'
* the numbers of the text
preserve
quietly import delimited using "out/mcs_summary.csv", clear
foreach d in srs svy {
    quietly summarize coverage if design == "`d'"
    local c0 = r(min)
    local c1 = r(max)
    quietly summarize se_sd if design == "`d'"
    display as text cond("`d'" == "srs", "random sampling: ", "survey design:   ") "coverage " ///
        as result %5.3f `c0' "--" %5.3f `c1' as text ", SE/SD " as result %5.3f r(min) "--" %5.3f r(max)
}
quietly summarize max_bias_sd
display as text "largest |bias| / Monte Carlo SD: " as result %5.2f r(max)
quietly summarize median_iter
display as text "median iterations: " as result r(min) "--" r(max)
quietly summarize far_draws
display as text "estimates more than 10 robust SD from the median: " as result r(max)
quietly summarize r
display as text "converged replications per cell: " as result r(min) "--" r(max)
restore

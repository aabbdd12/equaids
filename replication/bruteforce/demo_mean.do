*! demo_mean.do -- the influence function by brute force, on the simplest
*! estimator: the mean.  Pedagogical: 10 households, every step printed.
*!
*! theta(w) = sum_h w_h x_h / sum_h w_h.  Move the weight of household h to
*! 1 + eps and 1 - eps, all the others at 1, and re-run the estimator:
*!     U_h = [theta(1 + eps) - theta(1 - eps)] / (2 eps)  ~  d theta / d w_h
*! For the mean, d theta / d w_h = (x_h - xbar) / n exactly: U_h is the
*! influence function divided by n, and
*!     SE_bf = sqrt( n/(n-1) sum_h U_h^2 ) = s / sqrt(n),
*! the textbook standard error.  The jackknife removes household h
*! altogether (weight 1 -> 0, a finite step instead of a derivative):
*!     theta_(-h) - theta = -(x_h - xbar) / (n - 1),
*! and for the mean it gives s / sqrt(n) too.  The brute force is the
*! "infinitesimal jackknife" (Jaeckel 1972; Efron 1982).
clear all
set more off
set seed 20260928
quietly set obs 10
quietly gen double x = round(exp(rnormal(3, 0.5)), 0.1)
quietly gen double w = 1
local n = _N
quietly summarize x
local xbar = r(mean)
local s = r(sd)

di as txt _n "household" _col(12) "x_h" _col(22) "U_h brute force" _col(40) "(x_h - xbar)/n" _col(58) "jackknife: theta(-h) - theta"
local eps 1e-3
mata: U = J(`n', 1, .); JK = J(`n', 1, .)
forvalues h = 1/`n' {
	* the brute force: two re-runs with the weight of h moved by +/- eps
	quietly replace w = 1 + `eps' in `h'
	quietly summarize x [iw = w], meanonly
	local tp = r(mean)
	quietly replace w = 1 - `eps' in `h'
	quietly summarize x [iw = w], meanonly
	local tm = r(mean)
	quietly replace w = 1 in `h'
	local U = (`tp' - `tm') / (2 * `eps')
	* the jackknife: household h left out
	quietly summarize x if _n != `h', meanonly
	local jk = r(mean) - `xbar'
	mata: U[`h'] = `U'; JK[`h'] = `jk'
	di as txt %5.0f `h' _col(10) as res %7.1f x[`h'] _col(22) %12.6f `U' _col(40) %12.6f (x[`h'] - `xbar') / `n' ///
		_col(58) %12.6f `jk'
}
mata: st_numscalar("__sebf", sqrt(`n' / (`n' - 1) * sum(U:^2)))
mata: st_numscalar("__sejk", sqrt((`n' - 1) / `n' * sum((JK :- mean(JK)):^2)))
quietly jackknife m = r(mean), notable nodots: summarize x
tempname V
matrix `V' = e(V)
di as txt _n "standard error of the mean"
di as txt "  textbook       s / sqrt(n)                    " as res %12.8f `s' / sqrt(`n')
di as txt "  brute force    sqrt(n/(n-1) sum U_h^2)        " as res %12.8f __sebf
di as txt "  jackknife      sqrt((n-1)/n sum (dtheta)^2)   " as res %12.8f __sejk
di as txt "  Stata's jackknife prefix                      " as res %12.8f sqrt(`V'[1, 1])

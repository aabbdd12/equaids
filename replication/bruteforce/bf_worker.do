*! bf_worker.do -- the brute-force influence function, one block of
*! households, in its own Stata process (several run in parallel).
*!   cd <repository>/replication/bruteforce
*!   do bf_worker.do spec first last tag k [outdir]
*! spec.do defines three programs:
*!   _bf_setup     loads the data, creates the weight variable w = 1 (and
*!                 whatever the estimator needs);
*!   _bf_theta     runs the estimator with [iw = w] and leaves the vector of
*!                 estimates in the Mata row vector __bf_theta;
*!   _bf_analytic  (used by bf_combine.do only) leaves the analytic standard
*!                 errors in __bf_se and their names in __bf_names.
*! For each household h in first..last, the whole estimator is re-run with
*! its weight at 1 + eps and at 1 - eps (eps = $BF_EPS if the spec sets it,
*! else 1e-3); U_h = d theta / d w_h by central
*! differences.  The block U is saved to raw/bf_<tag>_<k>.mmat (outdir: another
*! folder), with CV the convergence flags of the two re-runs of each household.
args spec first last tag k outdir
clear all
set more off
local BF = subinstr("`c(pwd)'", "\", "/", .)
global BF_ROOT "`BF'/../.."
capture confirm file "$BF_ROOT/src/equaids.ado"
if _rc {
	display as error "run this script from the replication/bruteforce/ directory"
	exit 601
}
if "`outdir'" == "" local outdir "raw"
capture mkdir "`BF'/`outdir'"
run "`BF'/`spec'.do"
_bf_setup
timer clear 1
timer on 1
_bf_theta
mata: T0 = __bf_theta; U = J(`last' - `first' + 1, cols(T0), .); CV = J(`last' - `first' + 1, 2, .)
local eps = cond("$BF_EPS" == "", 1e-3, 0$BF_EPS)
forvalues h = `first'/`last' {
	quietly replace w = 1 + `eps' in `h'
	_bf_theta
	local cp = e(converged)
	mata: TP = __bf_theta
	quietly replace w = 1 - `eps' in `h'
	_bf_theta
	mata: U[`h' - `first' + 1, .] = (TP - __bf_theta) :/ (2 * `eps'); CV[`h' - `first' + 1, .] = (`cp', `=e(converged)')
	quietly replace w = 1 in `h'
}
timer off 1
quietly timer list 1
di as txt "block `k': households `first' to `last', " %6.1f r(t1) " s"
mata: B0 = `first'
mata: mata matsave `BF'/`outdir'/bf_`tag'_`k' U T0 B0 CV, replace

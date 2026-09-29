*! gn_equaids110.do -- the Gauss-Newton column of Table tab:numif: the robust
*! standard errors reported by equaids 1.1.0, the previous release (A^-1 in
*! the influence function of the coefficients, Sigma held fixed, the expected
*! information of the probits, the filled prices taken as data), on the three
*! specifications of the brute force.  The estimator did not change between
*! the releases, only its variance: the estimates are checked against those
*! of the current version.
*!
*!   cd <repository>/replication/bruteforce
*!   do gn_equaids110.do
*!
*! Downloads the ado-files of the release v1.1.0 from GitHub into v110/ once
*! (or put them there by hand), then writes raw/gn110_<tag>.csv (name,
*! standard error), read by numif_summary.do.  equaids 1.1.0 has no
*! "households" family (its households were the reference household): the
*! Gauss-Newton column covers the other estimates.  Each specification runs
*! in a call of this file with its tag (do gn_equaids110.do poi).
args tag
version 14.2
set more off
local BF = subinstr("`c(pwd)'", "\", "/", .)
capture confirm file "`BF'/../../src/equaids.ado"
if _rc {
	display as error "run this script from the replication/bruteforce/ directory"
	exit 601
}
global BF_ROOT "`BF'/../.."
if "`tag'" == "" {
	local URL "https://raw.githubusercontent.com/aabbdd12/equaids/v1.1.0/src"
	capture mkdir "`BF'/v110"
	foreach f in equaids.ado equaids_estat.ado _equaids_tabstars.ado _equaids_pimpute.ado {
		capture confirm file "`BF'/v110/`f'"
		if _rc copy "`URL'/`f'" "`BF'/v110/`f'"
	}
	capture mkdir "`BF'/raw"
	foreach t in poi sel selimp {
		do "`BF'/gn_equaids110.do" `t'
	}
	clear all
	capture adopath - "`BF'/v110"
	global BF_SRC
	global BF_IMP
	global GN_B12
	global BF_ROOT
	exit
}

* ---- one specification ----
clear all
global BF_IMP
* the estimates of the current version, for the check
global BF_SRC "`BF'/../../src"
run "`BF'/spec_`tag'.do"
_bf_setup
quietly equaids $BF_SPEC noelastse
* (clear all drops matrices, not global macros)
mata: st_global("GN_B12", invtokens(strofreal(st_matrix("e(b)"), "%21.0g")))
clear all
capture adopath - "`BF'/../../src"
global BF_IMP
global BF_SRC "`BF'/v110"
run "`BF'/spec_`tag'.do"
_bf_setup
quietly equaids $BF_SPEC
mata:
// the full names of a stripe, as rowfullnames/colfullnames give them
string rowvector _gn_full(string matrix S)
{
	string rowvector r
	real scalar i
	r = J(1, rows(S), "")
	for (i = 1; i <= rows(S); i++) r[i] = (S[i, 1] == "" ? S[i, 2] : S[i, 1] + ":" + S[i, 2])
	return(r)
}
void _gn_write(string scalar csv)
{
	string rowvector n, rn
	real rowvector   s
	string scalar    f, t
	real scalar      k, j, fh
	n = _gn_full(st_matrixcolstripe("e(b)"))
	s = sqrt(diagonal(st_matrix("e(V)")))'
	if (st_global("e(selection)") != "") {
		n = n, "psi" :+ strofreal(1..cols(st_matrix("e(sel_psi)")))
		s = s, sqrt(diagonal(st_matrix("e(V_sel_psi)")))'
	}
	for (k = 1; k <= 3; k++) {
		f = ("", "m", "h")[k]
		for (j = 1; j <= 3; j++) {
			t = ("x", "u", "c")[j]
			rn = _gn_full(st_matrixrowstripe("e(V_elas_" + t + f + ")"))
			n = n, (t + f + ":") :+ rn
			s = s, sqrt(diagonal(st_matrix("e(V_elas_" + t + f + ")")))'
		}
	}
	if (fileexists(csv)) unlink(csv)
	fh = fopen(csv, "w")
	fput(fh, "name,se")
	for (k = 1; k <= cols(n); k++) fput(fh, sprintf(`""%s",%21.15g"', n[k], s[k]))
	fclose(fh)
	printf("{txt}  %g standard errors -> %s\n", cols(n), csv)
}
end
mata: st_numscalar("__d", mreldif(st_matrix("e(b)"), strtoreal(tokens(st_global("GN_B12")))))
display as text "`tag': equaids 1.1.0, N = " as result e(N) ///
	as text ", estimates of 1.1.0 and of the current version: max rel. diff. " as result %9.2e __d
mata: _gn_write("`BF'/raw/gn110_`tag'.csv")

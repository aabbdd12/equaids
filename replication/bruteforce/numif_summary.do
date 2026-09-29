*! numif_summary.do -- Table tab:numif of the note: the robust standard errors
*! divided by the brute-force ones, for equaids (current version) and for the
*! Gauss-Newton variance of equaids 1.1.0.
*!
*!   cd <repository>/replication/bruteforce
*!   do numif_summary.do
*!
*! The brute force is read from raw/bf_<tag>_<k>.mmat (bf_parallel.sh or
*! bf_worker.do recompute it: hours), the Gauss-Newton standard errors from
*! raw/gn110_<tag>.csv (gn_equaids110.do); the standard errors of the current
*! version are recomputed here (one estimation per specification, seconds), so
*! that a change of the variance formula is scored against the same brute
*! force.  SE_bf = sqrt(n/(n-1) sum_h U_h^2).  Writes out/numif_summary.csv.
*! Each specification runs in a call of this file with its tag.
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
	capture mkdir "`BF'/out"
	tempname fh
	file open `fh' using "`BF'/out/numif_summary.csv", write text replace
	file write `fh' "spec,households,estimates,equaids_min,equaids_max,equaids_maxdev,gn_estimates,gn_min,gn_max,gn_coef_min,gn_coef_max" _n
	file close `fh'
	display as text _n "Table tab:numif: reported robust SE / brute-force SE"
	foreach t in poi sel selimp {
		do "`BF'/numif_summary.do" `t'
	}
	clear all
	capture adopath - "`BF'/../../src"
	global BF_SRC
	global BF_IMP
	global BF_ROOT
	exit
}

* ---- one specification ----
clear all
global BF_IMP
global BF_SRC "`BF'/../../src"
run "`BF'/spec_`tag'.do"
_bf_setup
_bf_analytic
capture confirm file "`BF'/raw/gn110_`tag'.csv"
local gn = (_rc == 0)
if `gn' {
	quietly import delimited using "`BF'/raw/gn110_`tag'.csv", clear varnames(1) stringcols(1)
	mata: __gn_n = st_sdata(., "name")'; __gn_s = st_data(., "se")'
}
else mata: __gn_n = J(1, 0, ""); __gn_s = J(1, 0, .)
mata: UA = J(0, 0, .)
forvalues k = 1/9 {
	capture confirm file "`BF'/raw/bf_`tag'_`k'.mmat"
	if !_rc {
		mata: mata matuse "`BF'/raw/bf_`tag'_`k'", replace
		mata: UA = (rows(UA) ? UA \ U : U)
	}
}
mata:
void _numif(string scalar tag, string scalar BF, real matrix UA, real rowvector se,
	string rowvector names, string rowvector block, string rowvector gnn, real rowvector gns)
{
	real rowvector sebf, rat, ok, rg, cg, okg
	real scalar    j, n, fh, i
	n    = rows(UA)
	sebf = sqrt(n / (n - 1) :* colsum(UA:^2))
	ok   = (sebf :> 0) :& (se :> 0) :& (se :< .)
	rat  = se :/ sebf
	// Gauss-Newton: matched by name
	rg = J(1, cols(names), .); cg = J(1, cols(names), 0)
	for (j = 1; j <= cols(gnn); j++) {
		for (i = 1; i <= cols(names); i++) {
			if (names[i] == gnn[j]) {
				rg[i] = gns[j] / sebf[i]
				cg[i] = (substr(block[i], 1, 12) == "coefficients")
				break
			}
		}
	}
	okg = ok :& (rg :< .)
	printf("{txt}%-8s households {res}%5.0f{txt}, estimates {res}%4.0f{txt}: equaids {res}%6.4f{txt}--{res}%6.4f{txt} (largest deviation {res}%7.1e{txt})",
		tag, n, sum(ok), min(select(rat, ok)), max(select(rat, ok)), max(select(abs(rat :- 1), ok)))
	if (sum(okg)) {
		printf("{txt}; Gauss-Newton ({res}%g{txt}) {res}%5.3f{txt}--{res}%5.3f{txt}, coefficients {res}%5.3f{txt}--{res}%5.3f\n",
			sum(okg), min(select(rg, okg)), max(select(rg, okg)),
			min(select(rg, okg :& cg)), max(select(rg, okg :& cg)))
	}
	else printf("\n")
	fh = fopen(BF + "/out/numif_summary.csv", "a")
	fput(fh, sprintf("%s,%g,%g,%12.6f,%12.6f,%12.3e,%g,%12.6f,%12.6f,%12.6f,%12.6f", tag, n, sum(ok),
		min(select(rat, ok)), max(select(rat, ok)), max(select(abs(rat :- 1), ok)), sum(okg),
		(sum(okg) ? min(select(rg, okg)) : .), (sum(okg) ? max(select(rg, okg)) : .),
		(sum(okg) ? min(select(rg, okg :& cg)) : .), (sum(okg) ? max(select(rg, okg :& cg)) : .)))
	fclose(fh)
}
end
mata: _numif("`tag'", "`BF'", UA, __bf_se, __bf_names, __bf_block, __gn_n, __gn_s)

*! bf_combine.do -- stacks the blocks of bf_worker.do and compares the
*! analytic standard errors with the brute-force ones.
*!   cd <repository>/replication/bruteforce
*!   do bf_combine.do spec tag K [indir]
*! The blocks are read from raw/ (indir: another folder).
*! SE_bf = sqrt( n/(n-1) sum_h U_h^2 ), the robust standard error in which
*! every step of the procedure moves with the weight of each household.
*! The ratios (analytic / brute force) are summarised by block (the Mata
*! string vector __bf_block of the spec, one label per estimate; one block if
*! absent), then the largest deviations are listed.  The whole table goes to
*! out/bf_<tag>.csv.
args spec tag K indir
clear all
set more off
local BF = subinstr("`c(pwd)'", "\", "/", .)
global BF_ROOT "`BF'/../.."
capture confirm file "$BF_ROOT/src/equaids.ado"
if _rc {
	display as error "run this script from the replication/bruteforce/ directory"
	exit 601
}
if "`indir'" == "" local indir "raw"
capture mkdir "`BF'/out"
run "`BF'/`spec'.do"
_bf_setup
_bf_analytic
capture mata: __bf_tmp = __bf_block
if _rc mata: __bf_block = J(1, cols(__bf_se), "all")
mata: UA = J(0, 0, .); CA = J(0, 2, .)
forvalues k = 1/`K' {
	mata: mata matuse `BF'/`indir'/bf_`tag'_`k', replace
	mata: UA = (rows(UA) ? UA \ U : U)
	capture mata: CA = CA \ CV
}
mata: st_numscalar("__bf_nnc", sum(!CA)); st_numscalar("__bf_nrr", 2 * rows(CA))
if __bf_nrr > 0 di as txt "re-runs not converged (stalled at the numerical floor): " __bf_nnc " of " __bf_nrr
mata:
real scalar __bf_med(real colvector x)
{
	real colvector v
	real scalar m
	v = sort(x, 1); m = rows(v)
	return(mod(m, 2) ? v[(m + 1) / 2] : (v[m / 2] + v[m / 2 + 1]) / 2)
}
void __bf_report(real matrix UA, real rowvector se, string rowvector names,
	string rowvector block, real scalar K, string scalar csv)
{
	real scalar n, j, b, fh, nw
	real rowvector sebf, rat, ok, sel, dev
	real colvector o
	string rowvector bl
	n    = rows(UA)
	sebf = sqrt(n / (n - 1) :* colsum(UA:^2))
	ok   = (sebf :> 0) :& (se :> 0) :& (se :< .)
	rat  = se :/ sebf
	rat  = rat :/ ok
	printf("\n{txt}brute force: %g households, %g blocks; analytic SE / brute-force SE\n", n, K)
	printf("{txt}  %-40s %6s %9s %9s %9s\n", "", "n", "median", "min", "max")
	bl = J(1, 0, "")
	for (j = 1; j <= cols(block); j++) if (!anyof(bl, block[j])) bl = bl, block[j]
	for (b = 1; b <= cols(bl); b++) {
		sel = (block :== bl[b]) :& ok
		if (!sum(sel)) continue
		printf("{txt}  %-40s {res}%6.0f %9.4f %9.4f %9.4f\n", bl[b], sum(sel),
			__bf_med(select(rat, sel)'), min(select(rat, sel)), max(select(rat, sel)))
	}
	printf("{txt}  %-40s {res}%6.0f %9.4f %9.4f %9.4f\n", "all", sum(ok),
		__bf_med(select(rat, ok)'), min(select(rat, ok)), max(select(rat, ok)))
	if (sum(!ok)) printf("{txt}  (%g estimates without variance left out)\n", sum(!ok))
	dev = abs(rat :- 1) :* ok
	dev = editmissing(dev, 0)
	o = order(-dev', 1)
	nw = min((8, sum(ok)))
	printf("{txt}  largest deviations:\n")
	for (j = 1; j <= nw; j++) {
		printf("{txt}    %-36s {res}%14.8g %14.8g %9.4f\n", block[o[j]] + " " + names[o[j]],
			se[o[j]], sebf[o[j]], rat[o[j]])
	}
	fh = fopen(csv, "w")
	fput(fh, "block,name,analytic,bruteforce,ratio")
	for (j = 1; j <= cols(rat); j++) {
		fput(fh, sprintf(`""%s","%s",%21.15g,%21.15g,%21.15g"', block[j], names[j], se[j], sebf[j], rat[j]))
	}
	fclose(fh)
	st_numscalar("__bf_maxdev", max(select(abs(rat :- 1), ok)))
}
end
capture erase "`BF'/out/bf_`tag'.csv"
mata: __bf_report(UA, __bf_se, __bf_names, __bf_block, `K', "`BF'/out/bf_`tag'.csv")

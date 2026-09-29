*! spec_sel.do -- brute force of equaids with the correction for the
*! non-buyers (Shonkwiler and Yen): the perfect DGP of the note's Monte Carlo
*! (replication/dgp_nonbuyers.do, rho = 0.3, beta = (-0.10, 0.08), every
*! parameter strongly identified), 1,000 households, AIDS with Ray's scaling
*! (hs), probits on (1, ln p - ln x, hs, q_i), iweights w = 1.
*! Two variants, one change at a time:
*!   spec_sel     the prices of the non-buyers observed (the DGP's prices):
*!                the probits and the system only;
*!   spec_selimp  ($BF_IMP = 1) the prices of the non-buyers missing and
*!                filled by pimpute(grp), as in the Monte Carlo.
*! theta = e(b), e(sel_psi) (the system and the probits stacked), then the
*! four families of elasticities (market, households, reference, household
*! mean).
capture program drop _bf_setup
program define _bf_setup
	* $BF_SRC: the folder of equaids (default the current version, src/)
	if "$BF_SRC" == "" global BF_SRC "$BF_ROOT/src"
	adopath ++ "$BF_SRC"
	run "$BF_SRC/equaids.ado"
	run "$BF_ROOT/replication/dgp_nonbuyers.do"
	_dgp_nonbuyers 1000 0.3 -0.10 0.08 20260928
	if "$BF_IMP" == "1" local imp "pimpute(grp)"
	else {
		* the DGP's prices of the non-buyers, before they were set missing
		quietly replace p1 = exp(lp1)
		quietly replace p2 = exp(lp2)
	}
	quietly gen double w = 1
	global BF_SPEC "w1 w2 w3 [iw = w], prices(p1 p2 p3) expenditure(x) noquadratic demographics(hs) anot(0) `imp' selection selvars(w1: q1 ; w2: q2) vce(robust) notable nolog tolerance(1e-12) iterate(1000)"
	global BF_EPS 1e-2
end
mata:
real rowvector __bf_rm(string scalar m) return(vec(st_matrix(m)')')
real rowvector __bf_vec()
{
	real rowvector t
	string scalar f
	real scalar k
	t = st_matrix("e(b)"), st_matrix("e(sel_psi)")
	for (k = 1; k <= 4; k++) {
		f = ("", "w", "m", "h")[k]
		t = t, __bf_rm("e(elas_x" + f + ")"), __bf_rm("e(elas_u" + f + ")"),
			__bf_rm("e(elas_c" + f + ")")
	}
	return(t)
}
end
capture program drop _bf_theta
program define _bf_theta
	quietly equaids $BF_SPEC noelastse novariance
	if !e(converged) di as err "not converged"
	mata: __bf_theta = __bf_vec()
end
capture program drop _bf_analytic
program define _bf_analytic
	quietly equaids $BF_SPEC
	mata: __bf_se = sqrt(diagonal(st_matrix("e(V)")))', sqrt(diagonal(st_matrix("e(V_sel_psi)")))'
	local cn : colfullnames e(b)
	mata: __bf_names = tokens(st_local("cn")); __bf_block = J(1, cols(__bf_names), "coefficients e(b)")
	* (under version 14, colsof() takes a matrix, not e())
	tempname tp tb
	matrix `tp' = e(sel_psi)
	matrix `tb' = e(b)
	local np = colsof(`tp')
	local nb = colsof(`tb')
	local ns = `np' - `nb'
	mata: __bf_names = __bf_names, "psi" :+ strofreal(1..`np'); ///
		__bf_block = __bf_block, J(1, `np', "psi: system + probits")
	foreach f in "" w m h {
		local fam = cond("`f'" == "", "market", cond("`f'" == "w", "households", cond("`f'" == "m", "reference", "household mean")))
		foreach t in x u c {
			local tn = cond("`t'" == "x", "expenditure", cond("`t'" == "u", "uncompensated", "compensated"))
			mata: __bf_se = __bf_se, sqrt(diagonal(st_matrix("e(V_elas_`t'`f')")))'
			local rn : rowfullnames e(V_elas_`t'`f')
			mata: __bf_names = __bf_names, ("`t'`f':") :+ tokens(st_local("rn")); ///
				__bf_block = __bf_block, J(1, cols(tokens(st_local("rn"))), "`fam', `tn'")
		}
	}
end

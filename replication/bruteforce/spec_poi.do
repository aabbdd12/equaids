*! spec_poi.do -- brute force of equaids without selection: Poi's (2012) food
*! data, 4,048 households, QUAIDS without demographics (the data have none), iweights w = 1.
*! theta = e(b), then for each family (market, households, reference, household mean)
*! the expenditure, uncompensated and compensated elasticities (row-major).
*! tolerance(1e-12): the Newton decrement below 1e-12 puts the estimate
*! within 1e-6 standard errors of the optimum (Sigma stalls near 1e-13, so
*! 1e-14 is not always reachable); with a step eps = 1e-2 of the weight the
*! estimates move by about 3e-4 standard errors, so each U_h is exact to
*! about 0.3%, and the error averages out in sum_h U_h^2.
capture program drop _bf_setup
program define _bf_setup
	* $BF_SRC: the folder of equaids (default the current version, src/)
	if "$BF_SRC" == "" global BF_SRC "$BF_ROOT/src"
	adopath ++ "$BF_SRC"
	run "$BF_SRC/equaids.ado"
	webuse food, clear
	quietly gen double w = 1
	global BF_SPEC "w1-w4 [iw = w], prices(p1-p4) expenditure(expfd) vce(robust) notable nolog tolerance(1e-12) iterate(1000)"
	global BF_EPS 1e-2
end
mata:
real rowvector __bf_rm(string scalar m) return(vec(st_matrix(m)')')
real rowvector __bf_vec()
{
	real rowvector t
	string scalar f
	real scalar k
	t = st_matrix("e(b)")
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
	local M = e(ngoods)
	local sh `e(lhs)'
	mata: __bf_se = sqrt(diagonal(st_matrix("e(V)")))'
	local cn : colfullnames e(b)
	mata: __bf_names = tokens(st_local("cn")); __bf_block = J(1, cols(__bf_names), "coefficients")
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

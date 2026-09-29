*! decomposition.do -- the paragraph after Table tab:numif: where the gaps of
*! the Gauss-Newton standard errors (equaids 1.1.0) come from.  The pieces of
*! the influence function are added one at a time and each is scored against
*! the brute force.
*!
*!   cd <repository>/replication/bruteforce
*!   do decomposition.do
*!
*! The iterated FGNLS estimate solves the stacked estimating equations
*!     sum_h psi_h(phi) = 0,  phi = (theta, vech Sigma [, alpha]),
*!     psi_h = ( G_h' Sigma^-1 u_h ; vech(u_h u_h') - vech Sigma
*!               [; lambda_h s_h, the scores of the probits] ),
*! and moving the weight of household h gives exactly
*!     d phi / d w_h = -H^-1 psi_h,   H = d sum_h psi_h / d phi',
*! H by central differences of the scores here (equaids has it analytic).
*!
*! Poi's data (spec_poi), the coefficients e(b):
*!   GN     A^-1 psi_theta, A = sum G' Sigma^-1 G (equaids 1.1.0);
*!   Hobs   -H_tt^-1 psi_theta: the observed Hessian, Sigma held fixed;
*!   full   -H^-1 psi: the observed Hessian and the estimation of Sigma.
*! The non-buyers (spec_sel, spec_selimp), psi = (theta, alpha) = e(sel_psi):
*!   GN     equaids 1.1.0 (raw/gn110_<tag>.csv, gn_equaids110.do);
*!   probit the observed Hessian of the probits only, the system as in 1.1.0
*!          (Gauss-Newton, Sigma fixed, the cross term through the residuals);
*!   full   the exact influence function, the filled prices taken as data:
*!          with pimpute(), its gap to the brute force is the weight of the
*!          imputation of the prices.
*! Reads the brute force from raw/bf_<tag>_<k>.mmat; writes
*! out/decomposition.csv (ranges of SE / brute-force SE by block).  Each
*! specification runs in a call of this file with its tag.
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
	file open `fh' using "`BF'/out/decomposition.csv", write text replace
	file write `fh' "spec,block,variance,estimates,min,max" _n
	file close `fh'
	foreach t in poi sel selimp {
		do "`BF'/decomposition.do" `t'
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
quietly equaids $BF_SPEC noelastse
local a0 = e(anot)
local sel = ("`e(selection)'" != "")
local imp = ("$BF_IMP" == "1")
matrix TH  = e(b_free)
matrix SIG = e(Sigma)
matrix EB  = e(b)
if `sel' matrix PSI = e(sel_psi)
local nb = colsof(EB)
quietly gen byte __es = e(sample)
local goods = cond(`sel', "w1 w2 w3", "w1 w2 w3 w4")
local prices = cond(`sel', "p1 p2 p3", "p1 p2 p3 p4")
local xvar = cond(`sel', "x", "expfd")
local j = 0
foreach p of local prices {
	local ++j
	quietly gen double __lp`j' = ln(`p')
}
* the log prices as filled at estimation (pimpute)
if `imp' {
	quietly _equaids_pimpute __lp1 __lp2 __lp3, touse(__es) wt(w) groups(grp)
}
quietly gen double __lx = ln(`xvar')
mata: UA = J(0, 0, .)
forvalues k = 1/9 {
	capture confirm file "`BF'/raw/bf_`tag'_`k'.mmat"
	if !_rc {
		mata: mata matuse "`BF'/raw/bf_`tag'_`k'", replace
		mata: UA = (rows(UA) ? UA \ U : U)
	}
}
mata:
es = st_data(., "__es")
W  = st_data(., tokens(st_local("goods")), "__es")
LP = st_data(., "__lp" :+ strofreal(1..cols(W)), "__es")
LX = st_data(., "__lx", "__es")
end
if `sel' {
	mata: Z = st_data(., "hs", "__es"); Q = st_data(., ("q1", "q2"), "__es")
}
else mata: Z = J(rows(LX), 0, .)
preserve
quietly import delimited using "`BF'/raw/gn110_`tag'.csv", clear varnames(1) stringcols(1)
mata: GNn = st_sdata(., "name")'; GNs = st_data(., "se")'
restore

mata:
real matrix _unvech(real rowvector v, real scalar L)
{
	real matrix S
	real scalar i, j, k
	S = J(L, L, .); k = 0
	for (j = 1; j <= L; j++) {
		for (i = j; i <= L; i++) {
			k++
			S[i, j] = v[k]; S[j, i] = v[k]
		}
	}
	return(S)
}
// PH, PHI of the two probits from alpha
void _phphi(real rowvector al, real matrix SD1, real matrix SD2, real matrix PH, real matrix PHI)
{
	real colvector x1, x2
	real scalar p1
	p1 = cols(SD1)
	x1 = SD1 * al[|1 \ p1|]'
	x2 = SD2 * al[|p1 + 1 \ cols(al)|]'
	PH  = normal(x1), normal(x2)
	PHI = normalden(x1), normalden(x2)
}
// the stacked psi_h: (theta, vech Sigma [, alpha]); al empty without selection
real matrix _psis(real rowvector t, real matrix Sig, real rowvector al,
	real matrix W, real matrix LP, real colvector LX, real matrix Z,
	real scalar a0, real scalar qd, real matrix D, real matrix SD1,
	real matrix SD2, real rowvector cix)
{
	real matrix U, US, Pp, Ps, PH, PHI, Pa, SD
	pointer(real matrix) rowvector G
	real colvector dd, xb, q, lam
	real scalar i, j, k, L, N, p1
	if (cols(al)) {
		_phphi(al, SD1, SD2, PH, PHI)
		(void) _eq_model(t, W, LP, LX, Z, a0, qd, D, U, G, 1, PH, PHI, cix)
	}
	else (void) _eq_model(t, W, LP, LX, Z, a0, qd, D, U, G, 1)
	L = cols(U); N = rows(U)
	US = U * invsym(Sig)
	Pp = J(N, cols(t), 0)
	for (i = 1; i <= L; i++) Pp = Pp + (*G[i]) :* US[., i]
	Ps = J(N, L * (L + 1) / 2, .)
	k = 0
	for (j = 1; j <= L; j++) {
		for (i = j; i <= L; i++) {
			k++
			Ps[., k] = U[., i] :* U[., j] :- Sig[i, j]
		}
	}
	if (!cols(al)) return((Pp, Ps))
	p1 = cols(SD1)
	Pa = J(N, 0, .)
	for (i = 1; i <= 2; i++) {
		SD = (i == 1 ? SD1 : SD2)
		xb = SD * (i == 1 ? al[|1 \ p1|] : al[|p1 + 1 \ cols(al)|])'
		dd = W[., i] :> 0
		q  = 2 :* dd :- 1
		lam = q :* normalden(q :* xb) :/ normal(q :* xb)
		Pa = Pa, lam :* SD
	}
	return((Pp, Ps, Pa))
}
// H = d sum psi / d phi' by central differences
real matrix _hess(real rowvector th, real matrix Sg, real rowvector al,
	real matrix W, real matrix LP, real colvector LX, real matrix Z,
	real scalar a0, real scalar qd, real matrix D, real matrix SD1,
	real matrix SD2, real rowvector cix)
{
	real matrix H
	real rowvector vs, tp, tm, vp, vm, ap, am
	real scalar P, L, nv, na, nt, k, h
	P = cols(th); L = rows(Sg); nv = L * (L + 1) / 2; na = cols(al); nt = P + nv + na
	vs = vech(Sg)'
	H = J(nt, nt, .)
	for (k = 1; k <= P; k++) {
		h = 1e-6 * max((1, abs(th[k])))
		tp = th; tp[k] = th[k] + h
		tm = th; tm[k] = th[k] - h
		H[., k] = (colsum(_psis(tp, Sg, al, W, LP, LX, Z, a0, qd, D, SD1, SD2, cix)) -
			colsum(_psis(tm, Sg, al, W, LP, LX, Z, a0, qd, D, SD1, SD2, cix)))' :/ (2 * h)
	}
	for (k = 1; k <= nv; k++) {
		h = 1e-6 * max((abs(vs[k]), 1e-4))
		vp = vs; vp[k] = vs[k] + h
		vm = vs; vm[k] = vs[k] - h
		H[., P + k] = (colsum(_psis(th, _unvech(vp, L), al, W, LP, LX, Z, a0, qd, D, SD1, SD2, cix)) -
			colsum(_psis(th, _unvech(vm, L), al, W, LP, LX, Z, a0, qd, D, SD1, SD2, cix)))' :/ (2 * h)
	}
	for (k = 1; k <= na; k++) {
		h = 1e-6 * max((1, abs(al[k])))
		ap = al; ap[k] = al[k] + h
		am = al; am[k] = al[k] - h
		H[., P + nv + k] = (colsum(_psis(th, Sg, ap, W, LP, LX, Z, a0, qd, D, SD1, SD2, cix)) -
			colsum(_psis(th, Sg, am, W, LP, LX, Z, a0, qd, D, SD1, SD2, cix)))' :/ (2 * h)
	}
	return(H)
}
// A = sum G' Sigma^-1 G, the Gauss-Newton information
real matrix _gnA(real rowvector th, real matrix Sg, real matrix W, real matrix LP,
	real colvector LX, real matrix Z, real scalar a0, real scalar qd, real matrix D,
	| real matrix PH, real matrix PHI, real rowvector cix)
{
	real matrix U, S, A, Hi
	pointer(real matrix) rowvector G
	real scalar i, j, L, N
	if (args() > 9) (void) _eq_model(th, W, LP, LX, Z, a0, qd, D, U, G, 1, PH, PHI, cix)
	else (void) _eq_model(th, W, LP, LX, Z, a0, qd, D, U, G, 1)
	L = cols(U); N = rows(U)
	S = invsym(Sg)
	A = J(cols(th), cols(th), 0)
	for (i = 1; i <= L; i++) {
		Hi = J(N, cols(th), 0)
		for (j = 1; j <= L; j++) Hi = Hi + S[i, j] :* *G[j]
		A = A + cross(*G[i], Hi)
	}
	return(A)
}
// one line of the report and of out/decomposition.csv
void _line(string scalar csv, string scalar tag, string scalar blk, string rowvector lab,
	real matrix R)
{
	real scalar k, fh
	fh = fopen(csv, "a")
	for (k = 1; k <= cols(lab); k++) {
		printf("{txt}  %-22s %-7s {res}%4.0f %8.4f %8.4f\n", blk, lab[k], rows(R), min(R[., k]), max(R[., k]))
		fput(fh, sprintf("%s,%s,%s,%g,%12.6f,%12.6f", tag, blk, lab[k], rows(R), min(R[., k]), max(R[., k])))
	}
	fclose(fh)
}
// Poi: the coefficients
void _do_poi(real matrix W, real matrix LP, real colvector LX, real matrix Z,
	real matrix UA, string rowvector GNn, real rowvector GNs)
{
	real matrix D, Sg, PS0, H, A, Pt, IFg, IFo, IFf, R
	real rowvector th, seg, seo, sef, sbf
	string rowvector nm
	real scalar P, M, qd, a0, N, nb, c, k
	th = st_matrix("TH"); P = cols(th); M = cols(W); qd = 1; a0 = strtoreal(st_local("a0"))
	N  = rows(W); D = _eq_dmat(P, M, 0, qd); Sg = st_matrix("SIG")
	nb = strtoreal(st_local("nb"))
	PS0 = _psis(th, Sg, J(1, 0, .), W, LP, LX, Z, a0, qd, D, J(0, 0, .), J(0, 0, .), J(1, 0, .))
	H   = _hess(th, Sg, J(1, 0, .), W, LP, LX, Z, a0, qd, D, J(0, 0, .), J(0, 0, .), J(1, 0, .))
	A   = _gnA(th, Sg, W, LP, LX, Z, a0, qd, D)
	Pt  = PS0[., 1..P]
	IFg = Pt * invsym(A)
	IFo = -Pt * luinv(H[1..P, 1..P])'
	IFf = -PS0 * luinv(H)'
	IFf = IFf[., 1..P]
	c   = N / (N - 1)
	seg = sqrt(diagonal(c :* cross(IFg * D', IFg * D')))'
	seo = sqrt(diagonal(c :* cross(IFo * D', IFo * D')))'
	sef = sqrt(diagonal(c :* cross(IFf * D', IFf * D')))'
	sbf = sqrt(rows(UA) / (rows(UA) - 1) :* colsum(UA[., 1..nb]:^2))
	printf("{txt}\npoi: the Gauss-Newton standard errors computed here / those of equaids 1.1.0: max |ratio - 1| = %9.2e\n",
		max(abs(seg :/ GNs[1..nb] :- 1)))
	nm = st_matrixcolstripe("EB")[., 2]'
	R  = (seg :/ sbf)', (seo :/ sbf)', (sef :/ sbf)'
	printf("\n{txt}%-14s %10s %10s %10s   (SE / brute-force SE)\n", "", "GN", "Hobs", "full")
	for (k = 1; k <= nb; k++) printf("{txt}%-14s {res}%10.4f %10.4f %10.4f\n", nm[k], R[k, 1], R[k, 2], R[k, 3])
	printf("\n{txt}  %-22s %-7s %4s %8s %8s\n", "block", "SE", "n", "min", "max")
	_line(st_local("BF") + "/out/decomposition.csv", "poi", "coefficients", ("GN", "Hobs", "full"), R)
}
// the non-buyers: psi = (theta, alpha)
void _do_sel(real matrix W, real matrix LP, real colvector LX, real matrix Z, real matrix Q,
	real matrix UA, string rowvector GNn, real rowvector GNs)
{
	real matrix D, Sg, SD1, SD2, PS0, H, IF, IFf, IFa, PH, PHI, A, IFp, R, zmask
	real rowvector cix, th, psi, al, ia, sef, sep, sbf, seg
	real scalar M, qd, a0, K, C, Pt, P, N, L, nv, na, nt, nb, c, np, k
	string scalar tag, csv
	tag = st_local("tag"); csv = st_local("BF") + "/out/decomposition.csv"
	cix = (1, 2); zmask = (1, 0 \ 0, 1)
	M = 3; qd = 0; a0 = 0; K = 1; C = 2
	th = st_matrix("TH"); Pt = cols(th); P = Pt - C
	D  = _eq_dmat(P, M, K, qd); Sg = st_matrix("SIG")
	psi = st_matrix("PSI"); al = psi[|Pt + 1 \ cols(psi)|]
	N  = rows(W); L = M - 1; nv = L * (L + 1) / 2; na = cols(al); nt = Pt + nv + na
	nb = strtoreal(st_local("nb"))
	SD1 = _eq_seldesign(LP, LX, Z, Q, zmask[1, .])
	SD2 = _eq_seldesign(LP, LX, Z, Q, zmask[2, .])
	PS0 = _psis(th, Sg, al, W, LP, LX, Z, a0, qd, D, SD1, SD2, cix)
	H   = _hess(th, Sg, al, W, LP, LX, Z, a0, qd, D, SD1, SD2, cix)
	IF  = -PS0 * luinv(H)'
	IFf = IF[., (1..Pt, Pt+nv+1..nt)]
	// the probits with the observed Hessian, the system as in 1.1.0
	ia  = (Pt+nv+1..nt)
	IFa = -PS0[., ia] * luinv(H[ia, ia])'
	PH = .; PHI = .
	_phphi(al, SD1, SD2, PH, PHI)
	A   = _gnA(th, Sg, W, LP, LX, Z, a0, qd, D, PH, PHI, cix)
	IFp = (PS0[., 1..Pt] + IFa * H[1..Pt, ia]') * invsym(A)
	IFp = IFp, IFa
	c   = N / (N - 1)
	sef = sqrt(diagonal(c :* cross(IFf, IFf)))'
	sep = sqrt(diagonal(c :* cross(IFp, IFp)))'
	np  = cols(psi)
	sbf = sqrt(rows(UA) / (rows(UA) - 1) :* colsum(UA[., nb+1..nb+np]:^2))
	seg = select(GNs, substr(GNn, 1, 3) :== "psi")
	R   = (seg :/ sbf)', (sep :/ sbf)', (sef :/ sbf)'
	printf("\n{txt}%s: %g households, theta %g, vech Sigma %g, alpha %g\n", tag, N, Pt, nv, na)
	printf("{txt}%-8s %10s %10s %10s   (SE / brute-force SE)\n", "psi", "GN", "probit", "full")
	for (k = 1; k <= np; k++) printf("{txt}%-8s {res}%10.4f %10.4f %10.4f\n", "psi" + strofreal(k), R[k, 1], R[k, 2], R[k, 3])
	printf("\n{txt}  %-22s %-7s %4s %8s %8s\n", "block", "SE", "n", "min", "max")
	_line(csv, tag, "system (theta)", ("GN", "probit", "full"), R[1..Pt, .])
	_line(csv, tag, "probits (alpha)", ("GN", "probit", "full"), R[Pt+1..np, .])
}
end
if `sel' mata: _do_sel(W, LP, LX, Z, Q, UA, GNn, GNs)
else     mata: _do_poi(W, LP, LX, Z, UA, GNn, GNs)

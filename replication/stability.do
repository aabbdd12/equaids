* stability.do -- Table tab:stab (Section 7.3, D7 of equaidsdiag): each
* demographic left out in turn, on the three data sets of the note, and the
* numbers of the paragraph after the table (rho, its z, the smallest m0(z)).
*   food  Poi's data have no demographic: two are built, zc moving with the
*         relative price of good 1 and zn pure noise (seed 1);
*   mex   Mexican cereals, hhsize and isMale, sampling weights;
*   lr    the data of Lecocq and Robin (2015), nbpers and rural
*         (get_lr_data.do downloads them).
* QUAIDS, default alpha_0.  For each demographic: the changes marked (|z|
* above the Bonferroni critical value) out of 2M, the largest |z|, its
* correlations with ln x and with the relative log prices (D3), and in the
* full model its rho, the z of rho and the smallest m0(z); the time of D7.
* Writes out/stability.csv.
*
*   cd <repository>/replication
*   do stability.do            (the three data sets, about 4 minutes)
*   do stability.do mex        (one of them)
args data
version 14.2
set more off
local ROOT = subinstr("`c(pwd)'", "\", "/", .) + "/.."
capture confirm file "`ROOT'/src/equaids.ado"
if _rc {
    display as error "run this script from the replication/ directory"
    exit 601
}
if "`data'" == "" {
    capture mkdir out
    tempname fh
    file open `fh' using "out/stability.csv", write text replace
    file write `fh' "data,demographic,marked,changes,zmax,corr_x,corr_p,rho,z_rho,m0_min,zcrit,seconds" _n
    file close `fh'
    foreach d in food mex lr {
        do stability.do `d'
    }
    exit
}

clear all
adopath ++ "`ROOT'/src"
if "`data'" == "food" {
    webuse food, clear
    set seed 1
    gen double zc = p1 / p4 + .1 * runiform()
    gen double zn = runiform()
    local CMD "w1-w4, prices(p1-p4) expenditure(expfd) demographics(zc zn)"
}
else if "`data'" == "mex" {
    use "`ROOT'/examples/mexico_2014_cereals.dta", clear
    local CMD "wcorn wwheat wrice wother wcomp [pw=sweight], prices(pcorn pwheat price pother pcomp) expenditure(hh_current_inc) demographics(hhsize isMale)"
}
else {
    do get_lr_data.do
    local CMD "s1-s7, prices(p1-p7) expenditure(somtot) demographics(nbpers rural)"
}

timer clear 1
timer on 1
equaidsdiag `CMD' stability
timer off 1
quietly timer list 1
local secs = r(t1)
matrix S = r(stability)
matrix D = r(demo)
local zcrit = r(stab_zcrit)
local demos : rownames D

* the full model: rho, its z and the smallest m0(z)
quietly equaids `CMD' notable nolog
matrix B = e(b)
matrix V = e(V)
local m0min = e(m0_min)

display as text _n "`data': D7 in " %5.1f `secs' " s, Bonferroni critical value " %4.2f `zcrit'
display as text %-10s "" %10s "marked" %10s "max |z|" %10s "corr ln x" %10s "corr p" ///
    %10s "rho" %8s "z" %10s "min m0"
mata:
S = st_matrix("S"); D = st_matrix("D"); B = st_matrix("B"); V = st_matrix("V")
cs = st_matrixcolstripe("B"); dn = st_matrixcolstripe("D")[., 2]'
dr = st_matrixrowstripe("D")[., 2]'
zc = strtoreal(st_local("zcrit")); K = rows(D); M = rows(S) / K
fh = fopen("out/stability.csv", "a")
for (k = 1; k <= K; k++) {
    Z  = abs(S[|(k - 1) * M + 1, 3 \ k * M, 3|]), abs(S[|(k - 1) * M + 1, 6 \ k * M, 6|])
    j  = selectindex((cs[., 1] :== "rho") :& (cs[., 2] :== "rho_" + dr[k]))
    rho = B[j]; zr = rho / sqrt(V[j, j])
    cx = D[k, selectindex(dn :== "corr_x")]; cp = D[k, selectindex(dn :== "corr_p")]
    printf("{txt}%-10s {res}%6.0f of %2.0f %10.2f %10.2f %10.2f %10.3f %8.1f %10.3f\n", dr[k],
        sum(Z :> zc), 2 * M, max(Z), cx, cp, rho, zr, strtoreal(st_local("m0min")))
    fput(fh, sprintf("%s,%s,%g,%g,%9.4f,%9.4f,%9.4f,%12.6f,%9.3f,%9.4f,%6.3f,%8.1f", st_local("data"),
        dr[k], sum(Z :> zc), 2 * M, max(Z), cx, cp, rho, zr, strtoreal(st_local("m0min")), zc,
        strtoreal(st_local("secs"))))
}
fclose(fh)
end

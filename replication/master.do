* master.do -- runs the fast replication scripts of the technical note (about
* ten minutes).  The long computations (bootstrap replications, brute force,
* Monte Carlo, timings) are shipped as raw results: their summaries run here.
*
*   cd <repository>/replication
*   do master.do
version 14.2
set more off
local here "`c(pwd)'"
foreach s in types_elasticities remark_alpha0 pimpute_effect dgp_describe ///
    diagnostics stability example_design example_engel selection {
    display as text _n "{hline 78}" _n "`s'.do" _n "{hline 78}"
    do `s'.do
    quietly cd "`here'"
}
do mc_nonbuyers.do combine 4
quietly cd "`here'/bootstrap"
do boot_summary.do
quietly cd "`here'/bruteforce"
do numif_summary.do
do decomposition.do
quietly cd "`here'/montecarlo"
do mcs_summary.do
quietly cd "`here'"

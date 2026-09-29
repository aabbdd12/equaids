*! spec_selimp.do -- spec_sel.do with the prices of the non-buyers missing
*! and filled by pimpute(grp) (see spec_sel.do).
global BF_IMP 1
run "spec_sel.do"

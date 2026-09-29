*! spec_mean.do -- the brute-force specification of the simplest estimator,
*! the mean of x over 2,000 households (for bf_worker.do / bf_combine.do).
capture program drop _bf_setup
program define _bf_setup
	clear
	set seed 20260928
	quietly set obs 2000
	quietly gen double x = exp(rnormal(3, 0.5))
	quietly gen double w = 1
end
capture program drop _bf_theta
program define _bf_theta
	quietly summarize x [iw = w], meanonly
	mata: __bf_theta = st_numscalar("r(mean)")
end
capture program drop _bf_analytic
program define _bf_analytic
	quietly summarize x
	mata: __bf_se = st_numscalar("r(sd)") / sqrt(st_numscalar("r(N)")); __bf_names = "mean of x"
end

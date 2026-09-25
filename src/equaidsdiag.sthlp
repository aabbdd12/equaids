{smcl}
{* *! version 0.1.0  25sep2026}{...}
{* the threshold of D4 is provisional: to be set after the audit on three data sets (audit/diag/)}{...}
{vieweralsosee "equaids" "help equaids"}{...}
{viewerjumpto "Syntax" "equaidsdiag##syntax"}{...}
{viewerjumpto "Description" "equaidsdiag##description"}{...}
{viewerjumpto "What it checks" "equaidsdiag##checks"}{...}
{viewerjumpto "Stored results" "equaidsdiag##results"}{...}
{viewerjumpto "Example" "equaidsdiag##example"}{...}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col:{cmd:equaidsdiag} {hline 2}}Diagnose an AIDS/QUAIDS specification before estimating it{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 20 2}
{cmd:equaidsdiag} {it:shares} {ifin} [{it:weight}]{cmd:,} {it:equaids_options}
[{opt sens:itivity} {opt a0l:ist(numlist)}]

{p 4 4 2}
The syntax is that of {helpb equaids}: the same shares, {opt prices()} or
{opt lnprices()}, {opt expenditure()} or {opt lnexpenditure()},
{opt demographics()}, {opt noquadratic}, {opt anot()} and weights. The options
that belong to the estimator alone ({opt vce()}, {opt elasticities()},
{opt dec()} and the others) are accepted, ignored, and listed.

{synoptset 22}{...}
{synopthdr}
{synoptline}
{synopt:{opt sens:itivity}}estimate the model at several values of alpha_0 (D6){p_end}
{synopt:{opt a0l:ist(numlist)}}the values of alpha_0; default: the smallest log expenditure minus 0.1, and 1, 2 and 4 below; implies {opt sensitivity}{p_end}
{synopt:{opt stab:ility}}estimate the model again without each demographic in turn (D7){p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:equaidsdiag} reports whether a specification is fit to be estimated, and
{bf:it does not estimate it} (except under {opt sensitivity} and {opt stability}). Sections D0 to
D5 are computed from the data and from the regressors at the starting point:
alpha at the mean shares and the other parameters at zero, so that the price
index is the Stone index plus alpha_0, and deflated expenditure is
l = ln x - alpha_0 - sum_k wbar_k ln p_k. A diagnostic that needed the model
to converge would be silent exactly when it does not.


{marker checks}{...}
{title:What it checks}

{dlgtab:Households lost}

{p 4 4 2}
The households of the {cmd:if}/{cmd:in} sample dropped for a missing share,
price, expenditure or demographic, or a nonpositive weight, variable by
variable. On surveys where prices are unit values, the missing prices of the
non-buyers can remove a large part of the sample.

{dlgtab:D0. Budget shares}

{p 4 4 2}
Mean and aggregate shares (the share of each good in total expenditure),
zero shares, shares outside [0,1], shares that do not sum to one. Many zeros
mean censoring, which the estimator does not model, and prices of the
non-buyers that have been filled in.

{dlgtab:D1. Relative prices}

{p 4 4 2}
The standard deviation of each relative log price ln(p_k/p_M) (below 0.05, the
price parameters are weakly identified), the largest correlation between two
relative prices, and extreme prices (more than 5 robust standard deviations
from the median).

{dlgtab:D2. Expenditure and alpha_0}

{p 4 4 2}
The range of ln x, alpha_0 (the rule of {helpb equaids} or {opt anot()}), and
the households whose deflated expenditure l is negative: when alpha_0 is above
the log expenditure of most households, the price index exceeds their
expenditure.

{dlgtab:D3. Demographic variables}

{p 4 4 2}
Distinct values, range, mean, standard deviation; for a binary variable, the
share of the rarer modality and the effective size N p (1-p); integer
variables with few values (categorical variables should enter as
indicators); negative values, which can bring Ray's m0 = 1 + rho'z near 0; the largest
correlation with a relative log price and the correlation with ln x (a note
above 0.3 in absolute value: omitted, such a demographic biases the price or
the expenditure elasticities).

{dlgtab:D4. Conditioning of the regressors at the starting point}

{p 4 4 2}
The condition indexes of Belsley, Kuh and Welsch (1980) of the columns of the
Jacobian at the starting point: the constant, the relative log prices, l,
l^2 (QUAIDS), the demographics and their products with l; columns
uncentered and scaled to unit length, so that the constant, and with it the
level of l, is taken into account. A component with an index above 30 on
which two or more columns carry more than half of their variance is a near
dependency: noted between 30 and 100 (moderate), a warning above 100
(strong). The two levels were calibrated on three data sets (Poi's food
data, Mexican cereals, the data of Lecocq and Robin): with the default
alpha_0 the index of the quadratic block is 23 to 36 and the standard errors
agree with the bootstrap; every ill-conditioned estimate had an index above
100. The index of the block (1, l, l^2) alone is reported: when
alpha_0 lies far from the log expenditures, l varies little relative to its
level and l^2 is almost a linear function of l, so that the quadratic
coefficients are weakly identified and the sandwich standard errors of the
coefficients unreliable. On Poi's data with his alpha_0 = 10 this index is
331, the information matrix of the estimate nearly singular, and the
bootstrap standard deviations of the coefficients two to four times the
sandwich standard errors; with the default alpha_0 it is 23.

{dlgtab:D5. Small goods}

{p 4 4 2}
Goods with less than 1% of total expenditure: QUAIDS can predict their
shares near zero or negative, which makes the mean of the household
elasticities unstable (not the aggregate elasticities, the default of
{helpb equaids}).

{dlgtab:D6. Sensitivity to alpha_0 (option sensitivity)}

{p 4 4 2}
The model estimated at each value of alpha_0: convergence, iterations, log
likelihood, reciprocal condition number of the information matrix, and the
aggregate expenditure elasticities.

{dlgtab:D7. Stability to each demographic (option stability)}

{p 4 4 2}
The model is estimated with all the demographics, then again without each of
them in turn, on the same sample and at the same alpha_0. For every good, the
table gives the change dE of the aggregate expenditure elasticity and of the
aggregate own-price elasticity when the demographic is left out, its robust
standard error, and z = dE / s.e. The standard error is that of the
difference of the two estimators: both are estimated on the same households,
so that the difference of their influence functions, household by household,
is the influence function of dE. A change is marked when |z| exceeds the
Bonferroni critical value for the 2M changes of a demographic,
invnormal(1 - 0.05/(4M)) (2.73 with four goods). The standard errors are
robust whatever the design: D7 asks whether the elasticities move, not how
precise they are.

{p 4 4 2}
A marked change means that the elasticities depend on the demographic. Most
often it belongs in the model, through its correlation with expenditure or
prices (D3): leaving it out moves its effect onto the expenditure and price
terms. On the data of Lecocq and Robin (seven goods, 25,776 households),
household size is correlated 0.30 with ln x and only 0.19 with the prices;
leaving it out moves the expenditure elasticities by up to z = 21, with rho
well identified (z = 6.8) and 1 + rho'z at least 1.33. Check also its
rho: when 1 + rho'z nears zero for a few households, the fit is driven by
them and the change is large as well; on Poi's food data a uniform noise
variable reaches rho = -1.1 (z = -10), min(1 + rho'z) = 0.02, and moves one
own-price elasticity by z = -3.1. D7 and D6 are the only sections that
estimate.


{marker results}{...}
{title:Stored results}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(N)}}number of observations used{p_end}
{synopt:{cmd:r(N_lost)}}households lost to missing values{p_end}
{synopt:{cmd:r(N_warn)}}number of warnings{p_end}
{synopt:{cmd:r(anot)}}alpha_0{p_end}
{synopt:{cmd:r(cond_max)}}largest condition index (D4){p_end}
{synopt:{cmd:r(cond_quad)}}condition index of (1, l, l^2){p_end}
{synopt:{cmd:r(n_l0neg)}}households with l <= 0 at the start{p_end}
{synopt:{cmd:r(stab_zmax)}}largest |z| of D7{p_end}
{synopt:{cmd:r(stab_nsig)}}number of marked changes in D7{p_end}
{synopt:{cmd:r(stab_zcrit)}}critical value of D7 (Bonferroni){p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:r(shares)}}D0 table{p_end}
{synopt:{cmd:r(bkw)}}condition indexes and variance-decomposition proportions{p_end}
{synopt:{cmd:r(demo)}}D3 table{p_end}
{synopt:{cmd:r(sensitivity)}}D6 table{p_end}
{synopt:{cmd:r(stability)}}D7 table: dE_x, se_x, z_x, dE_ii, se_ii, z_ii, one row per demographic and good{p_end}


{marker example}{...}
{title:Example}

{phang2}{cmd:. webuse food}{p_end}
{phang2}{cmd:. equaidsdiag w1-w4, prices(p1-p4) expenditure(expfd)}{p_end}
{phang2}{cmd:. equaidsdiag w1-w4, prices(p1-p4) expenditure(expfd) anot(10)}{p_end}
{phang2}{cmd:. equaidsdiag w1-w4, prices(p1-p4) expenditure(expfd) sensitivity}{p_end}

{pstd}With demographic variables {it:z1} and {it:z2}:{p_end}
{phang2}{cmd:. equaidsdiag} {it:shares}{cmd:, prices(}{it:prices}{cmd:) expenditure(}{it:x}{cmd:) demographics(}{it:z1 z2}{cmd:) stability}{p_end}


{title:Reference}

{phang}
Belsley, D. A., E. Kuh, and R. E. Welsch. 1980. {it:Regression Diagnostics:
Identifying Influential Data and Sources of Collinearity}. New York: Wiley.


{title:Author}

{pstd}Abdelkrim Araar, Universit{c e'} Laval / PEP, aabd@ecn.ulaval.ca{p_end}

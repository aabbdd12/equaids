{smcl}
{* *! version 1.0.0  25sep2026}{...}
{vieweralsosee "equaidsdiag" "help equaidsdiag"}{...}
{vieweralsosee "[R] demandsys" "help demandsys"}{...}
{viewerjumpto "Syntax" "equaids##syntax"}{...}
{viewerjumpto "Description" "equaids##description"}{...}
{viewerjumpto "Options" "equaids##options"}{...}
{viewerjumpto "Remarks" "equaids##remarks"}{...}
{viewerjumpto "Postestimation" "equaids##postest"}{...}
{viewerjumpto "Stored results" "equaids##results"}{...}
{viewerjumpto "Examples" "equaids##examples"}{...}
{viewerjumpto "References" "equaids##references"}{...}
{title:Title}

{p2colset 5 16 18 2}{...}
{p2col:{cmd:equaids} {hline 2}}AIDS and QUAIDS demand systems, with aggregate elasticities and survey-design inference{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:equaids} {it:shares} {ifin} [{it:{help equaids##weight:weight}}]{cmd:,}
{c -(}{opt pr:ices(varlist)} | {opt lnpr:ices(varlist)}{c )-}
{c -(}{opt exp:enditure(varname)} | {opt lnexp:enditure(varname)}{c )-}
[{it:options}]

{p 8 8 2}
{it:shares} are the budget shares of the {it:M} goods (at least three), which
must sum to one; the prices are given in the same order.

{synoptset 30 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{p2coldent:* {opt pr:ices(varlist)}}prices of the goods; or {opt lnpr:ices()}, their logarithms{p_end}
{p2coldent:* {opt exp:enditure(varname)}}total expenditure; or {opt lnexp:enditure()}, its logarithm{p_end}
{synopt:{opt noqu:adratic}}estimate the AIDS instead of the QUAIDS{p_end}
{synopt:{opt demo:graphics(varlist)}}demographic variables, entering by Ray's scaling{p_end}
{synopt:{opt anot(#)}}impose the constant alpha_0 of the price index; default: the smallest log expenditure minus 0.1{p_end}

{syntab:Variance}
{synopt:{opt vce(robust)}}robust; the default{p_end}
{synopt:{opt vce(cluster} {it:clustvar}{cmd:)}}clustered{p_end}
{synopt:{opt vce(svy)}}survey design declared by {helpb svyset}: strata, PSUs, finite-population correction, pweight{p_end}
{synopt:{opt vce(conventional)}}conventional (coefficients){p_end}
{synopt:{opt l:evel(#)}}confidence level; default {cmd:level(95)}{p_end}

{syntab:Elasticities and reporting}
{synopt:{opt elas:ticities(type)}}summary of the household elasticities: {cmd:aggregate} (the default), {cmd:means} or {cmd:household}{p_end}
{synopt:{opt compens:ated}}add the compensated (Hicksian) price elasticities{p_end}
{synopt:{opt checks}}show the aggregation identities{p_end}
{synopt:{opt det:ail}}{opt compensated} and {opt checks}{p_end}
{synopt:{opt noelastse}}do not compute the standard errors of the elasticities{p_end}
{synopt:{opt sn:ames(namelist)}}short names of the goods in the tables{p_end}
{synopt:{opt dec(#)}}decimals displayed; default 4{p_end}
{synopt:{opt dislas(0|1)}}show the last good; default 1{p_end}
{synopt:{opt dregres(0|1)}}display the coefficients; default 0{p_end}
{synopt:{opt st:ars}}significance stars on the estimates{p_end}
{synopt:{opt saveres(filename)}}write the tables to a file: {cmd:.docx}, {cmd:.tex}, {cmd:.xlsx}, {cmd:.csv} or {cmd:.md}{p_end}
{synopt:{opt notab:le}}do not display the tables{p_end}

{syntab:Estimation}
{synopt:{opt tol:erance(#)}}tolerance on the Newton decrement; default 1e-6{p_end}
{synopt:{opt iter:ate(#)}}maximum number of iterations; default 300{p_end}
{synopt:{opt from(matname)}}starting values: a row vector of free parameters, as {cmd:e(b_free)}{p_end}
{synopt:{opt nolog}}suppress the iteration log{p_end}
{synoptline}
{p 4 6 2}* one of each pair is required.{p_end}

{marker weight}{...}
{p 4 6 2}{opt aweight}s, {opt fweight}s, {opt pweight}s and {opt iweight}s are
allowed; see {help weight}. Analytic and sampling weights are normalized to sum
to the number of observations. With {cmd:vce(svy)} the weight comes from
{helpb svyset} and no weight may be given.{p_end}

{p 4 6 2}
The display options ({opt elasticities()}, {opt compensated}, {opt checks},
{opt dec()}, {opt dislas()}, {opt dregres()}, {opt stars}, {opt saveres()},
{opt notable}) can be given again on replay, {cmd:equaids} typed alone, without
re-estimating.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:equaids} estimates the almost ideal demand system (AIDS) of Deaton and
Muellbauer (1980) and its quadratic extension (QUAIDS) by Banks, Blundell and
Lewbel (1997), with demographic variables entering by the scaling of Ray
(1983), as in Poi (2012). For household {it:h} and good {it:i},

{p 8 8 2}
w_i = alpha_i + sum_j gamma_ij ln p_j + (beta_i + eta_i'z) l + lambda_i l^2 / (b(p) c(p,z)),

{p 8 8 2}
l = ln x - ln m0(z) - ln a(p),  ln a(p) = alpha_0 + sum_k alpha_k ln p_k + 1/2 sum_k sum_l gamma_kl ln p_k ln p_l,

{pstd}
with b(p) = prod_k p_k^beta_k, c(p,z) = prod_k p_k^(eta_k'z) and m0(z) = 1 +
rho'z. Adding-up, homogeneity and symmetry are imposed. AIDS is the case
lambda = 0.

{pstd}
The estimator is iterated feasible generalized nonlinear least squares, the
Gaussian quasi-maximum-likelihood estimator of the system, that of Poi (2012)
and of {helpb demandsys}. {cmd:equaids} computes it in Mata by Gauss-Newton
steps with an analytic Jacobian, and reproduces the estimates of these
commands under the same model and alpha_0.

{pstd}
The default summary of the elasticities is the {it:aggregate} elasticity: the
elasticity of the total demand of the population, the household elasticities
weighted by each household's expenditure on the good. Its standard errors are
analytic and include the sampling of the households as well as the estimation
of the coefficients. All the variances (robust, clustered, survey design) are
built from the same influence functions. See {help equaids##remarks:Remarks}.

{pstd}
{helpb equaidsdiag} diagnoses a specification before estimating it.


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{opt prices(varlist)} or {opt lnprices(varlist)} give the prices of the goods,
in the order of the shares, as levels or as logs; {opt expenditure(varname)}
or {opt lnexpenditure(varname)} give total expenditure (on the goods of the
system).

{phang}
{opt noquadratic} estimates the AIDS (lambda = 0).

{phang}
{opt demographics(varlist)} adds demographic variables by Ray's scaling: they
shift the intercepts of l (through m0 = 1 + rho'z) and the slopes (eta). The
scaling requires m0 > 0: counts and indicator variables (z >= 0) keep m0 away
from 0; a categorical variable should enter as indicators.

{phang}
{opt anot(#)} imposes the constant alpha_0 of the price index, which is not
identified with the other parameters. By default alpha_0 is the smallest log
expenditure of the estimation sample minus 0.1, so that deflated expenditure
is positive (Banks, Blundell and Lewbel 1997), a rule of the data recomputed
on every sample. A value far from the log expenditures makes l^2 almost a
linear function of l and the quadratic terms weakly identified; see
{helpb equaidsdiag}. Use {opt anot()} to reproduce results that impose a value
(Poi's example uses 10).

{dlgtab:Variance}

{phang}
{opt vce(robust)}, the default, sums the squares of the influence functions
of the households (factor N/(N-1)). {opt vce(cluster} {it:clustvar}{cmd:)}
sums them by cluster (factor G/(G-1)).

{phang}
{opt vce(svy)} reads the design from {helpb svyset}: the influence functions
are summed by PSU, and the PSU totals are combined within strata with the
finite-population correction (first-stage linearization, as Stata's
{cmd:svy}). The tests use the design degrees of freedom, PSUs minus strata.
Strata with a single PSU follow {cmd:svyset}'s {opt singleunit()}:
{cmd:certainty}, {cmd:scaled} and {cmd:centered} as in Stata; with
{cmd:missing}, the default, the model is estimated and the standard errors
are missing, with a note. A single PSU or cluster in the sample is refused.
The linearized variance is valid as the number of PSUs per stratum grows:
with few PSUs per stratum the intervals are slightly too narrow. In a
simulation of equaids with 5 to 16 PSUs per stratum, the 95% intervals of
the aggregate elasticities covered 92.5% to 94.4% of the time; with 10 to 32
PSUs per stratum, 93.1% to 94.8%.

{phang}
{opt vce(conventional)} gives the conventional variance of the coefficients
(no weights, homoskedastic errors); it cannot be combined with pweights.

{dlgtab:Elasticities and reporting}

{phang}
{opt elasticities(aggregate)}, the default, reports the aggregate
elasticities, with standard errors. {opt elasticities(means)} evaluates the
household elasticities at the weighted means of ln p, ln x and z, as Poi
(2012); {opt elasticities(household)} averages the household elasticities
with the sampling weights, as {helpb demandsys} (which does not weight the
average). These two have no standard errors. The household elasticities
divide by the predicted shares: {cmd:equaids} warns when some are near zero
or outside [0,1].

{phang}
{opt compensated} adds the table of the compensated price elasticities;
{opt checks} displays the aggregation identities of the aggregate
elasticities (Engel, Cournot, homogeneity, symmetry of the compensated
matrix), which hold exactly and are always checked; {opt detail} does both.

{phang}
{opt noelastse} skips the standard errors of the elasticities (the robust
variance of the coefficients is still computed).

{phang}
{opt snames()}, {opt dec()}, {opt dislas()}, {opt dregres()}, {opt stars},
{opt saveres()} and {opt notable} control the tables, as in {cmd:easi} and
{cmd:duvm}. {opt saveres()} writes all the tables to one file whose extension
gives the format.

{dlgtab:Estimation}

{phang}
{opt tolerance(#)} is the tolerance on the Newton decrement g'A^-1 g, the
squared distance to the optimum in the metric of the information: the
default, 1e-6, puts the estimate within 0.001 standard errors of the optimum,
in all directions jointly. The covariance of the residuals must be stable to
the same relative order.

{phang}
{opt iterate(#)} is the maximum number of iterations; default 300. The
iterations stop earlier, without convergence, when neither the log
likelihood nor the scaled gradient has progressed over 20 iterations; the
reasons are reported.

{phang}
{opt from(matname)} gives starting values for the free parameters (a row
vector such as {cmd:e(b_free)}); by default alpha is started at the mean
shares and the other parameters at zero.


{marker remarks}{...}
{title:Remarks}

{pstd}
{bf:The aggregate elasticity.} With x_h f_ih the expenditure of household h
on good i (f the predicted share), the aggregate expenditure elasticity is
E_i = 1 + sum_h w_h x_h mu_ih / sum_h w_h x_h f_ih, with mu_i = dw_i/d ln x;
the price elasticities are built in the same way. It is the elasticity of the
total demand of the population, what a simulation of a price or tax change on
aggregate demand needs; it divides by no household share, so that goods with
many small or zero shares do not make it unstable; and, as a ratio of weighted
totals, its inference under a survey design is standard. The mean of the
household elasticities, by contrast, divides by each predicted share: on
small goods a few households with shares near zero drive it.

{pstd}
{bf:Standard errors.} The influence function of an aggregate elasticity has
two terms, the sampling of the households (the summary is taken over a
sample) and the estimation of the coefficients; both are analytic. They agree
with the bootstrap of the households and with the Rao-Wu bootstrap of a
survey design; the design variance equals Stata's own linearization.

{pstd}
{bf:Diagnostics.} Before estimating, {cmd:equaids} refuses what is not
identified and notes near collinearity, weak variation of relative prices,
rare modalities of the demographics, zero shares and extreme prices. After
estimating, it explains a failure to converge and notes an optimum near the
boundary of Ray's scaling, households whose expenditure is below the
estimated cost a(p), and predicted shares outside [0,1].
{cmd:estat diagnostics} displays these diagnostics again after estimation;
{helpb equaidsdiag} gives the full diagnosis of a specification before it is
estimated.

{pstd}
{bf:Users of WELCOM's wquaids.} {opt hweight(w)} becomes {cmd:[pw=w]},
{opt model(2)} becomes {opt noquadratic}, and {opt xfil()} becomes
{opt saveres()}.


{marker postest}{...}
{title:Postestimation}

{p 8 15 2}
{cmd:estat diagnostics}{p_end}

{p 8 15 2}
{cmd:estat engel} {ifin} [{cmd:,} {it:engel_options}]{p_end}

{synoptset 26 tabbed}{...}
{synopthdr:engel_options}
{synoptline}
{synopt:{opt atm:eans}}log prices and demographics at their weighted means; the default{p_end}
{synopt:{opt asob:served}}the fitted shares of the households, smoothed{p_end}
{synopt:{opt obs:erved}}with {opt asobserved}: add the smoothed observed shares{p_end}
{synopt:{opt n(#)}}number of points of the grid; default {cmd:n(100)}{p_end}
{synopt:{opt bw:idth(#)}}bandwidth of the smoother; {opt asobserved} only{p_end}
{synopt:{opt trim(#)}}percent trimmed from each tail of the grid; default {cmd:trim(1)}{p_end}
{synopt:{opt l:evel(#)}}confidence level of the band{p_end}
{synopt:{opt noci}}omit the confidence band{p_end}
{synopt:{opt noturn}}do not mark the turning points{p_end}
{synopt:{opt lnx}}log expenditure on the horizontal axis instead of its percentiles{p_end}
{synopt:{opt data(filename)}}save the plotted curves as a dataset{p_end}
{synopt:{opt sav:ing(filename)}}save the graph{p_end}
{synopt:{opt nodraw}}compute but do not draw{p_end}
{synoptline}

{pstd}
{cmd:estat engel} traces the budget share of each good against total
expenditure, one panel per good, as {cmd:estat engel} after {cmd:easi} and
{cmd:duvm}; the other options are passed to {helpb graph combine}.

{pstd}
{bf:At the means} (the default), the model's share is evaluated over a grid of
log expenditure (the weighted percentiles of ln x, tails trimmed), with the
log prices and the demographics at their weighted means, the point of
{opt elasticities(means)}. It is the exact function of the estimate: linear
in ln x for AIDS, quadratic for QUAIDS. The band is the delta method with the
analytic Jacobian and the variance of the estimate, {cmd:e(V_free)}: robust,
by cluster or by design as estimated, with t on the design degrees of freedom
under {cmd:vce(svy)}. For QUAIDS, the turning point of each curve, where
d w_i / d l = beta_i + eta_i'z + 2 lambda_i l / (b c) = 0, is marked by a
vertical line when it falls inside the grid, and returned in
{cmd:r(turn)} (ln x and its percentile). A curve that turns inside the range
of the data, or that crosses zero for a small good, shows where the quadratic
form bends the shares.

{pstd}
{bf:As observed}, the fitted shares of the households, with their own prices
and demographics, are smoothed against ln x by a local linear smoother (the
bandwidth of {helpb lpoly} by default). With {opt observed}, the observed
shares are smoothed with the same bandwidth and drawn dashed: where the two
curves part, the functional form does not follow the data. The comparison is
fair only as observed, since the observed shares vary with prices and
demographics along ln x; {opt observed} is therefore refused at the means.

{pstd}
{cmd:estat engel} stores {cmd:r(n)} and, at the means, {cmd:r(turn)}; with
{opt asobserved}, {cmd:r(bwidth)}.


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:equaids} stores the following in {cmd:e()}:

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(ll)}}log likelihood{p_end}
{synopt:{cmd:e(anot)}}alpha_0{p_end}
{synopt:{cmd:e(ngoods)}, {cmd:e(ndemos)}}numbers of goods and of demographics{p_end}
{synopt:{cmd:e(converged)}, {cmd:e(iter)}}convergence (1 = yes) and iterations{p_end}
{synopt:{cmd:e(nrgrad)}, {cmd:e(sigdif)}}final Newton decrement and change of Sigma{p_end}
{synopt:{cmd:e(stalled)}}1 if the iterations stopped for lack of progress{p_end}
{synopt:{cmd:e(rcond)}}reciprocal condition number of the scaled information matrix{p_end}
{synopt:{cmd:e(m0_min)}, {cmd:e(n_lneg)}, {cmd:e(n_shout)}}diagnostics of the estimate{p_end}
{synopt:{cmd:e(N_clust)}}number of clusters ({cmd:vce(cluster)}){p_end}
{synopt:{cmd:e(N_strata)}, {cmd:e(N_psu)}, {cmd:e(df_r)}}design ({cmd:vce(svy)}){p_end}
{synopt:{cmd:e(chk_engel)}, ...}residuals of the aggregation identities{p_end}
{synopt:{cmd:e(time)}}execution time in seconds{p_end}

{p2col 5 22 26 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:equaids}{p_end}
{synopt:{cmd:e(model)}}{cmd:QUAIDS} or {cmd:AIDS}{p_end}
{synopt:{cmd:e(vce)}}{cmd:robust}, {cmd:cluster}, {cmd:svy} or {cmd:conventional}{p_end}
{synopt:{cmd:e(anot_rule)}}rule of alpha_0, or {cmd:user}{p_end}
{synopt:{cmd:e(elasticities)}}summary reported{p_end}
{synopt:{cmd:e(data_notes)}}notes of the data diagnostics{p_end}

{p2col 5 22 26 2: Matrices}{p_end}
{synopt:{cmd:e(b)}, {cmd:e(V)}}coefficients (Poi's parameterization) and their variance{p_end}
{synopt:{cmd:e(b_free)}, {cmd:e(V_free)}}free parameters and their variance{p_end}
{synopt:{cmd:e(Sigma)}}covariance of the residuals{p_end}
{synopt:{cmd:e(elas_x)}, {cmd:e(elas_u)}, {cmd:e(elas_c)}}aggregate expenditure, uncompensated and compensated elasticities (rows: goods, columns: prices){p_end}
{synopt:{cmd:e(V_elas_x)}, {cmd:e(V_elas_u)}, {cmd:e(V_elas_c)}}their variances{p_end}
{synopt:{cmd:e(elas_xm)}, ...}the same at the means ({cmd:m}) and household mean ({cmd:h}){p_end}
{synopt:{cmd:e(aggshare)}}aggregate budget shares{p_end}
{synopt:{cmd:e(vif)}, {cmd:e(demo_stats)}}data diagnostics{p_end}


{marker examples}{...}
{title:Examples}

{pstd}Poi's data (four food groups){p_end}
{phang2}{cmd:. webuse food}{p_end}
{phang2}{cmd:. equaids w1-w4, prices(p1-p4) expenditure(expfd) snames(meat fruitveg bread dairy)}{p_end}
{phang2}{cmd:. equaids, compensated checks stars}{p_end}
{phang2}{cmd:. equaids, elasticities(household)}{p_end}

{pstd}Survey design{p_end}
{phang2}{cmd:. svyset psu [pw=weight], strata(stratum)}{p_end}
{phang2}{cmd:. equaids w1-w5, prices(p1-p5) expenditure(x) demographics(hhsize) vce(svy)}{p_end}

{pstd}Engel curves{p_end}
{phang2}{cmd:. estat engel}{p_end}
{phang2}{cmd:. estat engel, lnx level(90) data(curves, replace)}{p_end}
{phang2}{cmd:. estat engel, asobserved observed}{p_end}

{pstd}Diagnose before estimating{p_end}
{phang2}{cmd:. equaidsdiag w1-w4, prices(p1-p4) expenditure(expfd) sensitivity}{p_end}


{marker references}{...}
{title:References}

{phang}
Araar, A. 2026. Estimating AIDS and QUAIDS demand systems with survey data:
the equaids Stata module. Technical note, Zenodo.
{browse "https://doi.org/10.5281/zenodo.22959991":doi:10.5281/zenodo.22959991}.

{phang}
Banks, J., R. Blundell, and A. Lewbel. 1997. Quadratic Engel curves and
consumer demand. {it:Review of Economics and Statistics} 79: 527-539.

{phang}
Deaton, A., and J. Muellbauer. 1980. An almost ideal demand system.
{it:American Economic Review} 70: 312-326.

{phang}
Poi, B. P. 2012. Easy demand-system estimation with quaids. {it:Stata Journal}
12: 433-446.

{phang}
Ray, R. 1983. Measuring the costs of children. {it:Journal of Public
Economics} 22: 89-102.


{title:Author}

{pstd}Abdelkrim Araar, Universit{c e'} Laval / PEP, aabd@ecn.ulaval.ca{p_end}

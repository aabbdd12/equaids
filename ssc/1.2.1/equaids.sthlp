{smcl}
{* *! version 1.2.1  01oct2026}{...}
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
{p2col:{cmd:equaids} {hline 2}}AIDS and QUAIDS demand systems, with elasticities of the households, the individuals or the market, and survey-design inference{p_end}
{p2colreset}{...}

{p 4 4 2}{txt}Package {cmd:equaids}, version {res}1.2.1{txt} (01/10/2026) {c |} Stata {res}14.2{txt} or later {c |} first release {res}1.0.0{txt} (25/09/2026){p_end}


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

{syntab:Missing prices and non-buyers}
{synopt:{opt pimp:ute(varlist)}}fill a missing price by the mean log price of the same group, the groups tried in order{p_end}
{synopt:{opt sel:ection}}correct for the households that do not buy (Shonkwiler and Yen 1999){p_end}
{synopt:{opt selg:oods(namelist)}}the shares corrected; default: those with zeros, the last good excepted; implies {opt selection}{p_end}
{synopt:{cmd:selvars(}[{it:good}{cmd::}] {it:varlist} [{cmd:;} ...]{cmd:)}}variables of the probits only, for every corrected good or for one good; implies {opt selection}{p_end}

{syntab:Variance}
{synopt:{opt vce(robust)}}robust; the default{p_end}
{synopt:{cmd:vce(cluster} {it:clustvar}{cmd:)}}clustered{p_end}
{synopt:{opt vce(svy)}}survey design declared by {helpb svyset}: strata, PSUs, finite-population correction, pweight{p_end}
{synopt:{opt vce(conventional)}}conventional (coefficients){p_end}
{synopt:{cmd:vce(bootstrap} [{cmd:,} {it:{help equaids##bootopts:boot_opts}}]{cmd:)}}bootstrap of the whole procedure{p_end}
{synopt:{opt l:evel(#)}}confidence level; default {cmd:level(95)}{p_end}

{syntab:Elasticities and reporting}
{synopt:{opt elas:ticities(type)}}{cmd:households} (the default), {cmd:individuals}, {cmd:market}, {cmd:reference} or {cmd:hhmean}{p_end}
{synopt:{opt hhs:ize(varname)}}size of the household; with it, the elasticities are those of the individuals{p_end}
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

{marker bootopts}{...}
{synoptset 24}{...}
{synopthdr:boot_opts}
{synoptline}
{synopt:{opt r:eps(#)}}number of replications; default 200{p_end}
{synopt:{opt seed(#)}}random-number seed{p_end}
{synopt:{opt str:ata(varname)}}resample within strata{p_end}
{synopt:{opt psu(varname)}}resampling unit; default the household{p_end}
{synopt:{opt svy}}take PSU, strata and weight from {helpb svyset}{p_end}
{synoptline}

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
The elasticities are means over the households of the elasticities of each
household, at its own prices, expenditure and demographics: those of the
{it:households} by default, each household counting for its weight; of the
{it:individuals} with {opt hhsize()}, each household counting for its weight
times its size; of the {it:market}, each household counting for its
expenditure, the elasticities of total demand. A reference household at the
means and the unweighted mean of the household elasticities are available
too. Their standard errors include the sampling of the households as well as
the estimation of the coefficients. All the variances (robust, clustered, survey design) are
built from the same influence functions. See {help equaids##remarks:Remarks}.

{pstd}
Survey data record zero shares for the households that do not buy a good, and
no price for them. {opt pimpute()} fills the missing prices from the other
households of the same group, and {opt selection} corrects the system for the
non-buyers by the two-step method of Shonkwiler and Yen (1999); see
{help equaids##selection:the selection of the buyers}.

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

{dlgtab:Missing prices and non-buyers}

{phang}
{opt pimpute(varlist)} fills the missing prices. Without it, a household with
a missing price leaves the sample; with unit values as prices, that drops the
households that do not buy some good, and a correction for the non-buyers is
then impossible. Each missing log price is replaced by the weighted mean of
the log prices of the households of the sample in the same group of the first
variable (for example the PSU); when nobody in the group has a price, the
group of the next variable is used (for example urban/rural), and so on.
Households still without a price leave the sample, and the means are computed
again until the sample no longer changes: the donors are the households of the
final sample. The weight is that of the estimation (times {opt hhsize()} under
{cmd:elasticities(individuals)}). A note reports, good by good, the prices
filled at each level. The standard errors include the imputation: a household
with a price moves the filled prices of its group, and its influence function
carries that effect (checked against a brute-force influence function). Under
{cmd:vce(cluster)} or {cmd:vce(svy)} at the level of the first grouping
variable the effect cancels within clusters, and only the prices filled at a
wider level contribute. {cmd:vce(bootstrap)} fills the prices again on every
replication. {cmd:estat engel} fills them in the same way on {cmd:e(sample)}.

{phang}
{opt selection} corrects for the households that do not buy (Shonkwiler and
Yen 1999). For every corrected good {it:i}, a probit of purchase (w_i > 0) is
estimated first, with the weights of the estimation, on

{p 12 12 2}
s_i = (1, ln(p_1/x), ..., ln(p_M/x), z, the variables of {opt selvars()} for {it:i}),

{pmore}
and the system is then estimated on all the households, buyers or not, with
the expected shares

{p 12 12 2}
E[w_i] = Phi(s_i'a_i) f_i + delta_i phi(s_i'a_i),

{pmore}
f_i the share of the model (AIDS or QUAIDS) and delta_i a new parameter, the
covariance of the error of the share with that of the probit; delta_i = 0 is
no selection. The last good closes the system, E[w_M] = 1 - sum of the
others: it must be bought by every household (place last a good that everybody
buys, such as the rest of the budget), and it is not corrected. The probit is
written in prices relative to expenditure, so that it is homogeneous of degree
zero, as are the expected shares: the Engel and Cournot aggregation and
homogeneity hold exactly for the market elasticities; symmetry holds for the
shares f of the model, not for the expected shares. The elasticities are those
of the expected demand of all the households, buyers and non-buyers. The
standard errors include the estimation of the probits (their influence
functions are stacked into those of the system). {opt selection} is refused
with {cmd:vce(conventional)}.

{phang}
{opt selgoods(namelist)} restricts the correction to the goods listed (names
of the shares). By default, every good with zero shares in the sample is
corrected, the last good excepted; a good listed that every household buys is
left uncorrected, with a note.

{phang}
{cmd:selvars(}[{it:good}{cmd::}] {it:varlist} [{cmd:;} ...]{cmd:)} gives the
variables of the probits only: they move the purchase of the good but not the
share of those who buy (exclusion restrictions). Segments are separated by
{cmd:;}. A segment without a good name goes to the probit of every corrected
good; {it:good}{cmd::} {it:varlist} to the probit of that good only; the two
add up. With the goods {cmd:wcorn}, {cmd:wwheat} and {cmd:wrice} corrected:

{p2colset 12 48 50 2}{...}
{p2col:{cmd:selvars(dist)}}{cmd:dist} in the three probits{p_end}
{p2col:{cmd:selvars(wrice: perc_ocupa)}}{cmd:perc_ocupa} in the probit of {cmd:wrice} only{p_end}
{p2col:{cmd:selvars(dist ; wrice: perc_ocupa)}}{cmd:dist} in the three, {cmd:perc_ocupa} for {cmd:wrice}{p_end}
{p2col:{cmd:selvars(wcorn: rururb ; wrice: perc_ocupa)}}{cmd:rururb} for {cmd:wcorn}, {cmd:perc_ocupa} for {cmd:wrice}, nothing more for {cmd:wwheat}{p_end}
{p2colreset}{...}

{pmore}
Rules, checked before estimating: a good named must be one of the shares, not
the last one, and corrected (in {opt selgoods()} when it is given); one good
per segment; a variable of {opt selvars()} must not be in the model (a share,
price, expenditure or demographic). Households with a missing value of these
variables leave the sample, with a note. Without {opt selvars()}, the
correction is identified by the nonlinearity of the probit only.

{dlgtab:Variance}

{phang}
{opt vce(robust)}, the default, sums the squares of the influence functions
of the households (factor N/(N-1)). {cmd:vce(cluster} {it:clustvar}{cmd:)}
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
PSUs per stratum, 93.2% to 94.8%.

{phang}
{cmd:vce(bootstrap} [{cmd:,} {it:boot_opts}]{cmd:)} estimates the full sample,
then resamples the households (or the PSUs of {opt psu()}, within the strata of
{opt strata()}; {opt svy} takes both and the weight from {helpb svyset}) with
replacement, and on each replication runs the whole procedure again: the
prices filled by {opt pimpute()}, the probits of {opt selection}, alpha_0 by
its rule, the system. The households keep their weights. The variances of the
coefficients, of delta and of the four types of elasticities are the variances
of the replications; the replications of the parameters are kept for
{cmd:estat engel}, which draws the curve of each one. A replication that
does not converge is dropped and counted ({cmd:e(N_reps_ok)}); the replications
iterate up to 1,000 times unless {opt iterate()} is given. Use it when
{cmd:equaids} warns that the linearized standard errors are unreliable (see
{help equaids##selection:Remarks}).

{phang}
{opt vce(conventional)} gives the conventional variance of the coefficients
(no weights, homoskedastic errors); it cannot be combined with pweights.

{dlgtab:Elasticities and reporting}

{phang}
{opt elasticities(households|individuals|market|reference|hhmean)} sets which
elasticities are reported. The first three are means over the households of
the elasticities of each household, at its own prices, expenditure and
demographics; they differ by the weight of each household in the mean (see
{help equaids##remarks:Remarks}). The names are those of {cmd:easi} and
{cmd:duvm}.

{phang2}
{cmd:households}, the default, gives those of the households: each household
counts for its weight.

{phang2}
{cmd:individuals}, the default with {opt hhsize()}, gives those of the
individuals: each household counts for its weight times its size,
{opt hhsize()}, in the whole estimation (the coefficients, their variance and
the mean). With household data the behaviour is that of the household in all
the types; {cmd:individuals} weights each household by its size: it gives the
elasticities of the household of the average person, not those of a person
within the household. It is an estimation, not a display: after it, only
{cmd:individuals} can be displayed, and it cannot be displayed after another
estimation. The weight used is stored in {cmd:e(wexp)}, so that {cmd:estat}
uses the same.

{phang2}
{cmd:market} gives the elasticities of total demand: each household counts
for its weight times its expenditure, which makes each elasticity a ratio of
totals (the aggregate elasticities).

{phang2}
{cmd:reference} gives those of a reference household, at the weighted means of
ln p, ln x and z, as Poi (2012); these were the {cmd:households} of versions
1.0.0 and 1.1.0.

{phang2}
{cmd:hhmean} averages the household elasticities with the sampling weights
only, as {helpb demandsys} (which does not weight the average). The household
elasticities divide by the predicted shares: {cmd:equaids} warns when some are
near zero or outside [0,1].

{pmore}
All five have standard errors. {cmd:households}, {cmd:market},
{cmd:reference} and {cmd:hhmean} come from the same estimation and can be
displayed on replay ({cmd:equaids, elasticities(market)}). The older names are
accepted: {cmd:aggregate} for {cmd:market}, {cmd:means} for {cmd:reference},
{cmd:household} for {cmd:hhmean}. Version 1.2.0 changed the default (it was
{cmd:market}) and the meaning of {cmd:households} and {cmd:individuals}
(they were evaluated at the means).

{phang}
{opt hhsize(varname)} gives the size of the household, positive. With it and
without {opt elasticities()}, the elasticities are those of the individuals.
It does not enter the model (add it to {opt demographics()} for that); with
another type of elasticities it is not used.

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
{bf:The types of elasticities.} With f_ih the predicted share of good i for
household h and mu_ih = df_ih/d ln x, the expenditure elasticity of type v is
E_i = 1 + sum_h v_h mu_ih / sum_h v_h f_ih, with v_h = w_h ({cmd:households}),
w_h n_h ({cmd:individuals}) or w_h x_h ({cmd:market}); the price
elasticities are built in the same way. Each is the mean of the elasticities
of the households weighted by v_h f_ih: no household is chosen as a
reference, and the mean does not depend on how the means of prices and
expenditure are taken. The market elasticity is the elasticity of the total
demand of the population, what a simulation of a price or tax change on
aggregate demand needs. None divides by a household share, so that goods with
many small or zero shares do not make them unstable; and, as ratios of
weighted totals, their inference under a survey design is standard. The
unweighted mean of the household elasticities ({cmd:hhmean}), by contrast,
divides by each predicted share: on small goods a few households with shares
near zero drive it.

{pstd}
{bf:Standard errors.} The influence function of a households, individuals or
market elasticity has two terms, the sampling of the households (the summary is taken
over a sample) and the estimation of the coefficients; both are analytic. They
agree with the bootstrap of the households and with the Rao-Wu bootstrap of a
survey design; the design variance equals Stata's own linearization. For the
reference household ({cmd:reference}) the sampling term is that of the means
of ln p, ln x and z, and for {cmd:hhmean} that of the
mean of the household elasticities; their derivatives with respect to the
coefficients and to the means are taken by central differences.

{pstd}
The influence function of the coefficients is the exact derivative of the
whole procedure with respect to the weight of each household: the observed
Jacobian of the estimating equations of the coefficients and of Sigma jointly
(the fixed point of iterated FGNLS), and with {opt selection} the probits with
their observed Hessian. A brute-force influence function (the weight of each
household moved, the whole estimation re-run) reproduces the standard errors
to four decimals. The Gauss-Newton information alone, which the robust
variance of {helpb demandsys} uses, leaves out the residuals times the second
derivatives of the shares and the estimation of Sigma: on Poi's food data its
standard errors are up to 8% smaller.

{marker selection}{...}
{pstd}
{bf:The selection of the buyers.} With {opt selection}, a table under the
header reports, for each corrected good, the percentage of buyers, the
pseudo-R2 of its probit (McFadden), the households predicted with probability
0 or 1 ({it:Perfect}: separation), the variance inflation {it:VIF(beta)} of the
expenditure coefficient due to the selection term, 1/(1-rho^2) with rho the
correlation of phi_i and of the column of beta_i in the gradient of the
share, given the other columns, and delta_i with its z. Above 10, phi_i is
almost collinear with ln x: the correction multiplies the standard error of
the expenditure coefficient by more than 3 and rests on the curvature of the
probit; the good can be left uncorrected with {opt selgoods()}, or a
variable of the probit only added in {opt selvars()}. The probits are
estimated good by good (a probit for each good, not a multivariate probit): the
two-step estimator is consistent, not efficient.

{pstd}
{bf:When to trust the analytic standard errors.} The analytic (linearized)
standard errors assume every parameter well identified and the estimate inside
its domain. In a simulation of the model with every parameter strongly
identified (700 samples of 3,000 households, zero shares, missing prices of
the non-buyers filled by {opt pimpute()}, a count demographic, a variable of
the probit only per good), they were 0.95 to 1.02 times the dispersion of the
estimates for delta, rho and the market elasticities. On survey data they can
diverge from the bootstrap when identification is weak: when the probits
barely separate buyers from non-buyers and no variable of the probit only is
available, delta is almost collinear with the rest of the share; and the
correction can move the estimate of Ray's scaling towards its boundary
m0(z) = 1 + rho'z > 0 (a count demographic such as the household size is the
usual case). The estimator is then close to irregular, the linearization
understates the uncertainty, and the bootstrap draws spread widely. The
coefficient rho itself is often weakly identified by Ray's scaling; its
standard error should not be read too literally, while the elasticities are
less affected. {cmd:equaids} reports the distance of the estimate to the
boundary in standard errors ({cmd:e(m0_t)}) and warns below 3; with the
warning, or when in doubt, use {cmd:vce(bootstrap)}.

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
{opt elasticities(reference)}. It is the exact function of the estimate: linear
in ln x for AIDS, quadratic for QUAIDS. The band is the delta method with the
analytic Jacobian and the variance of the estimate, {cmd:e(V_free)}: robust,
by cluster or by design as estimated, with t on the design degrees of freedom
under {cmd:vce(svy)}. Under {cmd:vce(bootstrap)}, the band comes from the
replications themselves, as the standard errors of the elasticities do: the
curve of each replication, with its parameters and its alpha_0
({cmd:e(boot_b_free)}, {cmd:e(boot_anot)}), on the same grid and at the same
means, and the band is the curve plus or minus z times their standard
deviation. The delta method with the bootstrap variance of the parameters
fails where the estimator is not regular: near the boundary of Ray's scaling
(example 7), the parameters of the replications spread far along directions
that hardly move the shares, and the linearized band is far too wide. A
replication whose curve is not defined (m0(z) <= 0 at the means) is left out
and counted. For QUAIDS, the turning point of each curve, where
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
{bf:After selection}, the curves are those of the expected shares,
Phi f + delta phi, at the means of the variables of the probits too; the band
includes the estimation of the probits ({cmd:e(V_sel_psi)}; under
{cmd:vce(bootstrap)}, the replications {cmd:e(boot_sel_psi)}), and no turning
point is marked, the expected share not being quadratic in ln x.

{pstd}
{cmd:estat engel} stores {cmd:r(n)} and, at the means, {cmd:r(turn)} and the
kind of band, {cmd:r(band)} ({cmd:delta} or {cmd:bootstrap}; with
{cmd:bootstrap}, {cmd:r(band_reps)}, the replications used); with
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
{synopt:{cmd:e(m0_t)}}distance of the estimate to the boundary m0(z) > 0, in standard errors{p_end}
{synopt:{cmd:e(N_reps)}, {cmd:e(N_reps_ok)}}bootstrap replications, successful ones ({cmd:vce(bootstrap)}){p_end}
{synopt:{cmd:e(N_clust)}}number of clusters ({cmd:vce(cluster)}){p_end}
{synopt:{cmd:e(N_strata)}, {cmd:e(N_psu)}, {cmd:e(df_r)}}design ({cmd:vce(svy)}){p_end}
{synopt:{cmd:e(chk_engel)}, ...}residuals of the aggregation identities{p_end}
{synopt:{cmd:e(time)}}execution time in seconds{p_end}

{p2col 5 22 26 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:equaids}{p_end}
{synopt:{cmd:e(model)}}{cmd:QUAIDS} or {cmd:AIDS}{p_end}
{synopt:{cmd:e(vce)}}{cmd:robust}, {cmd:cluster}, {cmd:svy}, {cmd:conventional} or {cmd:bootstrap}{p_end}
{synopt:{cmd:e(boot_design)}, {cmd:e(boot_seed)}}resampling unit and seed ({cmd:vce(bootstrap)}){p_end}
{synopt:{cmd:e(anot_rule)}}rule of alpha_0, or {cmd:user}{p_end}
{synopt:{cmd:e(elasticities)}}{cmd:households}, {cmd:individuals}, {cmd:market}, {cmd:reference} or {cmd:hhmean}{p_end}
{synopt:{cmd:e(hhsize)}}the household size variable ({opt hhsize()}){p_end}
{synopt:{cmd:e(wtype)}, {cmd:e(wexp)}}the weight used (times {opt hhsize()} under {cmd:individuals}){p_end}
{synopt:{cmd:e(data_notes)}}notes of the data diagnostics{p_end}
{synopt:{cmd:e(pimpute)}}the grouping variables of {opt pimpute()}{p_end}
{synopt:{cmd:e(selection)}}{cmd:shonkwiler-yen} with {opt selection}{p_end}
{synopt:{cmd:e(selgoods)}}the goods corrected; {cmd:e(selvars)} the {opt selvars()} specification{p_end}
{synopt:{cmd:e(sel_z_}{it:good}{cmd:)}}the variables of the probit of {it:good} only{p_end}
{synopt:{cmd:e(sel_vars)}}all the variables of {opt selvars()}{p_end}

{p2col 5 22 26 2: Matrices}{p_end}
{synopt:{cmd:e(b)}, {cmd:e(V)}}coefficients (Poi's parameterization) and their variance{p_end}
{synopt:{cmd:e(b_free)}, {cmd:e(V_free)}}free parameters and their variance{p_end}
{synopt:{cmd:e(Sigma)}}covariance of the residuals{p_end}
{synopt:{cmd:e(elas_x)}, {cmd:e(elas_u)}, {cmd:e(elas_c)}}market expenditure, uncompensated and compensated elasticities (rows: goods, columns: prices){p_end}
{synopt:{cmd:e(V_elas_x)}, {cmd:e(V_elas_u)}, {cmd:e(V_elas_c)}}their variances{p_end}
{synopt:{cmd:e(elas_xw)}, {cmd:e(V_elas_xw)}, ...}the same for the households
({cmd:w}; the individuals after {cmd:elasticities(individuals)}), for the reference
household ({cmd:m}: {cmd:reference}) and for the household mean ({cmd:h}: {cmd:hhmean}){p_end}
{synopt:{cmd:e(aggshare)}}aggregate budget shares (market){p_end}
{synopt:{cmd:e(shares_m)}, {cmd:e(shares_h)}}predicted shares at the means, and mean of the predicted shares{p_end}
{synopt:{cmd:e(vif)}, {cmd:e(demo_stats)}}data diagnostics{p_end}
{synopt:{cmd:e(sel_delta)}, {cmd:e(se_sel_delta)}}delta and its standard error by good (missing if not corrected){p_end}
{synopt:{cmd:e(sel_alpha)}}probit coefficients by good: constant, ln(p/x), demographics, {opt selvars()} (missing where a variable does not enter){p_end}
{synopt:{cmd:e(sel_diag)}}by good: percentage of buyers, pseudo-R2, perfectly predicted, VIF of beta{p_end}
{synopt:{cmd:e(sel_psi)}, {cmd:e(V_sel_psi)}}free parameters and probit coefficients, and their variance ({cmd:estat engel}){p_end}
{synopt:{cmd:e(boot_b_free)}, {cmd:e(boot_anot)}}the replications of {cmd:e(b_free)} and of alpha_0, one row each ({cmd:vce(bootstrap)}; {cmd:estat engel}){p_end}
{synopt:{cmd:e(boot_sel_psi)}}the replications of {cmd:e(sel_psi)} ({cmd:vce(bootstrap)} with {opt selection}; {cmd:estat engel}){p_end}


{marker examples}{...}
{title:Examples}

{pstd}
The examples use Poi's food data ({cmd:webuse food}), and the Mexican cereals and
simulated non-buyers, ancillary files of the package ({cmd:mexico_2014_cereals.dta},
{cmd:equaids_nonbuyers.dta}): {stata "ssc install equaids, all replace"} (or
{cmd:net get equaids}) copies them into the current folder; the links read them
from there, else from the SSC archive, else from GitHub, and write nothing to
disk. Each one runs
from its blue links: in the command window, in the dialog box (filled in; click
OK), or as a do-file opened in the Do-file Editor. The data in memory are not
lost: the command window and the do-file give them back at the end, even after
an error; the dialog box, which needs the example data in memory, refuses to
replace data that have unsaved changes, unless they are example data loaded
for a dialog box (by equaids, duvm or easi). Files written by the examples go to
Stata's temporary folder. The links call {cmd:equaids_examples} {it:#}
[{cmd:, db} | {cmd:do}].

{title:Example 1: The elasticities of four food groups (Poi's data)}

{phang2}{cmd:. webuse food, clear}{p_end}
{phang2}{cmd:. equaids w1-w4, prices(p1-p4) expenditure(expfd) snames(meat fruitveg bread dairy)}{p_end}
{phang2}{cmd:. equaids, compensated checks stars}{p_end}
{p 8 8 2}{txt}({stata "equaids_examples 1":example 1: click to run in command window}){p_end}
{p 8 8 2}{txt}({stata "equaids_examples 1, db":click to run in dialog box}){p_end}
{p 8 8 2}{txt}({stata "equaids_examples 1, do":open as a do-file}){p_end}

{title:Example 2: The types of elasticities}

{pstd}The elasticities of the households (the default), of the market, of the reference household and the mean of the household elasticities, all from the same estimation.{p_end}
{phang2}{cmd:. webuse food, clear}{p_end}
{phang2}{cmd:. equaids w1-w4, prices(p1-p4) expenditure(expfd) snames(meat fruitveg bread dairy) notable}{p_end}
{phang2}{cmd:. equaids, elasticities(households)}{p_end}
{phang2}{cmd:. equaids, elasticities(market)}{p_end}
{phang2}{cmd:. equaids, elasticities(reference)}{p_end}
{phang2}{cmd:. equaids, elasticities(hhmean)}{p_end}
{p 8 8 2}{txt}({stata "equaids_examples 2":example 2: click to run in command window}){p_end}
{p 8 8 2}{txt}({stata "equaids_examples 2, do":open as a do-file}){p_end}

{title:Example 3: Survey design (Mexican cereals)}

{pstd}The sample keeps the households with all their prices; ten strata then have a single PSU,
which {cmd:singleunit(centered)} handles (with {cmd:singleunit(missing)} their standard errors
would be missing).{p_end}
{phang2}{cmd:. use mexico_2014_cereals, clear}{p_end}
{phang2}{cmd:. svyset psu [pweight=sweight], strata(strata) vce(linearized) singleunit(centered)}{p_end}
{phang2}{cmd:. equaids wcorn wwheat wrice wother wcomp, prices(pcorn pwheat price pother pcomp) expenditure(hh_current_inc) demographics(hhsize isMale) vce(svy)}{p_end}
{p 8 8 2}{txt}({stata "equaids_examples 3":example 3: click to run in command window}){p_end}
{p 8 8 2}{txt}({stata "equaids_examples 3, db":click to run in dialog box}){p_end}
{p 8 8 2}{txt}({stata "equaids_examples 3, do":open as a do-file}){p_end}

{title:Example 4: The elasticities of the individual (Mexican cereals)}

{phang2}{cmd:. use mexico_2014_cereals, clear}{p_end}
{phang2}{cmd:. equaids wcorn wwheat wrice wother wcomp [aw=sweight], prices(pcorn pwheat price pother pcomp) expenditure(hh_current_inc) hhsize(hhsize)}{p_end}
{p 8 8 2}{txt}({stata "equaids_examples 4":example 4: click to run in command window}){p_end}
{p 8 8 2}{txt}({stata "equaids_examples 4, db":click to run in dialog box}){p_end}
{p 8 8 2}{txt}({stata "equaids_examples 4, do":open as a do-file}){p_end}

{title:Example 5: Engel curves}

{pstd}Run from its link, the example writes the curves to Stata's temporary folder.{p_end}
{phang2}{cmd:. webuse food, clear}{p_end}
{phang2}{cmd:. equaids w1-w4, prices(p1-p4) expenditure(expfd) notable}{p_end}
{phang2}{cmd:. estat engel}{p_end}
{phang2}{cmd:. estat engel, lnx level(90) data(curves, replace)}{p_end}
{phang2}{cmd:. estat engel, asobserved observed}{p_end}
{p 8 8 2}{txt}({stata "equaids_examples 5":example 5: click to run in command window}){p_end}
{p 8 8 2}{txt}({stata "equaids_examples 5, do":open as a do-file}){p_end}

{title:Example 6: Diagnose a specification before estimating it}

{phang2}{cmd:. webuse food, clear}{p_end}
{phang2}{cmd:. equaidsdiag w1-w4, prices(p1-p4) expenditure(expfd) sensitivity}{p_end}
{p 8 8 2}{txt}({stata "equaids_examples 6":example 6: click to run in command window}){p_end}
{p 8 8 2}{txt}({stata "equaids_examples 6, db":click to run in dialog box}){p_end}
{p 8 8 2}{txt}({stata "equaids_examples 6, do":open as a do-file}){p_end}

{title:Example 7: The non-buyers (Mexican cereals)}

{pstd}The prices of the non-buyers are filled from their PSU, else from the urban or rural area;
the four cereals with zero shares are corrected, with the share of employed members in the probits
only; the standard errors are those of a bootstrap of the whole procedure, household size in Ray's
scaling bringing the estimate near its boundary under the correction (a few minutes); then the Engel
curves of the expected shares. Example 8 shows a regular case, where the analytic standard errors
are valid.{p_end}
{phang2}{cmd:. use mexico_2014_cereals, clear}{p_end}
{phang2}{cmd:. equaids wcorn wwheat wrice wother wcomp [pw=sweight], prices(pcorn pwheat price pother pcomp)}
{cmd:expenditure(hh_current_inc) demographics(hhsize isMale) pimpute(psu rururb)}
{cmd:selection selvars(perc_ocupa) vce(bootstrap, reps(50) seed(1))}{p_end}
{phang2}{cmd:. estat engel}{p_end}
{p 8 8 2}{txt}({stata "equaids_examples 7":example 7: click to run in command window}){p_end}
{p 8 8 2}{txt}({stata "equaids_examples 7, db":click to run in dialog box}){p_end}
{p 8 8 2}{txt}({stata "equaids_examples 7, do":open as a do-file}){p_end}

{title:Example 8: The non-buyers with analytic standard errors (simulated data)}

{pstd}3,000 households simulated from the model of the correction (the design of the technical
note, Section 6.4): three goods, the first two bought by about 73% and 53% of the households, a
variable of each probit only (q1, q2), no price for the non-buyers, filled from their group of 25
households that share a price shock, and household size in Ray's scaling, far from the boundary
m0(z) > 0. Every parameter is strongly identified: the robust standard errors, clustered by the
group of the imputation, are valid, as a Monte Carlo of 700 such samples confirms (standard errors
within 5% of the standard deviations of the estimates), and are computed in seconds.{p_end}
{phang2}{cmd:. use equaids_nonbuyers, clear}{p_end}
{phang2}{cmd:. equaids w1 w2 w3, prices(p1 p2 p3) expenditure(x) noquadratic demographics(hs) anot(0)}
{cmd:pimpute(grp) selection selvars(w1: q1 ; w2: q2) vce(cluster grp)}{p_end}
{phang2}{cmd:. display "e(m0_t), distance to the boundary m0(z) > 0: " %4.1f e(m0_t) " standard errors"}{p_end}
{p 8 8 2}{txt}({stata "equaids_examples 8":example 8: click to run in command window}){p_end}
{p 8 8 2}{txt}({stata "equaids_examples 8, db":click to run in dialog box}){p_end}
{p 8 8 2}{txt}({stata "equaids_examples 8, do":open as a do-file}){p_end}


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
Ray, R. 1983. Measuring the costs of children.
{it:Journal of Public Economics} 22: 89-102.

{phang}
Shonkwiler, J. S., and S. T. Yen. 1999. Two-step estimation of a censored
system of equations. {it:American Journal of Agricultural Economics} 81:
972-982.


{title:Author}

{pstd}Abdelkrim Araar, Universit{c e'} Laval / PEP, aabd@ecn.ulaval.ca{p_end}
{pstd}Version 1.2.1. Requires Stata 14.2 or later. License: GPL-3.0-or-later.{p_end}

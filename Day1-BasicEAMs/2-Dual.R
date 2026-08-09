rm(list=ls())
library(EMC2)

# This script and the following one analyze data collected from "ajh" in a
# conflict task (the Flanker task) with an added speed vs. accuracy
# emphasis manipulation.

# The task is described in more detail, and instructions provided on how to
# collect your own data in this task in R, in Day1.pptx document. The experiment
# code also has some flexibility to define different types of conflict tasks
# (Stroop & Simon tasks).

# This script introduces methods to work with "dual" diffusion EAMs, and the
# following script covers race models.
#
# Most of the sampling here is relatively quick as we fit only a small data set
# from a single subject, but to work through without having to run the sampling
# for all models load pre-computed samples

# As a homework exercise after the class, apply the following analyses to your
# own data from a conflict task, or if you have already collected your own data,
# analyse that as you follow along!


#### Load and format the data ------------------------------------------------

print( load("data/Andrew_flanker.RData"))

# Format the data to be suitable for analysis by the EMC2 package.
names(dat)[1] <- "subjects"

# Here we have a flanker task, meaning an inherent congruency manipulation (CI)
# Additionally this flanker task included a speed-accuracy emphasis manipulation
# (E)
dat <- dat[,c("subjects","S","CI","E","R","rt")]
head(dat)

# Always good practice to inspect your data for irregularities.
# There were some fast responses due to double button presses, lets remove
dat <- dat[dat$rt>.2,]

# Your data may also contain very slow outlying responses, check for and
#     remove these if you are analyzing your own data.

# A general note:
# Your data must be a data frame with any number of columns but must at least
# have these columns: "R", "rt" and "subjects"
#
# These are reserved for columns in the data frame giving the response factor
# (R), a numeric response  time column (rt), and a factor column specifying the
# subjects.

# You should also avoid factors called "trials" as that name is used internally
# by EMC2.

# Now lets plot correct and error defective densities.
correct_fun <-  function(data) factor(data$S == data$R)
plot_density(dat, factors = c("CI", "E"),
    functions = list(correct = correct_fun), defective_factor = "correct")

# Very high accuracy, and no mistakes in the congruent factor!

# We will fit some standard evidence accumulation models to this data, the
# Diffusion Decision Model (DDM, and its simplified version the Wiener Diffusion
# model, WDM), and three race models, the LBA (Linear Ballistic Accumulator),
# RDM (Racing Diffusion Model) and LNR (log-normal race).

# The expectation is that these models, although appropriate for standard
# decision tasks with speed vs. accuracy manipulations, will have shortcomings
# with conflict data, enabling an illustration of how to critically evaluate
# a cognitive model.

#### WDM ----
### Design ----------------------------------------------------------

# A design and associated model are specified in EMC2 using the "design"
# function. We will now work through how to specify core arguments of the design
# function see:
?design

# a), THE MODEL: The DDM model is specified with this argument:
model <- DDM

# DDM is a function exported by the EMC2 package. The function is for
# internal use so you don't need to understand its structure now (although
# it will be addressed later in the course), but if you are curious you can look
# at it by printing it.
DDM
# The help file of the function provides an overview of the parameters
?DDM()

# b) THE FORMULA. design() specifies the mapping between model parameters and
# the different types of data defined by factors ("cells") using R's linear
# modeling language through a LIST OF FORMULAS, one for each model-parameter
# type, given as a list to the formula argument. We will build up the model
# arguments one-by-one, although typically we do it all in one go.

# 1) First we specify that thresholds (a) are affected by E (i.e., emphasis)
# with parameters:
#   a (the intercept), corresponding to the threshold in the accuracy condition
#   a_Espeed, the "speed minus accuracy" difference, expected to be negative as
#     speed emphasis reduces thresholds.
design_a <- design(formula = list(a ~ E), data= dat, model = DDM)

# This will print out 3 things:
# - A warning saying that all parameters omitted from the formula are assumed
#     constant
# - An overview of the sampled parameters
# - The design matrices, all parameters that do not have a formula mapping
#   simply get a `1`

# In this case, the only design matrix of interest is that of a.
# Note that the `a` parameter holds for all cells, whereas an additional
# `a_Espeed` parameter is added if the E level is `speed`.

# This can be further illustrated by using the mapped_pars function:
mapped_pars(design_a)

# Notice the exp?

# Parameter types in EMC2 have a default transformation used during sampling,
# which ensures they are unbounded:
#
# a) a log transformation for positive parameters e.g., variability parameters
#    like a, t0, s, st0 and sv, and
# b) a probit (pnorm/qnorm) transformation for doubly bounded parameters,
#    like Z, SZ.
#
# Only the mean rate (v) parameter is estimated without a transformation.

# 2) Next we specify the mapping for the drift rates (v).
# Here we are using response coding meaning that the lower boundary of the DDM
# is associated with the first response and the upper boundary with the second
# response. In this case, left is the first response:
levels(dat$R)

# Since we are using response coding we assume that left arrows elicit drift
# rates that are on average towards the lower boundary (negative), and right
# arrows instead drift rates that are towards the upper boundary (positive):

design_av <- design(formula = list(a ~ E, v ~ 0 + S),
                    data= dat, model = DDM)

# Here by specifying ~ 0 + S, we are telling EMC2 to use cell coding.

# Cell coding does not use the standard 'intercept' and 'effect'
# parameterization, instead it just assigns each cell of the design a unique
# value.

mapped_pars(design_av)

# Note that parameter types omitted from the formula are set at default values.
# Both the default values and scales of model parameters are given in the help
# file of the corresponding model ?DDM() in this case.

# Now we want to add parameters that have one value across all trials, in this
# case: non-decision time (t0) and bias (Z).

# EMC2 estimates the bias as a relative bias, meaning that 0.5 means equidistant
# from the lower and upper boundary. 0 is the lower boundary and 1 is the upper
# boundary

# To keep it on the unit interval Z is probit scaled and to keep t0
# positive it logarithmically scaled.

design_WDM_simple <- design(formula = list(a ~ E, v ~ 0 + S, t0 ~ 1, Z ~ 1),
                            data= dat, model = DDM)

# The "WDM" is a special case of the more generally used diffusion decision model
# (DDM) without the latter's between-trial variability parameters (i.e., they
# are set to zero: st0=log(0),sv=log(0),and SZ=qnorm(0)).


# The output of this model will be the same as above, except Z and t0 will be
# added to the sampled parameters output above.

# To illustrate our design a bit better, let's provide some example
# parameter values.
p_vector <- sampled_pars(design_WDM_simple)
p_vector[] <- c(log(.75), -.2, -2, 2, log(.25), qnorm(.5))

# These example values are some initial guesses for what the parameters could
# look like. We can illustrate the implied mapping of these parameters using the
# `mapped_pars` function again, this time also specifying the `p_vector`
# argument.
mapped_pars(design_WDM_simple, p_vector = p_vector)

# This functions shows the implied mapping of our formula to each cell of the
# design.

# Note that boundary separation is larger in the accuracy condition, and drift
# rates are negative for left stimuli and positive for right stimuli.

# The mapped pars also prints the z (rather than Z) parameter, which is the
# absolute starting point, as opposed to relative (Z). It's lower for speed
# compared to accuracy trials, because a is lower on speed trials.

# We can also plot the implied accumulation process of our design
plot(design_WDM_simple, p_vector = p_vector, factors = list(v = "S", a = "E"),
     plot_factor = "S")


# To complete our model we need to include the effect of congruent and
# incongruent stimuli on the drift rate. We expect that the drift rate is higher
# for congruent than for incongruent stimuli.

# Formalizing this expectation might appear tricky, since it means that the
# drift rate for left arrows should be less negative, and for right arrows
# should be less positive.

# One way to achieve this is to use nesting:
design_WDM <- design(formula = list(a ~ E, v ~ 0 + S/CI, t0 ~ 1, Z ~ 1),
                     data= dat, model = DDM)

# This nesting ("CI" nested within "S") estimates separate effect parameters for
# CI for each level of S. This nesting can be illustrated with the mapped_pars
# function again:
mapped_pars(design_WDM)

# Thus for left/right incongruent stimuli, we estimate one additional parameter
# which corresponds to how much the drift rate is affected by incongruent
# information. Let's again show that mapping using a example parameter values.

# This again prints the names of the sampled parameters
p_vector <- sampled_pars(design_WDM)
p_vector[] <- c(log(.75), -.2, -2, 2, .5, -.5, log(.25), qnorm(.5))

# Drift rates are lower for incongruent stimuli
mapped_pars(design_WDM, p_vector = p_vector)

# Lets inspect the final model graphically.
plot(design_WDM, p_vector = p_vector, factors = list(v = c("S", "CI"), a = "E"),
     plot_factor = "S")

# !UPCOMING!
# One final thing to complete the design we will use, we specify that the data
# were truncated at 0.2s. Note that if this was specified in the data file
# (with a column called LT with a value of 0.2 in each row) we would not have to
# specify it here. We can also add in information truncation and censoring
# later, as we describe below.
# design_WDM <- design(formula = list(a ~ E, v ~ 0 + S/CI, t0 ~ 1, Z ~ 1),
                     # data= dat, model = DDM,LT=.2)

design_WDM <- design(formula = list(a ~ E, v ~ 0 + S/CI, t0 ~ 1, Z ~ 1),
                     data= dat, model = DDM)

### Prior ----

# In order to fit the model, we need to make a prior distribution for each
# parameter (where we don't, defaults will be applied). In EMC2, priors are
# always normal (which is one of the reasons why parameters are transformed to
# be unbounded).

# The following prior for intercept parameters is based on the survey of DDM
# fits by Matzke, D., & Wagenmakers, E. J. (2009). Psychonomic Bulletin & Review,
# 16(5), 798–817.
#
# For effects (i.e., differences between conditions) we assume a mean of zero
# (i.e., no direction assumed) and standard deviations with a vague
# commitment to smaller values.

# Note that prior values are specified on the same scale as sampling (i.e.,
# unbounded).

# A prior for a single-subject fit consists of two elements, prior means
# (theta_mu_mean) and a variance-covariance matrix (theta_mu_var). Here we
# assume that the covariances are zero.

# First specify the means
pmean <- c(a=log(.17),a_Espeed=0,v_Sleft=-2.25,v_Sright=2.25,
  'v_Sleft:CIincongruent'=0,'v_Sright:CIincongruent'=0,t0=log(.45),Z=qnorm(.5))

# For the standard deviations we can omit the parameter names,
# but make sure the order follows that of pmean.
psd <- c(.7,.5,2.5,2.5,1,1,.4,.4)

# We can use the "prior" function to combine this information, indicting that
# the prior is at the single subject level.
priorWDM <- prior(design_WDM,pmean=pmean,psd=psd,type="single")

priorWDM

# To see the 0 covariances as well:
summary(priorWDM)

# Note that the standard deviations (psd) we supplied above have been
# translated to variances.

# We can use the "plot_prior" function to sample from and plot the prior.
plot(priorWDM,designWDM,map=FALSE,layout=c(2,4))

# You can also save the prior samples instead of (just) plotting them.
# The prior samples are returned as a coda mcmc.list object
priorSamples <- plot(priorWDM,designWDM,map=FALSE,layout=c(2,4), do_plot = FALSE)
head(priorSamples[[1]])

# NB: In EMC2, alpha just means that the plotted the parameter is a
#     subject-level parameter (and not e.g. a population-level parameter).

# It is more informative to plot the prior on the natural scale and mapped to
# the design cells over which they vary, which is the default (map=TRUE).
plot(priorWDM,designWDM,layout=c(2,4))

# NB: For your paradigm you will sometimes have to tinker with psd to make sure
#     that the main part of the prior encompasses a reasonable data range on the
#     natural scale, using the last command to guide your judgement graphically.

### Fitting ----

# Using make_emc(), we make an object to be used to perform fitting which
# combines the data, design and prior. The sampling type is "single" because we
# are fitting a non-hierarchical model with just one participant.


# By default, three independent chains are sampled, which is useful in order to
# check whether sampling has succeeded, see below.
sWDM <-  make_emc(dat,design_WDM,type="single",prior=priorWDM)

# The likelihood speed-up factor is related to the compression factor.
# RTs are measured with some degree of measurement error caused by hardware
# constraints such as monitor refresh rate and keyboard polling.
# EMC2 tries to respect this measurement error by grouping RTs in small bins
# (by default 1/60 s). Consequently, likelihoods have to be only calculated once
# for each unique combination of  parameters and discretized RT.


# We can now run the fit. This goes through several stages and if all goes well
# gets a default 1000 converged iterations (3000 samples) in the final stage.
#
# 1) "preburn", default 50 iterations with settings that try to quickly search
# for a parameter region with good fit.

# 2) "burn", the sampler runs repeatedly, adding 100 iterations then checking
# the rhat (gelman.diag) for each parameter until the mean is < 1.1 by default.
# Older samples are thrown away if that makes the rhat better.
#
# 3) "adapt", in this stage the sampler estimates an approximate model of the
# posterior in order to make sampling more efficient (less autocorrelated)
# in the next stage, which is evident by the end of this stage, with the chains
# looking more like "fat flat hairy caterpillars".
#
# 4) Finally in the sample stage the efficient sampling model is used to get
# the samples we will use for further analysis, with rhat < 1.1.

# By default the fitting function will report progress in each stage and the
# total run time. The default is also to use parallel processing, allocating one
# core for each chain, so with the default number of chains 3 cores will be
# used. As the chains are independent this (almost) triples the speed of fitting.
sWDM <- fit(sWDM)
# save(sWDM, file = "samples/sWDM.RData")
load("samples/sWDM.RData")

# NB1: You can adjust the number of cores used with the cores_for_chains
#     argument, but will get no speed-up using more cores than there are chains.

# NB2: Once a fit is done you can add more samples by simply running it again.
#      It will automatically start from where it left off (including if it
#      stopped before finishing the sample stage). The number of extra
#      samples added in the sample stage is determined by the "iter" argument
#      (default 1000).

### Convergence ----

# The following function checks whether the final sampling stage has converged,
# using "Rhat" values from coda (see ?gelman.diag). In its basic form Rhat
# test whether the chains are "mixed". EMC2 also splits chains into the first
# and second half before calculating Rhat (so by default it is based on 6
# chains), which makes it sensitive to "stationarity", that is, the distribution
# of a chain being unchanging (flat) over the course of sampling (e.g., a
# constant mean, variance etc.).

# check() also reports the effective sample size (ESS) for each parameter
# (see ?effectiveSize from coda), the number of independent samples obtained
# after taking account of autocorrelation (typically less than the nominal
# number of samples, were default 3 x 1000).


# By default check provides a chain plot of the parameter with the worst Rhat.
# Each chain (set of mcmc samples) is shown in a different colour.
check(sWDM)

# With practice you will be able to recognize when chains are converged, with
# each chain looking like a "flat fat hairy caterpillar" (flat becasue it does
# not move over iterations, fat because most samples are in the middle,  and
# hairy because the distributions have tails).

# To look at chain plots for all parameters we can simply plot the samples.
# This provides a good example of fat,flat, hair caterpillars and good chain
# mixing. On some sampling runs there is an indication of bi-modality in t0.
plot(sWDM, layout = c(2,4))

# We can see why EAMs are "sloppy" models (i.e., their parameters are often
# highly correlated) with a "pairs" plot. Some correlations are due to the
# way the model is parameterised (e.g., v_Sright and v_Sright:CIincongruent),
# others are more structural, e.g., higher thresholds (a) are associated with
# lower non-decision time (t0)).
pairs_posterior(sWDM)

# NB: The other models examined here often have even higher correlations, making
#     them difficult for standard samplers to work with, although sloppiness is
#     not itself an intrinsically bad thing.

# We can also plot the posterior (black) and prior (red) densities together to
# see how much we have learned from the data (i.e., how much the information
# provided by the data has "updated" the prior). Visually updating is indicated
# by the posterior "dominating" (being much higher) and "contracting" (being
# less variable) than the prior.
plot_pars(sWDM,layout=c(2,4))

# A contraction statistic is also provided, with values near 1 indicating strong
# updating of the prior. In this case we see that contraction was least for the
# rate parameters, although it is still quite high. Also, we see our prior
# greatly underestimated the absolute magnitude of the posterior (i.e., rates
# for congruent were much higher, and the effect of conflict on reducing rates
# much greater) than assumed by the prior.

# By default this plot is on the scale of the prior to check whether the priors
# encompass the posterior and hence are unlikely to have been very influential.
# Here we see that this was not the case for v_Sleft and v_Sright, but that the
# data largely "overwhelmed" the prior, so that this is likely to have lead to
# only a small bias.

# We can also use map = TRUE to map the posterior samples back to the cells of
# the experimental design and to the natural scale (the two must occur together),
# which makes interpretation easier. Here we also see that the prior under-
# estimated caution.
plot_pars(sWDM,layout=c(2,4), map = TRUE)

# In this case (with scales that make it hard to see in detail) it can
# also be useful to look at the plot on the scale of the posterior. Setting
# use_prior_lim = FALSE, we see that all estimates are fairly well localized.
plot_pars(sWDM,layout=c(2,4),use_prior_lim=FALSE,map=TRUE)


# Credible intervals (by default 95%) can also be obtained to see how well the
# estimates are localized in tabular form.
credint(sWDM,map=TRUE)

# The generic summary function provides the same information along with
# convergence and efficiency information.
summary(sWDM)

### Posterior testing ----

# The credible function acts much like a Bayesian version of the standard t.test
# function. For example, we might ask whether responding is unbiased, i.e.,
# whether Z = 0.5 (although this example is somewhat redundant with the output
# of summary).

# We see a minor bias to the bottom
credible(sWDM,"Z",map=TRUE,mu=.5)

# The same test can be made on the probit scale, using the default null value of
# mu = 0. As well as providing a credible interval, the output reports a
# probability for a directional test, by default the probability of samples
# less than the null value. Suppose we had reason to believe bias would be
# towards the right response (i.e., positive) then we might be interested in the
# probability of samples greater than zero. These results indicate that
# credibility of that directional hypothesis is low.
credible(sWDM,"Z",alternative="greater")

# Credible intervals cannot prove the null, but that can be done using a
# Savage-Dickey test which approximates a Bayes Factor by the ratio of the
# heights of prior to posterior densities at the null value.

# A large ratio (prior >  posterior) rejects the NULL (i.e., after being
# informed by the data the NULL value is less likely than under the prior).

# A small ratio (prior <  posterior) favors the NULL (i.e., after being informed
# by the data the NULL value is more likely than under the prior).

# Here the latter holds, as is shown graphically by the "hypothesis" function,
# with the Bayes Factor (i.e., the ratio) being returned. Note that there is some
# variability because every time you run hypothesis(), new samples are generated
# from the prior.
layout(1)
hypothesis(sWDM,"Z",selection="alpha")

# Conventionally Bayes Factors of greater than 10 (or less than 0.1) are
# considered "strong" evidence, values between 1/3 and 3 are considered
# "equivocal" and in between as providing "positive" evidence.

# We can also test hypotheses not implicit in the way the model is
# parameterized. For example, suppose we wished to test if the (absolute)
# conflict effects on rates differ for left vs. right arrows.

# To do this we need to provide a function that calculates the necessary
# difference on a parameter vector. This function is then applied to the
# posterior samples, and the credible interval of the difference returned.
fun <- \(x) diff(abs(x[c("v_Sleft:CIincongruent","v_Sright:CIincongruent")]))

# Here the results are equivocal, i.e., we do not have clear evidence either
# for no difference or for a greater effect for right.
credible(sWDM,x_fun = fun)
# Here the inverse is to make the numerical result > 1 (easier to interpret)
1/hypothesis(sWDM,fun = fun,selection="alpha", map =FALSE)

### Goodness of fit ----

# In order to test whether the model fits the data we need to obtain "posterior
# predictive" samples. These are obtained by randomly selecting a posterior
# sample and using it to simulate data in the same design as the empirical data.

# By default 50 data sets are obtained. More samples might sometimes be
# needed to get clear indications about the uncertainty of the model's
# predictions about the data.
ppWDM <- predict(sWDM)

# Global fit is shown by "defective" cumulative distribution function (CDF)
# plots. A CDF shows the cumulative probability RT. The CDF is not defective
# (i.e., it asymptotes on the right at a probability of 1) when only one
# response is made, as is the case here (i.e., perfect accuracy) for some
# conditions

# Observed data are shown in black and posterior predictives in green.
# The green ribbons indicate the 95% uncertainty of the model predictions

# Note that in contrast to frequentist approaches, where uncertainty is in
# terms of hypothetical future replications, in the Bayesian approach "the data
# are the data" (so shown with no uncertainty) and it is the model that is
# uncertain (i.e., its predictions have a distribution).
plot_cdf(dat,ppWDM,layout=c(2,2), factors = c("S", "CI", "E"))

# For the incongruent condition there are two CDFs, one for correct and one for
# error responses, and so data CDFs asymptote at values below 1 (i.e., the
# probability of correct and error responses).

# We make this plot easier to read, first dropping the less interesting stimulus
# factor, and using the correct function we defined above as the defective factor:
plot_cdf(dat,ppWDM,layout=c(2,2), factors = c("E", "CI"),
         functions = list(correct = correct_fun),
         defective_factor = "correct")

# Generally we see that the WDM provides a poor model of the data. Data
# quantiles often barely overlap with predicted quantiles and the asymptotic
# levels are very different in the incongruent speed condition.

### Parameter recovery ----

# It is important to know if model parameters can well estimated from a given
# design. We will check this by simulating data from the present design and then
# fitting and comparing to the generating parameters.

# make_data can take a fitted emc object and simulate based on posterior median
# estimates.
sim <- make_data(sWDM,design_WDM,data=dat)
emc <- fit(make_emc(sim,design_WDM,type="single"),prior=priorWDM)

# The recovery function plots the results and shows that this works for the WDM.
recovery(emc,credint(sWDM)[[1]][,"50%"])

#### DDM ----

# First let's fit the DDM by also estimating the three between-trial variability
# parameters (previously set to zero in the WDM).

# As before, we will assume a conventional model in which E selectively affects
# thresholds and CI selectively affects rates, as with the WDM.

# designDDM <- design(model=DDM,data=dat,LT=.2,
#   formula=list(a~E,v~0+S/CI,Z~1,t0~1,st0~1,sv~1,SZ~1)
# )
designDDM <- design(model=DDM,data=dat,LT=.2,
  formula=list(a~E,v~0+S/CI,Z~1,t0~1,st0~1,sv~1,SZ~1)
)


# The priors for the extra parameters are again based on Matzke & Wagenmakers.
pmean <- c(a=log(0.17),a_Espeed=0,
          v_Sleft=-2.25,v_Sright=2.25,
         'v_Sleft:CIincongruent'=0,'v_Sright:CIincongruent'=0,
         t0=log(.45),Z=0,st0=qnorm(.03),sv=log(1.2),SZ=qnorm(0.38))
psd <- c(.7,.5,2.5,2.5,1,1,.4,.4,0.9,0.7,.75)
priorDDM <- prior(designDDM,pmean=pmean,psd=psd,type="single")
plot(priorDDM,designDDM,layout=c(2,6))

# Fitting is much slower (~50x) than for the WDM because numerical integration
# over the st0 and SZ parameters is required to obtain the likelihood (although
# EMC2 has its own C++ code for this likelihood)

# Hence, dropping st0 and/or SZ speeds fitting greatly. st0 is often not
# required (and we will drop it in later lessons) and SZ only when
# errors are faster relative to correct responses (as occurs when fast
# responding is emphasized, as is the case here so we continue to
# estimate this parameter throughout).
sDDM <-  make_emc(dat,designDDM,type="single",prior=priorDDM)
sDDM <- fit(sDDM)
# save(sDDM, file = "samples/sDDM.RData")
load("samples/sDDM.RData")

# Although slower, sampling works quite well thanks to EMC2's robust and
# efficient sampler.
check(sDDM)

# It is also clear that the data do not strongly constrain the sv and SZ
# parameters.
plot_pars(sDDM,layout=c(2,6))
plot_pars(sDDM,layout=c(2,6),use_prior_lim=FALSE,map=TRUE)

# We will look at goodness of fit after introducing the other EAMs in the next
# lesson. As the basics of posterior parameter testing were covered earlier in
# this lesson, we will largely leave posterior testing of the parameters of the
# models fit in this lesson as an exercise, but make a few comments on parameter
# values to provide some orientation.

# Briefly looking at DDM parameter estimates relative to the WDM we see rates
# are larger and and thresholds greater. sv and SZ being greater than
# zero causes extra errors, so to maintain the same accuracy faster rates and
# higher thresholds are required.
credint(sDDM)









rm(list=ls())
library(EMC2)

# First, as in 1-Dual.R, the data needs to be loaded.
print(load("data/Andrew_flanker.RData"))
# Format the data to be suitable for analysis by the EMC2 package.
names(dat)[1] <- "subjects"
dat <- dat[,c("subjects","S","CI","E","R","rt")]
dat <- dat[dat$rt>.2,]

#### Race models ----

# EMC2 provides three "race" models, the Racing Diffusion Model (RDM), the
# Linear Ballistic Accumulator (LBA) and the Lognormal Race (LNR).

### RDM ----

# First we fit the RDM. All RDM parameters are positive and so estimated on a
# log scale. Negative rates, although conceptually possible, mean that sometimes
# and accumulator may never reach its threshold, and it also invalidates the
# easy analytic solution for the likelihood we use.

# As in the DDM we set moment-to-moment variability to one and we also use the
# simple form with no start-point variability (A) by fixing it to zero (although
# there is a fast analytic solution when A>0 so fitting is still fast).

## Threshold parameterization ----------------------------------------------

# For race models EMC2 automatically makes a "latent response" (lR) factor that
# indexes each accumulator by its corresponding response.
levels(dat$R)

# What it does under the hood, is for every available response a accumulator is
# created. For every trial, these accumulators then race towards their
# thresholds. The first to win determines the response.

# Note that unlike the DDM, race models can accommodate more than two responses
# by having more than two Rlevels, although here we focus on only the binary
# response case.

# We show below how lR can be used to allow for response bias as race models
# do not have a separate response-bias parameter as does the DDM (i.e., Z).

# Matching the WDM/DDM models we assume:

# Thresholds (B) are affected by E and we also allow for response bias that is
# the same for speed and accuracy by an additive combination with lR.

# !UPCOMING!
# designRDM <- design(data=dat,model=RDM,LT=.2,
#                   formula=list(B~E+lR),
#                   matchfun=function(d)d$S==d$lR)

designRDM <- design(data=dat,model=RDM,
                  formula=list(B~E+lR),
                  matchfun=function(d)d$S==d$lR)

mapped_pars(designRDM)

# Here the intercept (B) estimates the threshold in accuracy-left, which acts as
# a baseline condition.

# B_Espeed is the change from accuracy to speed (equally for both left and
# right responses due to the additive linear model). Hence, it quantifies the
# effect of speed emphasis, with negative values indicating a lower threshold
# in speed and positive values a higher threshold for speed. The conventional
# expectation is the former.

# B_lRright is the change from left to right (equally for speed and
# accuracy conditions due to additivity). Hence, it quantifies the bias effect
# (negative values indicating a bias to left and positive a bias to right).

# Given accuracy-left is the baseline note how both B_Espeed and B_lRright are
# added to determine speed-right.


## Drift rate parameterization ---------------------------------------------

# A new argument to design used with race models, "matchfun", is supplied with
# a function that indicates which accumulator (indexed by the lR factor)
# corresponds to the a correct response for different stimuli.

# In simple designs where the levels of a "stimulus" factor (S here) and the
# levels of lR (i.e., the levels of the R factor) are the same the matchfun
# has the following form.
matchfun=function(d)d$S==d$lR

# This automatically makes available a second way of addressing the
# accumulators, a "latent match" factor (lM) with levels that can be used to
# set different rates for matching accumulators (i.e., corresponding to a
# correct response in a given condition) and mismatching accumulators
# (corresponding to incorrect responses).

# Race models have one accumulator for each response, each with its own
# parameters. For accumulation rates it is useful to re-code these parameters in
# terms of the average rate across accumulators (intercept), and d = difference
# between accumulators.

# This re-coding is achieved using a contrast matrix that is supplied to the
# contrasts argument of design, so we create it here for use in the race
# models examined below. The values of 1/2 and -1/2 are chosen so the "d"
# parameter magnitude corresponds to the magnitude of the difference between
# accumulators (AD = average/difference).
ADmat <- cbind(d = c(-1/2,1/2))
ADmat

# We can then use this contrast for the lM factor. In this case, d is positive
# when responding is above chance and equals the match - mismatch rate.

# We will call the d parameter rate "quality" as it indicates how well a
# participant can discriminate between correct and incorrect responses.

# We will call the intercept either "urgency" or more generically rate
# "quantity" as it quantifies the overall level of evidence driving a response.

# # !UPCOMING!
# designRDM <- design(data=dat,model=RDM,LT=.2,
#                   formula=list(v~lM),
#                   matchfun=function(d)d$S==d$lR,
#                   contrasts=list(lM=ADmat))

designRDM <- design(data=dat,model=RDM,
                  formula=list(v~lM),
                  matchfun=function(d)d$S==d$lR,
                  contrasts=list(lM=ADmat))

mapped_pars(designRDM)

# Here, the v parameter is the average rate (over match and mismatching
# accumulators, i.e., urgency or quantity)

# Linear models are applied to an augmented data frame, which is built internally
# by make_data and make_emc using the internal EMC2:::design_model function.
# The augmented data frame is expanded to have a row for each trial for each
# accumulator or latent response (lR).

# For illustration lets look at how it makes lR and lM:
head(EMC2:::design_model(dat,designRDM))

# This function also adds a column "winner" indicating which response
# was made (i.e., the lR row corresponding to R)

## Congruency --------------------------------------------------------------

# Next we'll take into account how congruency affects the drift rate.

# We assume rates vary with CI, both in terms of the difference between
# match and mismatch (analogous to the DDM v parameter) and in terms of urgency
# (which does not have an analogue in the standard DDM).

# Here we use the "*" operator of the linear modeling language, which is
# shorthand where CI*lM = CI + lM + CI:lM.

# v_CIincongruent is the incongruent condition effect on urgency/quantity (i.e.,
# the difference from the congruent condition).

# v_lMd is the average quality for the congruent condition, here we would expect
# a positive value because congruent accuracy is well above chance.

# v_lMd:CIincongruent is the incongruent condition effect on quality
# (i.e., the difference from the congruent condition). Here we would expect a
# negative value because accuracy is lower in the incongruent condition.

# You can see the last contrast is the product of the preceding two, which is
# typical for interaction (":") parameters.

# # !UPCOMING!
# designRDM <- design(data=dat,model=RDM,LT=.2,
#                   formula=list(v~lM*CI),
#                   matchfun=function(d)d$S==d$lR,
#                   contrasts=list(lM=ADmat))

designRDM <- design(data=dat,model=RDM,
                  formula=list(v~lM*CI),
                  matchfun=function(d)d$S==d$lR,
                  contrasts=list(lM=ADmat))

mapped_pars(designRDM)


## The complete design ----------------------------------------------------------------

# Next we combine all the parameter choices we included above into one full
# RDM model. Additionally we estimate one non-decision time (t0) across
# conditions.

# # !UPCOMING!
# designRDM <- design(data=dat,model=RDM,LT=.2,
#                     formula=list(v~lM*CI, B ~ E + lR, t0 ~ 1),
#                     matchfun=function(d)d$S==d$lR,
#                     contrasts=list(lM=ADmat))

designRDM <- design(data=dat,model=RDM,
                    formula=list(v~lM*CI, B ~ E + lR, t0 ~ 1),
                    matchfun=function(d)d$S==d$lR,
                    contrasts=list(lM=ADmat))

## Priors ----

# Priors are based on our informal experience with this model to be
# be reasonable, but still fairly vague, largely just setting a broad scale.
pmean <- c(B=log(2),B_Espeed=0,B_lRright=0,v=log(2),
  v_lMd=1,v_CIincongruent=-0.3,'v_lMd:CIincongruent'=-0.5,t0=log(.2))

# With the aid of plot_prior we choose standard deviations so that the prior
# distributions cover a reasonable range.
psd <-  c(1,.5,.5,1,.5,.5,.5,.5)
priorRDM <- prior(designRDM,pmean=pmean,psd=psd,type="single")
plot(priorRDM,designRDM)

# Note that the priors have long tails, but these are truncated here
# automatically by plot_prior to make distributions easy to visualize.

# Particularly when parameters are not on the natural scale it can be hard to
# anticipate how your choices of prior means map to the design. The following
# function helps you do that.
mapped_pars(designRDM, p_vector = pmean)

# Suppose you wished the prior to reflect findings that thresholds are lower
# in speed than accuracy conditions by setting B_Espeed=log(.75), i.e., 25%
# less than B in the accuracy condition.
pmean['B_Espeed'] <- log(.75)
mapped_pars(designRDM, p_vector = pmean)

# For the accuracy conditions we see as before that the threshold is given by
# exponentiating the intercept (B)
exp(log(2))

# But now for the speed condition adding B_Espeed produces a 25% mlower value
exp(log(2) + log(.75))

# Although it is useful to understand these basic computations, when
# parameterizations become more complicated it is easier and less error prone
# to simply use mapped_pars to guide your choices.

# To finish let's illustrate our mappings with a design plot
plot_design(designRDM, factors = list(v = c("CI", "lM"), B = "E"),
            p_vector = pmean,plot_factor = "E")


## Fitting ----

# Returning to fitting we see that it is fairly fast.
sRDM <-  make_emc(dat,designRDM,type="single",prior=priorRDM)
sRDM <- fit(sRDM)
# save(sRDM,file="samples/sRDM.RData")
# But you can also load it to avoid computation
print(load("samples/sRDM.RData"))

# The results look good.
check(sRDM)
plot(sRDM)

# Priors are dominated and estimates are fairly well localized, although less
# so for the "FALSE" rates due to the low error rate. As is typical, t0
# estimates are less than the DDM.
plot_pars(sRDM,layout=c(2,4))
plot_pars(sRDM,layout=c(2,5),use_prior_lim=FALSE,map=TRUE)

# There are quite strong correlations between the v parameters. As for the DDM
# there is a negative correlation between the B and t0 parameters.
pairs_posterior(sRDM)

### LBA ----

# The LBA is very similar to the RDM, except that
# 1) A is always fit
# 2) between-trial rate variability (sv) replaced with trial rate variability (s)
# 3) Mean rates (v) are unbounded (i.e., can be negative) so are estimated on
#    the natural scale. All other parameters are on a log scale.

# NB: Although the mean of v can be negative, the default setting for the LBA
#     uses a normal rate distribution truncated below at zero. As a consequence
#     rates on any given trial are always positive even when v is negative.
#     This assumption ensures that the accumulators must eventually reach their
#     threshold so a response is guaranteed. Otherwise all rates on a trial
#     might be negative and so no response would occur.

# By convention sv varies with lM, with the intercept (v) fixed for
# identifiability. This typically results in sv for match being less than
# mismatch.

# # !UPCOMING!
# designLBA <- design(data=dat,model=LBA,LT=.2,
#                     matchfun=function(d)d$S==d$lR,
#   formula=list(B~E+lR,v~lM*CI,A~1,sv~lM,t0~1),
#   constants=c(sv=0),contrasts=list(lM=ADmat))
designLBA <- design(data=dat,model=LBA,
                    matchfun=function(d)d$S==d$lR,
  formula=list(B~E+lR,v~lM*CI,A~1,sv~lM,t0~1),
  constants=c(sv=0),contrasts=list(lM=ADmat))

# Again, priors are based on our informal experience with this model to be
# be reasonable but fairly vague.
pmean <- c(B=log(2),B_Espeed=log(1),B_lRright=log(1),v=3,
  v_lMd=2,v_CIincongruent=-.3,'v_lMd:CIincongruent'=-.5,A=log(.5),sv_lMd=log(1),t0=log(.2))
psd <- c(1,.25,.25,3,3,3,3,1,1,.5)

priorLBA <- prior(designLBA,pmean=pmean,psd=psd,type="single")
plot(priorLBA,designLBA,layout=c(3,4))
plot_design(designLBA, factors = list(v = c("CI", "lM"), B = "E"),
            p_vector = pmean, plot_factor = "E")

# Note that when A is estimated the threshold is b = B+A (this is also true for
# the RDM). This is done so that b > A (i.e., you cant start above the
# threshold) is easy to enforce by requiring B and A to be positive.

# The addition is calculated after the mapping, so to see the implied prior we
# must look on the natural scale and ask plot(<prior>) to also do the calculation.
# Both are achieved with the add_recalculated argument. We see the b prior
# is a little broader than the B prior due to the addition of A.
plot(priorLBA,designLBA,layout=c(2,4),add_recalculated = TRUE)

# Fitting is quick and works well
sLBA <-  make_emc(dat,designLBA,type="single",prior=priorLBA)
sLBA <- fit(sLBA)

# save(sLBA,file="samples/sLBA.RData")
# Or you can load
print(load("samples/sLBA.RData"))

# Convergence looks good
check(sLBA)
plot(sLBA)

#  We see that contraction is weaker for rate parameters, particularly quality.
plot_pars(sLBA,layout=c(2,5))

# Looking at the mapped results on the natural scale we see that the mismatch
# rate is very uncertain, reflecting the lack of error responses that do most
# to constrain it.
plot_pars(sLBA,layout=c(2,3),use_prior_lim=FALSE,map=TRUE)


# As commonly occurs when error rates are low, there are strong correlations
# between the different v parameters. As with the LBA B and t0 are negatively
# correlated. As is common with the LBA t0 is estimated as quite small.
pairs_posterior(sLBA)


### LNR ----

# Finally we fit the LNR model, a very simple race model that does not separately
# estimate rates and thresholds, although it can accommodate them conceptually
# and identify them in some cases.

# As in the LBA, between-trial variability is conventionally allowed to vary
# with lM, but in contrast to the other race models nothing needs to be fixed
# for identifiability.

# Here we fit a model where the Lognormal mean (m) is an additive combination
# of lM, CI and E main effects and interactions between lM and CI and lM and E
# (analogous to the d effects in the RDM and LBA). We also add in response
# bias (lR).

# Note that larger values of m result in slower performance, so they can be
# thought of as being directly proportional to thresholds and inversely
# proportional to rates.

# This means that if we interpret "d" parameters as analogous to quality they
# will be negative when accuracy is above chance (as then the value for match
# should be LESS than for mismatch). To make interpretation easier we flip the
# sign of the ADmat (i.e., we use -ADmat) so that

# # !UPCOMING!
# designLNR <- design(data=dat,model=LNR,LT=.2,
#   matchfun=function(d)d$S==d$lR,
#   contrasts=list(lM=-ADmat),
#   formula=list(m~lM*(CI + E)+lR,s~lM,t0~1))
designLNR <- design(data=dat,model=LNR,
  matchfun=function(d)d$S==d$lR,
  contrasts=list(lM=-ADmat),
  formula=list(m~lM*(CI + E)+lR,s~lM,t0~1))

sampled_pars(designLNR)

# The m parameter is the mean over accumulators in the baseline congruent-accuracy
# condition.

# m_CIincongruent is the difference between congruent and incongruent conditions.
# As responding is slower in the incongruent condition, this is expected to be
# negative.

# m_Espeed the difference between speed and accuracy conditions. As responding is
# slower in the accuracy condition, this is expected to be negative.

# m_lMd is analogous to rate quality, where bigger implies better discrimination.
# It estimates quality in the baseline congruent-accuracy condition.

# m_lMd:CIincongruent is the difference in quality between congruent
# and incongruent conditions. As accuracy is less in the incongruent condition
# this parameter is expected to be negative.

# m_lMd:Espeed is the difference in quality between accuracy and speed
# conditions. As accuracy is less in the speed condition this parameter is
# expected to be negative.

# m_lRright is response bias. In total we have the same number of parameters
# as the LBA.

# Priors are based on our informal experience with this model to be
# be reasonable but fairly vague.
pmean <- c(m=-.5,m_lMd=1,m_CIincongruent=0,m_Espeed=0,m_lRright=0,
  'm_lMd:CIincongruent'=0,'m_lMd:Espeed'=0,s=log(1),s_lMd=log(1),t0=log(.2))
psd <- c(2,2,2,2,2,2,2,1,1,1)


priorLNR <- prior(designLNR,pmean=pmean,psd=psd,type="single")
plot(priorLNR,designLNR,layout=c(3,4))

# Fitting is very quick as the LNR likelihood is very fast to compute,
# and convergence is excellent.
sLNR <-  make_emc(dat,designLNR,type="single",rt_resolution=.05,prior=priorLNR)
sLNR <- fit(sLNR)

# save(sLNR,file="samples/sLNR.RData")
load("samples/sLNR.RData")

check(sLNR)
plot(sLNR,layout=c(2,5))

# Contraction is quite good except for the quality (d) parameters, reflecting a
# lack of constraint due to low error rates. t0 is estimated at a similar level
# to the RDM.
plot_pars(sLNR,layout=c(2,5))
plot_pars(sLNR,layout=c(2,5),use_prior_lim=FALSE,map=TRUE)

# Consistent with our observations about v parameters in the other race models,
# there are high correlation among m parameters.
pairs_posterior(sLNR)

### More on posterior testing ----

# One aspect of posterior testing we had not covered previously is directly
# comparing estimated parameters. For example, suppose we wanted to test
# whether the conflict effect for quantity (m_CIincongruent) is bigger than the
# speed effect (m_Espeed).

# We could do this with the "fun" method described previously, but it is
# easier to just supply parameter names
credible(sLNR,c("m_Espeed","m_CIincongruent"))

# In this case the difference is not credible. Note that it is not possible to
# make the test in the same way with the hypothesis function, but the fun method
# still works to test equality.

# This test provides positive (and almost strong) evidence for no difference.
1/hypothesis(sLNR,fun=\(x) diff(x[c("m_Espeed","m_CIincongruent")]))

# As another example, for quality the conflict effect is credibly larger in
# magnitude, with strong evidence for the difference.
credible(sLNR,c("m_lMd:Espeed","m_lMd:CIincongruent"))
hypothesis(sLNR,fun=\(x) diff(x[c("m_lMd:Espeed","m_lMd:CIincongruent")]))

#### Comparing models ----

# EMC2 provides three methods of comparing models, all of which take account
# of both goodness-of-fit and simplicity.

# Two, DIC and BPIC, are relatively insensitive to the prior, with BPIC having
# twice the complexity penalty of DIC and so preferring simpler models.

# The third, marginal deviance (MD) is more dependent on the prior. Bayes
# Factors are MD ratios.

load("samples/sWDM.RData")
load("samples/sDDM.RData")

# All model selection methods favor the LBA, with the LNR second.
compare(list(WDM=sWDM,DDM=sDDM,RDM=sRDM,LNR=sLNR,LBA=sLBA))

# Note that this tables also provides measures of misfit (without any complexity
# penalty), minD and meanD. Use the smallest of the two (they are often very
# similar as here), and smaller indicates a better fit. In this case the LBA
# fits best. Note, however, that this is only just one way of deciding between
# models, other factors can also be taken into account.

# The table also provides "EffectiveN", which counts parameters while taking
# account of correlations among them to come up with something analogous to ESS
# for samples. It identifies the DDM as being most complex.

### Model fit ----

# To examine model fit we calculate post-predictives. We simulate twice the
# usual number to get smoother plots
ppWDM <- predict(sWDM,n_cores=9,n_post=200)
ppDDM <- predict(sDDM,n_cores=9,n_post=200)
ppRDM <- predict(sRDM,n_cores=9,n_post=200)
ppLNR <- predict(sLNR,n_cores=9,n_post=200)
ppLBA <- predict(sLBA,n_cores=9,n_post=200)

# For nice plots
correct_fun <- function(data) factor(data$S == data$R)

# First for the race models
plot_cdf(dat,ppLBA, factors = c("E", "CI"), functions = list(correct = correct_fun),
         defective_factor = "correct",xlim=c(.45,1.1))
plot_cdf(dat,ppRDM, factors = c("E", "CI"), functions = list(correct = correct_fun),
         defective_factor = "correct",xlim=c(.45,1.1))
plot_cdf(dat,ppLNR, factors = c("E", "CI"), functions = list(correct = correct_fun),
         defective_factor = "correct",xlim=c(.45,1.1))


# Now for the DDM/WDM
plot_cdf(dat,ppWDM, factors = c("E", "CI"), functions = list(correct = correct_fun),
         defective_factor = "correct",xlim=c(.45,1.1))
plot_cdf(dat,ppDDM, factors = c("E", "CI"), functions = list(correct = correct_fun),
         defective_factor = "correct",xlim=c(.45,1.1))

# Clear misfit of the speed - incongruent condition, marginally improved by the
# DDM. Let's compare the 'best' model of both the race category and the DDM
# category according to model selection.

plot_cdf(dat,list(LBA = ppLBA, DDM=ppDDM),factors = c("E", "CI"),
  functions = list(correct = correct_fun),defective_factor = "correct",
  posterior_args = list(col = c("red", "blue")))

# Note that all the above also works with plot_density
plot_density(dat,list(LBA = ppLBA, DDM=ppDDM),factors = c("E", "CI"),
  functions = list(correct = correct_fun),defective_factor = "correct",
  posterior_args = list(col = c("red", "blue")))
# LBA clearly better than DDM

# It is often hard to judge the fine-grained details of misfit from CDF plots.
# To address this we can also focus on how well the model accounts for effects
# on particular statistics like the mean and SD of RT and error rate.

# This can be done by supply plot_fit with a function to calculate each
# statistic (the same function you would use on the data) through the stat
# argument.
mrt <- function(d)mean(d$rt)
sdrt <- function(d)sd(d$rt)
err <- function(d)mean(d$R!=d$S)

# The resulting plots specify the data value as a vertical line superimposed
# on the posterior predictive distribution of the statistics.

# Again the LBA clearly fits better than the DDM
par(mfrow=c(1,3))
plot_stat(dat, list(LBA = ppLBA, DDM=ppDDM), stat_fun = mrt, layout = NULL,
          main = "mean RT",posterior_args = list(col = c("red", "blue")))
plot_stat(dat, list(LBA = ppLBA, DDM=ppDDM), stat_fun = sdrt, layout = NULL,
          main = "sd RT",posterior_args = list(col = c("red", "blue")))
plot_stat(dat, list(LBA = ppLBA, DDM=ppDDM), stat_fun = err, layout = NULL,
          main = "accuracy", posterior_args = list(col = c("red", "blue")))

# Just for illustrative purposes this could also be split by factors
E_acc <- plot_stat(dat, list(LBA = ppLBA, DDM=ppDDM), stat_fun = err,
          posterior_args = list(col = c("red", "blue")), factors = "S")

# We can save and print out the results in tabular form
E_acc

#### Exercises (see script 4-BasicExercisesAnswers for answers) ----

# 1) Try dropping the st0 parameter from fitting.
# a) By what factor does it increase the speed of fitting?
# b) Use model selection to test the support for including st0.
#
# 2) The RDM did less well than the LBA, but had two disadvantages,
# i) it did not allow s to vary with lM (analogous to sv~lM in the LBA)
# ii) it did not estimate A
#
# a) Fit a model addressing (i), how long does it take relative to the original?
# b) Fit a model addressing (i) and (ii), how long does that take?
# c) Compare the three models and the LBA, which wins?
# d) Which is more beneficial appears more beneficial, adding s~lM or adding A~1
# as well?

# 3) Both RDM and LBA suffer from problems with identifying mismatch rates due
# to low error rates. Try to address this by setting the mismatch rate to a
# constant value of 1 and freely estimate the match rates as before. For the RDM
# use the sv~lM. Note that this makes the race models more like the DDM, as an
# increase in match v increases quantity and quality together.
#
# a) Which model wins between the original and new RDM?
# b) How is the fit of the new RDM model affected?
# c) How are posterior correlations and updating in the new RDM model compared
#    to the old model?
# d) Which model wins between the original and new LBA?
# e) How is the fit of the new LBA model affected?
# f) How are posterior correlations and updating in the new RDM model compared
#    to the old model?
#
# 4) Compare the mapped parameters for the new models that are analogous.
# a) Non-decision time
# b) Rate variability
# c) Thresholds
# d) Mean rates
# e) How does the picture of performance provided by each model differ?



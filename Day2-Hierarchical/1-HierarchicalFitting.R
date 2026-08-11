rm(list = ls())
library(EMC2)
dat <- forstmann
# Note that this data was truncated at 0.25s and 1.5s

# This script is based on the LBA design of the paper introducing EMC2.
# https://osf.io/preprints/psyarxiv/2e4dq_v4

# We have kept relevant descriptions of that tutorial in here, but will
# not go into detail on that.

plot_cdf(dat, factors = c("E", "S"))


#### Design ----

# As the LBA is a race model it has one accumulator representing each
# possible response, in this case separate accumulators for the left and
# right response. EMC2 uses the levels of the "R" factor found in the data
# to construct a factor which denotes the accumulators "lR" (the latent
# response) with level names taken from Rlevels.
#
# The lR factor allows for response bias, analogous to Z in the DDM. In
# race models, response bias can be modeled by allowing different
# thresholds (B) for different responses. Here we will use B~lR to allow
# for different thresholds for the accumulator corresponding to left and
# right stimuli (e.g., a bias to respond left occurs if the left threshold
# is lower than the right threshold).


mapped_pars(design(formula = list(B ~ lR),
                   data = forstmann,
                   model = LBA))

#
# For race models, you often need to provide a "matchfun" function that
# marks which accumulator matches the correct response. Typically, this
# function takes the latent response (lR) and returns TRUE when it matches
# the stimulus (S):

matchfun=function(data)data$S==data$lR

# matchfun is used to automatically to create a latent match (lM) factor
# with levels "FALSE" (i.e., the stimulus does not match the accumulator)
# and "TRUE" (i.e., the stimulus does match the accumulator). This is
# added internally and can also be used in model formula, typically for
# parameters related to the rate of accumulation.

mapped_pars(design(formula = list(v ~ lM),
                   data = forstmann,
                   model = LBA,
                   matchfun = matchfun))


# When working with lM, it's often useful to design an "average and
# difference" contrast matrix. For binary responses, this has a simple
# form:

ADmat <- cbind(d = c(-1/2,1/2))
ADmat

mapped_pars(design(formula = list(v ~ lM),
                   data = forstmann,
                   model = LBA,
                   matchfun = matchfun,
                   contrasts = list(lM = ADmat)))


# When applied to drift rates (v~lM), the contrast parameter v_lMd
# represents the difference between correct and incorrect response drift
# rates. This "rate quality" parameter is analogous to the DDM drift rate.
# When accuracy exceeds chance, rate quality is positive.
#
# The ADmat intercept (v) represents the average accumulation rate across
# accumulators, measuring both stimulus magnitude and decision urgency
# effects. Unlike the DDM, which represents only evidence differences
# between possible responses, this "rate quantity" parameter captures the
# overall tendency to respond. The DDM does not directly represent this
# tendency.
#
# For simplicity, we omit rate bias, although it could be included by
# allowing v~S*lM.
#
# We also allow the drift rate standard deviation to vary as a function of
# match (sv~lM), as has become standard in more recent applications of
# the LBA.
#
mapped_pars(design(formula = list(v ~ lM, sv ~ lM),
                   data = forstmann,
                   model = LBA,
                   matchfun = matchfun,
                   contrasts = list(lM = ADmat)))

# In the paper we find using the hierarchical DDM that there's no evidence for
# a threshold or drift rate difference between neutral and accuracy
# emphasis. This can kind of be seen in the data as well.
correct_fun <- function(data) factor(data$S == data$R)
plot_cdf(dat, functions = list(correct = correct_fun), defective_factor = "correct",
             factors="E", layout = c(1,3))

# We thus construct our model to estimate one shared E parameter for
# the neutral and accuracy condition. To achieve this, we use the
# functions argument.
# Here the E2 function creates a new factor in the data, with levels speed
# and nonspeed.
E2 <- function(data) factor(data$E!="speed",labels=c("speed","nonspeed"))
forstmann$E2 <- E2(forstmann)
head(forstmann)

# Just as with the DDM, we fix one parameter to make the model
# identifiable; in this case, as is common, we chose to fix drift rate
# standard deviation intercept (sv). We also include the truncation.

design_hLBA <- design(data = dat,model=LBA,matchfun=matchfun,
                        formula=list(v~lM*E2,sv~lM,B~E2+lR,A~1,t0~1),
                        contrasts=list(lM = ADmat),
                        constants=c(sv=log(1)),
                        functions=list(E2=E2),LT=.25,UT=1.5)


#
# Whereas the drift rate intercept in the DDM represents response bias,
# which is unlikely to be affected by emphasis instructions, in the LBA
# the intercept codes for urgency effects, which are more likely to be
# affected by emphasis instructions. We therefore allow the intercept to
# also vary by E, yielding a V_E2nonspeed parameter and a v_LMd:E2nonspeed
# parameter. The latter parameter represents the difference in
# accumulation rate between correct and incorrect responses under
# neutral/accuracy relative to the speed emphasis condition.
#
# Let's examine how our parameters map to the design:

mapped_pars(design_hLBA)


#### Priors ----

# We now create a prior. Note that all parameters except v are estimated
# on a log scale (see ?LBA). For the LBA, t0 is typically smaller than for
# the DDM, so we set it to 0.2. We further set increasing thresholds and
# urgency for higher response caution conditions. As we have no response
# bias expectation, B_lRright is set to zero.

mu_mean <- c(v=1, v_E2nonspeed = -.2, v_lMd=1, "v_lMd:E2nonspeed"=.2,
             sv_lMd=log(1),B=log(1), B_E2nonspeed = log(1.5), B_lRright=0,
             A=log(0.25), t0=log(.2))
mu_sd <- c(v=1, v_lMd=0.5, "v_lMd:E2nonspeed"=0.5,
           sv_lMd=.5,B=0.3, B_E2nonspeed=0.3,B_lRright=0.3, A=0.4, t0=.5)

prior_hLBA <- prior(design_hLBA, type = 'standard',mu_mean=mu_mean,
                      mu_sd=mu_sd)

# We can also plot what our prior proposes for the individual accumulation paths.

plot_design(prior_hLBA, factors = list(v = c("E2", "lM"), B = "E2"),
            plot_factor = "E2")

# By default this plots the mapped prior on the group-level mean
plot(prior_hLBA, layout = c(2,3))
# You can also plot the non-mapped prior on the group-level mean
# Which is just a (multivariate) normal distribution
plot(prior_hLBA, layout = c(2,3), map = F)

# For now we've only set priors on the group-level means. But hierarchical
# models also require a prior on the variance (covariance matrix). The default
# setting in EMC2 is to fit the full variance covariance matrix (type = 'standard').
# The default settings for EMC2 lead to the following prior on the group-level variances
# (sigma2)
plot(prior_hLBA, layout = c(2,3), selection = "sigma2", N = 1e4, map = F)

# And the following prior on the covariances
plot(prior_hLBA, layout = c(2,3), selection = "covariance", N = 1e4, map = F)

# But more interpretable, the correlations:
plot(prior_hLBA, layout = c(2,3), selection = "correlation", N = 1e4, map = F)

# We can change these settings using the A and v hyper priors:
prior_new <- prior(design_hLBA, update = prior_hLBA, v = 3)

# This centers the correlation more around 0
plot(prior_new, layout = c(2,3), selection = "correlation", N = 1e4, map = F)

# Or we can adjust the variances, all variances:
prior_new <- prior(design_hLBA, update = prior_hLBA, A = 1)
plot(prior_new, layout = c(2,3), selection = "sigma2", N = 1e4, map = F)

# Or per parameter
prior_new <- prior(design_hLBA, update = prior_hLBA, A = c(t0 = .1))
plot(prior_new, layout = c(2,3), selection = "sigma2", N = 1e4, map = F)

# But for now we'll just use the default values

#### Fitting ----

emc <- make_emc(dat,design_hLBA, prior=prior_hLBA)
hLBA <- fit(emc, cores_per_chain = 4, fileName="tmp.RData")
hLBA <- save(hLBA, file = "samples/hLBA.RData")
hLBA <- get(load("samples/hLBA.RData"))

# The model shows good convergence and generally good efficiency. Note that
# we now have four outputs, the population mean, variance, and correlations and
# random effects (i.e., individual participant fits). Note that by default you 
# correlations are not printed. All look good.
check(hLBA,selection = c("mu", "sigma2","correlation","alpha"))

# Now we repeat this for the 5 parameter DDM, note instead of setting
# sv as a constant (which is typical for the LBA); for the DDM you typically
# set s as a constant, sv is not sufficient for scaling constraints.
# Furthermore, we don't assume response bias is affected by E2; so we
# set the intercept to 0.
Smat <- cbind(d = c(-1, 1))
design_hDDM <- design(data = dat,model=DDM,LT=.25,UT=1.5,
                      formula=list(v~S*E2,a~E2,Z~1,t0~1, sv ~ 1),
                      functions=list(E2=E2),
                      contrasts = list(S = Smat),
                      constants = c(v_E2nonspeed = 0))


mu_mean <- c(v=0, v_Sd=2, "v_Sd:E2nonspeed"=.2, a=log(1), a_E2nonspeed = log(1.5), Z=0,t0=log(.2), sv = log(.3))
mu_sd <- c(v=1, v_Sd=0.5, "v_Sd:E2nonspeed"=0.5, a = 1, a_E2nonspeed = .5, Z=0.3, t0=.5, sv = .3)

prior_hDDM <- prior(design_hDDM, type = 'standard',mu_mean=mu_mean,
                    mu_sd=mu_sd)



emc <- make_emc(dat,design_hDDM, prior=prior_hDDM)
hDDM <- fit(emc, cores_per_chain = 4, fileName="tmp.RData")
save(hDDM,file="samples/hDDM.RData")

# We will look at the results in the next session.







#
# =============================================================================
# 
#   STOP SIGNAL RESPONSE TIME (SSRT) ESTIMATION USING BEESTS IN EMC2
#   MBCN Summer School Day 2
#
# =============================================================================
# 
# This script demonstrates how to analyse stop signal data using the
# Bayesian Estimation of Ex-Gaussian Stop-signal Reaction Time distributions
# (BEESTS), implemented in the EMC2 package.
#
# Code adapted from Michelle Donzallaz 



# First load the packages - dev version of EMC for stop signal models at the moment
rm(list=ls())
remotes::install_github("ampl-psych/EMC2",ref="ss-on-dev", force=TRUE)
library(EMC2)




#
# NOTES ON MODEL PARAMETERS
#
# =============================================================================
# There are some key model parameters for stop signal models, which were covered
# in Dora's lecture. These are:
#
# - Go process (ex-Gaussian): 
#       - mu
#       - sigma
#       - tau
#   per default truncated at zero (exg_lb)
#
# - Stop process (ex-Gaussian): 
#       - muS
#       - sigmaS
#       - tauS
#   per default truncated at zero (exgS_lb)

# - Go trigger failure probability: gf (probit scale)

# - Stop trigger failure probability: tf (probit scale)

# mu, sigma, tau, muS, sigmaS, tauS all require log transformation;
# gf and tf require probit transformation

# A key point that we will return to later - although the mu parameter looks 
# like a drift rate, it actually represents the mean of the distribution, 
# _not_ the rate of evidence accumulation. 



#
# SIMULATE SINGLE SUBJECT STOP SIGNAL DATA 
#
# =============================================================================
# We will first simulate some stop-signal data with a simple left/right go response

ADmat <- cbind(d = c(-0.5, 0.5))



# MODEL DESIGN
# -------------
# The way you specify a model is quite consistent throughout EMC 

stopsignal_design <- design(model = SSEXG,                 
                            factors = list(subjects = 1, S = c("left", "right")),         
                            Rlevels = c("left", "right"),      
                            matchfun = function(d) as.numeric(d$S) == as.numeric(d$lR),
                            contrasts = list(lM = ADmat),
                            formula = list(
                              mu ~ lM, sigma ~ 1, tau ~ 1,    # go parameters
                              muS ~ 1, sigmaS ~ 1, tauS ~ 1,  # stop parameters
                              gf ~ 1,                         # go trigger failure     
                              tf ~ 1                          # stop trigger failure
                              ))


# NOTES ON DESIGN CALL 
#    - SSEXG is the ex-Gaussian stop signal model. It models both go RTs and SSRTs 
#       as ex-Gaussian distributions
#
#    - "factors" defines the experimental factors:
#         - subjects: single subject (1) for now 
#         - S is the stimulus direction - the direction of the go signal in this case 
#
#   - "Rlevels" is the possible response options (left or right hand responses)
#
#   - "matchfun" is a matching function that works out if the response direction 
#      matches the stimulus direction 
#
#   - "contrasts" are the worst bit
#      ADmat is an averaging/difference contrast matrix applied to lM (and lM is
#      the latent match factor, coding whether S and R are matched)
#
#   - "formula" is where you set up up which factors the parameters can vary by.
#      If you're using ~ 1, the parameter does not vary by any of the factors.
#      We're only setting the go mu to vary as a function of the stimulus-response match, 
#      hence mu ~ lM
#       
#




# Make the parameter vector 
p_vector <- sampled_pars(stopsignal_design, doMap = FALSE)
p_vector

p_vector[1:length(p_vector)] <- c(
  # mu, mu_lMd, sigma, tau
  log(0.6), log(0.8), log(0.06), log(0.3),
  # muS, sigmaS, tauS
  log(0.25), log(0.02), log(0.05),
  # gf, tf
  qnorm(0.15), qnorm(0.2)
)

# You can use mapped_pars to check if the specified design and parameters are 
# sensible 
mapped_pars(p_vector=p_vector, stopsignal_design)


# For race models to produce higher than chance accuracy, the matching runner
# should be faster than the mismatching runner (i.e., the "correct" runner needs
# be faster than the "incorrect" runner). In the mapped_pars call, you can see 
# that mu is larger for matching than for mismatching.
# Note that the value of the runner here is NOT the drift rate - this is a 
# measurement model not a process model.  
# A higher value of mu means the response would be SLOWER. You can see the 
# matching values in mapped_pars (where S = R) are smaller, so they are FASTER. 



# DEFINE THE SSD 
# -------------

RNGkind("L'Ecuyer-CMRG")
set.seed(4293)


# Define the SSD staircase procedure - these are all set to be very standard 
# - staircase = TRUE: use an adaptive staircase rather than fixed SSDs (better)
# - SSD0 = 0.25: starting SSD of 0.25 seconds
# - stairstep = 0.05: SSD increases/decreases by 50ms after each stop trial
#   depending on whether the participant successfully inhibited or not
# - stairmin = 0: SSD cannot go below 0ms
# - stairmax = Inf: no upper limit on SSD
# - p_stop = 0.25: 25% of trials are stop trials
ssd_staircase <- make_ssd(
  staircase = TRUE,
  SSD0 = 0.25,
  stairstep = .05,
  stairmin = 0,
  stairmax = Inf,
  p_stop = 0.25
)



# SIMULATE DATASET
# -------------

dat <- make_data(
  p_vector, 
  stopsignal_design, 
  n_trials = 400,
  functions = list(SSD = ssd_staircase)
)


# Required columns in the data: subjects, S, R, SSD, rt
# - SSD must be numeric (Inf on go-trials) and in seconds
# - rt is response time in seconds; NA for non-responses





#
# FIT A SINGLE SUBJECT'S DATA
#
# =============================================================================
# Fit the simulated data (which has one person)



# PRIOR
# -------------


prior_ss <- prior(
  stopsignal_design, 
  type = "single",
  pmean = c(log(0.6), 0, log(0.5), log(0.5), log(0.7), log(0.5), log(0.5), qnorm(0.5), qnorm(0.5)),
  psd = c(rep(1, 7), c(1, 1))
)

plot(prior_ss,map = FALSE)





# FIT THE MODEL
# -------------


single_sub_object <- make_emc(data = dat, 
                              design = stopsignal_design,
                              compress = F,     # rt binning - dev thing 
                              type = "single",  # one person
                              prior_list = prior_ss)
#single_sub_model <- fit(single_sub_object,fileName = "single_sub_model.RData")

#Takes a bit of time to run - you can just load the model
single_sub_model <- get(load("single_sub_model.RData"))

summary(single_sub_model)
check(single_sub_model)
credint(single_sub_model)


# 
# POSTERIOR PREDICTIVE CHECKS
# -------------

ppsamps <- predict(single_sub_model, n_post = 100)
head(ppsamps)
plot_cdf(dat, ppsamps, factors = "S")




# Back to powerpoint 





#
# REAL DATA
#
# =============================================================================
# Stop signal data can be very challenging in real life 



real_df <- read.csv('real_stopsignal_data.csv')
real_df$S <- factor(real_df$S)
real_df$R <- factor(real_df$R)
real_df$diagnosis <- factor(real_df$diagnosis)



# 
# ASSUMPTIONS
# -------------



# Stop signal models have a lot of assumptions - critical that your data meet 
# these assumptions in order for the model to be valid. There are lots of 
# assumptions, e.g.:
#     - Go trials should be very accurate (hopefully 95-98%)
#     - Go omissions should be low (they should respond on almost all go trials)
#     - Stop trial inhibition rate, p(inhibit|signal), should be ~0.5
#     - Stop respond RTs should be faster than mean go RT 




## Let's first look at the distributions of the failed stops 
## and the go responses
## Failed stop responses should be faster than the go response 


# Simulated dataset
plot_density(dat,
             factors = c("trial_type", "S"),
             functions = list(
               trial_type = \(d) factor(
                 ifelse(is.finite(d$SSD), "stop", "go"),
                 levels = c("go", "stop"))))


# Real dataset
plot_density(real_df,
  factors = c("trial_type", "S"),
  functions = list(
    trial_type = \(d) factor(
      ifelse(is.finite(d$SSD), "stop", "go"),
      levels = c("go", "stop"))))



# 
# VISUALISE SSDs
# -------------
# This is an extremely useful way to look at data in the stop signal task

# Plot staircase behavior across stop trials - simulated data 
layout(1)
stop_dat <- dat[is.finite(dat$SSD), ]
go_dat <- dat[!is.finite(dat$SSD), ]
plot(stop_dat$trials, stop_dat$SSD, ylim = c(0, 0.7),
     xlab = "Stop trials over time", ylab = "SSD")


#####
##### SSDs on the powerpoint 
#####


# 
# INHIBITION FUNCTION
# -------------
# Response probability should increase with SSD - if the stop signal appears 
# later, they should be more likely to response 
plot_ss_if(dat, use_global_quantiles = FALSE, probs = seq(0, 1, length.out = 5))



# 
# SSRT PLOTS
# -------------
# The signal-response RT (go responses on the stop trials) should increase 
# with SSD. As the stop signal appears later, their response times will go up. 
plot_ss_srrt(dat, use_global_quantiles = TRUE, probs = seq(0, 1, length.out = 5))







# 
# INDIVIDUAL DESIGN
# -------------


individual_design <- design(formula = list(mu ~ lM,
                                           sigma ~ 1, 
                                           tau ~ 1,
                                           muS ~ 1, 
                                           sigmaS ~ 1, 
                                           tauS ~ 1,
                                           gf ~ 1, 
                                           tf ~ 1),
                            data = real_df, model = SSEXG,
                            contrasts = list(lM = ADmat),
                            matchfun = function(d) as.numeric(d$S) == as.numeric(d$lR))


# Parameter vector 
p_vector <- c(
  # mu, mu_lMd, sigma, tau
  log(0.6), log(0.8), log(0.06), log(0.3),
  
  # muS, sigmaS, tauS
  log(0.25), log(0.02), log(0.05),
  
  # gf, tf
  qnorm(0.15), qnorm(0.2)
)



# 
# GROUP DESIGN
# -------------
# For now we're not going to use a prior



# This is essentially the maximal model - we could also set gf ~ diagnosis, but 
# this parameter tends to be harder to fit and is not particularly of interest so 
# I have left it out. 
full_group_design <- group_design(formula = list(mu ~ diagnosis, 
                                                 mu_lMd ~ diagnosis,
                                                 sigma ~ diagnosis, 
                                                 tau ~ diagnosis,
                                                 muS ~ diagnosis, 
                                                 sigmaS ~ diagnosis, 
                                                 tauS ~ diagnosis, 
                                                 gf ~ 1, 
                                                 tf ~ diagnosis),
                                  subject_design = individual_design, # individual level model 
                                  data = real_df)



full_object <- make_emc(real_df, 
                        individual_design, 
                        group_design = full_group_design)

#full_model <- fit(full_object, fileName = "full_clinical_model.RData")

full_clinical_model <- get(load("full_clinical_model.RData"))







# 
# CHECKS
# -------------


# Have some convergence issues - let's see which parameter is causing us grief
summary(full_clinical_model)


# SigmaS is causing us the issues - Rhat is 1.101 and ESS is very small 
# Note that these are on the log scale. Let's transform it back to look at sigma 

exp(c(-4.791, -4.286, -3.853))
# So 8–21 ms, which are extremely tight values 

# Let's look at check
check(full_clinical_model)
# Here we can see that it's not mixing well 

# Where to go from here? 
# You could try running it for longer, you could simplify the model
# Can't interpret this parameter at the moment - here you would need to take 
# some time to try to improve the model  





# Posterior summaries (map=TRUE for seconds and probability scale)
summary(full_clinical_model, map = TRUE)

# Note mu_lM true and false - lmTRUE is faster as we would expect













# ------------------------------------------------------------------------------
# Posterior Predictive Checks
# ------------------------------------------------------------------------------

# Simulate 100 datasets from posterior predictive distribution
ppsamps <- predict(single_sub_model, n_post = 100)
head(ppsamps)
plot_cdf(dat, ppsamps, factors = "S")

# ------------------------------------------------------------------------------
# Compute Mean SSRT
# ------------------------------------------------------------------------------

# In practice, one typically wants to compute mean SSRT.
# For an untruncated ex-Gaussian stop distribution this is muS + tauS.
# EMC2's ex-Gaussian stop distribution is lower-truncated at exgS_lb = 0 by
# default, so mean SSRT is computed as E[X | X > exgS_lb] for
# X = Normal(muS, sigmaS) + Exponential(mean = tauS).

# EMC2 computes this under the hood with ssrt_summary(),
ssrt_summary(single_sub_model)

# but we can also do it manually:

ssrt_exg_mean <- function(mu, sigma, tau, lb = 0) {
  n <- max(length(mu), length(sigma), length(tau), length(lb))
  mu <- rep(mu, length.out = n)
  sigma <- rep(sigma, length.out = n)
  tau <- rep(tau, length.out = n)
  lb <- rep(lb, length.out = n)

  mean <- rep(NA_real_, n)

  untruncated <- is.infinite(lb) & lb < 0
  mean[untruncated] <- mu[untruncated] + tau[untruncated]

  finite_lb <- is.finite(lb)
  if (any(finite_lb)) {
    idx <- which(finite_lb)
    b <- (lb[idx] - mu[idx]) / sigma[idx]
    k <- sigma[idx] / tau[idx]
    q <- pnorm(b, lower.tail = FALSE)
    phi <- dnorm(b)
    i0 <- exp(-(lb[idx] - mu[idx]) / tau[idx] + 0.5 * k^2) * pnorm(b - k)
    survival <- q + i0

    mean[idx] <- ((mu[idx] + tau[idx]) * q + sigma[idx] * phi +
                    (lb[idx] + tau[idx]) * i0) / survival
  }

  mean
}



pars <- parameters(single_ss, selection = "alpha", map = TRUE)
pars$ssrt_mean <- ssrt_exg_mean(pars$muS, pars$sigmaS, pars$tauS, exgS_lb=0)
hist(pars$ssrt_mean, xlab = "seconds", main = "Mean SSRT")
abline(v = mean(pars$ssrt_mean), lty = 3, lwd = 5, col = "darkblue")
mean_ssrt <- mean(pars$ssrt_mean)

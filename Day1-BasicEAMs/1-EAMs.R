rm(list=ls()) # clear R environment
# This script allows you to simulate data from five standard Evidence
# Accumulation models, the WDM, DDM, RDM, LBA and LNR so you can gain an
# understanding of the different parameter types by looking at the effect on
# the simulated data of changing their values.

# Getting Started ----

# # If not already done, install the EMC2 package, you could use CRAN
# install.packages("EMC2")
# or from our Github
remotes::install_github("ampl-psych/EMC2",ref="dev")

library(EMC2)

#### Dual diffusion models ----

# The help defines the parameters and the scales on which they are estimated.
# The corresponding object is a function containing a list that defines the
# DDM
?DDM

# The most general version is called the Diffusion Decison Model or Drift
# Diffusion Model, but lets first look at the wiener diffusion model, WDM.

### WDM ----

# We will go into designs in more detail in the next lesson. For now
# we set up a minimal design.
designWDM <- design(model = DDM,
                    factors=list(subjects=1,S=1:2), # S is the stimulus 1, or 2
                    Rlevels = 1:2,                  # Corresponding responses
                    # Single value of each WDM parameter except rates
                    # where there is on value for each stimulus
                    formula = list(a ~ 1, v ~ S, t0 ~ 1, Z ~ 1))

# Note that as our formula does not mention the three between-trial variability
# parameters of the full DDM they are set to their default values, which
# result in between-trial variability being absent.

# We can look at the rate mapping symbolically
mapped_pars(designWDM)

# Choose parameters for simulation, first getting a vector with the right name
parsWDM <- sampled_pars(designWDM)

# Then choosing some values. Note that we include transforms here to illustrate
# the scale on which sampling occurs
parsWDM[] <- c(log(1),-1,2,log(.3),qnorm(.5))

# Here are the values
parsWDM

# Here is how they map to the design numerically
mapped_pars(designWDM,parsWDM)

# and graphically
plot_design(designWDM,p_vector=parsWDM)

# Now lets make a large amount of data
datWDM <- make_data(parsWDM,designWDM,n_trials=1e4)

# Here we see the data format used by EMC2. Trials is added to simulated data
# but need not be in real data. Also, the LT/UT/LC/UC columns indicate
# whether an truncation (T) or Censoring (C) has been applied to the upper (U)
# or lower (L) end of the rt distribution. Again these need not be added to
# real data (when absent the defaults here, which indicate no censoring or
# truncation) apply.

head(datWDM)

# Note the use of factors
lapply(datWDM,levels)

# We can plot the data distributions as defective densities. Here we break the
# data down by stimulus. By default R is used as the "defective factor"
plot_density(datWDM,factors="S")

# However, it is often more informative to use accuracy as the defective factor,
# so make a function that scores the data
crctfun <- \(d) d$S==d$R

# We pass the function in a list and specify its name as a string
plot_density(datWDM,factors="S",functions=list(C=crctfun),defective_factor="C")

# We can also plot defective cumulative densities. The points indicate the
# 10th, 50th, and 90th percentiles
plot_cdf(datWDM,factors="S",functions=list(C=crctfun),defective_factor="C")

# Here is a little helper function to examine the statistics of the
# simulated data, feel free to augment as you think is useful
stats <- function(dat,digits=2)
{
  crct <- crctfun(dat)
  accuracy <- tapply(crct,dat$S,mean)
  mrt <- tapply(dat$rt,dat$S,mean)
  mrtc <- tapply(dat$rt[crct],dat$S[crct],mean)
  mrte <- tapply(dat$rt[!crct],dat$S[!crct],mean)
  qs <- do.call(rbind,tapply(dat$rt,dat$S,quantile,probs=c(.1,.5,.9)))
  qsc <- do.call(rbind,tapply(dat$rt[crct],dat$S[crct],quantile,probs=c(.1,.5,.9)))
  qse <- do.call(rbind,tapply(dat$rt[!crct],dat$S[!crct],quantile,probs=c(.1,.5,.9)))
  out <- list(accuracy=accuracy,all=cbind(mrt,qs),
                 correct=cbind(mrtc,qsc),error=cbind(mrte,qse))
  print(lapply(out,round,digits=digits))
  invisible(out)
}

stats(datWDM)

# Note that correct and error RT distributions are identical (up to sampling
# error).

# EXERCISE: try different parameter values, how does it affect the plots and
#           stats? e.g., how does response bias affect correct and error RT
#           distributions for each stimulus?


#### DDM ----

# To look at the full DDM we need to specify the three types of between-trial
# variability in the formula. Note that because the experimental design
# (i.e., independent and dependent variables) is the same for all of our
# examples we can just pass in the simulated WDM data and the design function
# will pull out the factors and Rlevels information.
designDDM <- design(model=DDM,data=datWDM,
  formula = list(a ~ 1, v ~ S, t0 ~ 1, Z ~ 1, sv ~ 1, SZ ~ 1, st0 ~ 1))

# Pick some parameters for simulation
parsDDM <- sampled_pars(designDDM)

# Note the scales of sv, SZ, and st0
parsDDM[] <- c(parsWDM,log(1),qnorm(.25),log(0.05))
parsDDM

# Look at the parameters
mapped_pars(designDDM,parsDDM)
plot_design(designDDM,p_vector=parsDDM)


datDDM <- make_data(parsDDM,designDDM,n_trials=1e4)
plot_density(datDDM,factors="S",functions=list(C=crctfun),defective_factor="C")
plot_cdf(datDDM,factors="S",functions=list(C=crctfun),defective_factor="C")

# The stats make it clear that the equivalence of correct and error RT
# distributions no longer holds
stats(datDDM)

# EXERCISE: Look at the effects of adding each type of between trial variability
#           separately.We provide an initial example to show how to use plot
#           overlays.

parsDDM1 <- parsDDM; parsDDM1[c("sv","st0")] <- c(log(0),st0=log(0))
datDDM1 <- make_data(parsDDM1,designDDM,n_trials=1e4)
stats(datDDM1)

# Pass both data sets to get a plot overlay
plot_density(list("all"=datDDM,"SZ"=datDDM1),col=c("black","red"),
             factors="S",functions=list(C=crctfun),defective_factor="C")

plot_cdf(list("all"=datDDM,"SZ"=datDDM1),col=c("black","red"),
              factors="S",functions=list(C=crctfun),defective_factor="C")

####  Race models ----

# For race models it is useful to pass a "matchfun", a function to indicate
# which accumulator matches the stimulus.Note that EMC2 automatically creates
# an accumulator factor, lR, with levels equal to Rlevels.
SeqR <- function(d)d$S==d$lR

# Passing the matchfun then creates a new factor that allows the parameters
# (typically rates) to be mapped to the matching vs. mismatching accumulator.
# As the matchfun produces a logical the latent match (lM) factor will have
# levels FALSE and TRUE

### RDM ----
?RDM

# We can now use lM and lR in the formula in the typical way, lR for thresholds
# to allow for response bias and lM for rates to allow for accurate responding.
designRDM <- design(data=datWDM,model=RDM,
                    matchfun=SeqR,
                    formula=list(v~lM, B ~ lR, t0 ~ 1))

# Note that all parameters are on a log scale. Hence, additive effects on the
# sampled scale are multiplicative on the natural scale.
mapped_pars(designRDM)


# Set parameters so match is greater than mismatch with no response bias
parsRDM <- sampled_pars(designRDM)
parsRDM[] <- c(log(1),log(2),log(2),log(1),log(.3))
parsRDM

# Looking at the parameters we can see match twice mismatch and both
# accumulators have the same threshold.
mapped_pars(designRDM,parsRDM)

# The graphical representation superimposes the two accumulators
plot_design(designRDM, factors = list(v = "lM"),p_vector = parsRDM)


# Simulate some data and look at it
datRDM <- make_data(parsRDM,designRDM,n_trials=1e4)
plot_density(datRDM,factors="S",functions=list(C=crctfun),defective_factor="C")
plot_cdf(datRDM,factors="S",functions=list(C=crctfun),defective_factor="C")
stats(datRDM)


# EXERCISE: Note that you can also specify start-point variability (A) and
#           moment to moment variability (s) parameters, try adding them.
#           Errors tend to be slow by default, can you get fast errors?


### LBA ----
?LBA

# Parameters are very similar to the RDM except between trial variability (sv)
# replaces within trial (s), and we always specify A. By default sv=1
designLBA <- design(data=datWDM,model=LBA,
  matchfun=SeqR,formula=list(v~lM, B ~ lR, t0 ~ 1, A~1))

# Note that rates are on the natural scale for the LBA
mapped_pars(designLBA)

# Set parameters so match is greater than mismatch with no response bias
parsLBA <- sampled_pars(designLBA)
parsLBA[] <- c(log(1),log(4),log(.75),log(1),log(.3),log(.25))
parsLBA

# plot_design(designLBA,p_vector=parsLBA,plot_factor="lM")

# We see that B is the gap between the top of the start-point noise (A) and the
# threshold (b).
mapped_pars(designLBA,parsLBA)

# The graphical representation shows the model's ballistic nature.
plot_design(designLBA, factors = list(v = "lM"),p_vector = parsLBA)


# Simulate some data and look at ti
datLBA <- make_data(parsLBA,designLBA,n_trials=1e4)
plot_density(datLBA,factors="S",functions=list(C=crctfun),defective_factor="C")
plot_cdf(datLBA,factors="S",functions=list(C=crctfun),defective_factor="C")
stats(datLBA)

# EXERCISE: Again play with parameter setting and see what they do. How do the
#           the effects of sv in the LBA compare to those in the RMD? Similarly
#           for error speed.

### LNR ----
?LNR

# In the LNR threhold and rates collapse, so we will put both lM and lR
# effects on the lognormal mean (m).
designLNR <- design(data=datWDM,model=LNR,matchfun=SeqR,
  formula=list(m~lM+lR,s~1,t0~1))

# The mean is unbounded
mapped_pars(designLNR)

# Set parameters so match is greater than mismatch with no response bias
parsLNR <- sampled_pars(designLNR)
parsLNR[] <- c(log(1),log(1/3),log(1),log(1),log(.3))
parsLNR

# plot_design(designLNR,p_vector=parsLNR,plot_factor="lM")

# We see that B is the gap between the top of the start-point noise (A) and the
# threshold (b).
mapped_pars(designLNR,parsLNR)

# Simulate some data and look at it
datLNR <- make_data(parsLNR,designLNR,n_trials=1e4)
plot_density(datLNR,factors="S",functions=list(C=crctfun),defective_factor="C")
plot_cdf(datLNR,factors="S",functions=list(C=crctfun),defective_factor="C")
stats(datLNR)

# EXERCISE: Play with parameter settings. Again errors tend to be slow by
#           default. Which can be adjusted by having variability differ between matching and
#           mismatching accumulators. Think about how a change in the threshold
#           vs. rate might affect LNR parameters.


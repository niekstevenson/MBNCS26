rm(list = ls())
load("data/data_ELP_subset.RData")
library(EMC2)

# This is a subset of the pre-processed ELP data, let's for instructive purposes
# first fit the model to the first subject. First inspect the data
head(data, 10)

# Fill in log-freq for non-words
data$log_freq[data$S == "non-word"] <- 0
data$R <- factor(data$R)
data$S <- factor(data$S)

# extract first sub
dat1 <- data[data$subjects == unique(data$subjects)[1],]
Smat <- cbind(d = c(-1, 1))

# Many more errors for words
plot_density(dat1, factors = "S")


# Null model only word-non-word effect
# For compression purposes it matters which columns are included
# It will compress based on all columns in the data. We can't just discard
# unused columns, since some columns might be needed in the construction of 
# functions. 
dat1_null <- dat1[,c("subjects", "S", "R", "rt")]
null_design <- design(list(v ~ S, a ~ 1, t0 ~ 1, Z ~ 1, sv ~ 1),
                      data = dat1_null, model = DDM,
                      contrasts = list(S = Smat))

# Now we construct a reasonable prior.
# We again used a design with drift bias and drift rate
pars <- sampled_pars(null_design)
pars[1:6] <- c(0, 2, log(1), .3, qnorm(.5), log(.3))
prior_null <- prior(null_design, pmean = pars, type = "single")
plot(prior_null, map = T)
plot(prior_null, map = F)

DDM_null <- make_emc(dat1_null, null_design, type = "single", prior_list = prior_null)
DDM_null <- fit(DDM_null)
save(DDM_null, file = "samples/DDM_null_s1.RData")
load("samples/DDM_null_s1.RData")

# Do the standard checks
check(DDM_null)
summary(DDM_null)
# sv poorly updated, for the rest the posterior is nicely peaked
plot_pars(DDM_null)

# Interestingly there seems to be quite a bit of drift bias towards responding 
# non-word, whereas the startpoint bias is towards responding word.
plot_pars(DDM_null, use_prior_lim = F, use_par = c("v", "Z"))

# Remember that non-word is the lower boundary and word is the upper boundary!
# Potentially this person (or people in general) 
# is very quick to see that something is not a word and is therefore faster at 
# accumulating evidence towards the non-word boundary. However, this person
# has an a-priori bias towards responding word. 

# See how well the model fits the data. 
pp_DDM_null <- predict(DDM_null, n_cores = 4)
# Not really a good model of the data though. Slower errors than 
# observed and underestimation of the number of errors. 
plot_cdf(DDM_null, pp_DDM_null, defective_factor = "correct", 
         functions = list(correct = function(d) d$S == d$R))

# We also know that this model fits the data poorly since it's ignoring one 
# important aspect of the data:
plot_cdf(dat1, factors = c("S", "log_freq"))
# The higher the log-frequency the lower the rt! 
# And fewer errors one makes (see plot_cdf results)
cor(dat1$rt, dat1$log_freq)

# Now we'll construct an alternative model, one that also takes
# the log frequency of the word stimuli into account. 
dat1_full <- dat1[,c("subjects", "S", "R", "rt", "log_freq")]
head(dat1_full)


# Now let's construct the alternative model that log-frequency affects drift rate
alt_design <- design(list(v ~ S*log_freq, a ~ 1, t0 ~ 1, Z ~ 1, sv ~ 1),
                     data = dat1_full, model = DDM,
                     contrasts = list(S = Smat),
                     constants = c(v_log_freq = 0))

# Here the covariate is randomly sampled, so having it here as 8.9
# for non-words is erroneous. We set v_log_freq = 0, since it doesn't make
# sense to assume that the bias parameter is affected by the log-frequency. 
# Having log-freq = 0 for the non-words means that log-freq only affects the 
# words

# Update prior
prior_alt <- prior(alt_design, update = prior_null, pmean = c(`v_Sd:log_freq` = 0))
prior_alt

# Far more unique trials, since the covariate for the word stimuli now almost
# per definition is unique
DDM_alt <- make_emc(dat1, alt_design, type = "single", prior_list = prior_alt)
DDM_alt <- fit(DDM_alt)
save(DDM_alt, file = "samples/DDM_alt_s1.RData")
load("samples/DDM_alt_s1.RData")


# Basic checking again
check(DDM_alt)
summary(DDM_alt)

# One interesting result here is that it seems that our model suggests
# that the default accumulation is to say non-word. Potentially, enforcing
# the same drift rate in both directions was a mistake and that quicker
# drift rate for non-words is now captured by the drift bias parameter. 

# Do some posterior predictives
pp_full <- predict(DDM_alt, n_cores = 10)
plot_cdf(DDM_alt, factors = c("S", "log_freq"))

# The model still fits the data poorly, but the tendency that words get
# easier the higher the log-frequency is well captured. 

# compare the two models to see the winning one
compare(list(null = DDM_null, alt = DDM_alt))
# To get to the BayesFactor BF10:
# exp(-.5*(MD:M1 - MD:M2))
exp(-.5*(755-1024))

# NB: If you run this repeatedly you will get slightly different values of 
#     MD for the two models, but the BF will always be huge.

# Similarly we can also get the BayesFactor using the Savage Dickey ratio
# with the `hypothesis` function. 
# Infinite support! That seems good? They should evaluate the same
# but both get super unstable in these far-end tails of the Bayes Factor
hypothesis(DDM_alt, parameter = "v_Sd:log_freq")

# Remember that the savage dickey ratio provides evidence for a nested hypothesis
# The nested hypothesis being that v_Sd:log_freq = 0. 
# Therefore for testing this hypothesis we only needed one model, the alternative
# model, the null_model is a nested version of the alternative model.
# The savage dickey ratio only applies for nested group-level tests
# This kind of comparison I did you can basically never do for a hierarchical 
# model, since the alternative model is that parameters vary on the group-level
# AS WELL AS between subjects.
# When you use the savage-dickey ratio, you are only testing if there's no 
# group-level difference. 

# Group-level inference ---------------------------------------------------
# To test if there's individual differences in an effect 
# of log-frequency on individual drift rates we should again construct
# a null model and alternative model and use `compare` to compare these. 
# However, usually the target of inference is at the group-level, for the
# group-level we can test nested hypotheses using `hypothesis` and the 
# Savage-Dickey ratio. 
# Thus we only construct an alternative model for all data:
data_full <- data[,c("subjects", "S", "R", "rt", "log_freq", "age", "uni")]
data_full$subjects <- factor(data_full$subjects)

# This is still only a subset of 30 participants of the total 
# of 800 participants that participated in the ELP. 
data$subjects <- factor(data$subjects)
table(data$subjects)
# Every subject completed more than 3000 trials!
# The total data set has 2.5 million trials, I'll show the results 
# after this lesson. 

alt_design <- design(list(v ~ S*log_freq, a ~ 1, t0 ~ 1, Z ~ 1, sv ~ 1),
                     data = data_full, model = DDM,
                     contrasts = list(S = Smat),
                     constants = c(v_log_freq = 0))

# This data set also records the age of the participants
# All participants were middle-aged, but there was still considerable spread
# (40-60) so more diverse than you're average uni students
hist(data$age)

# What if we're interested in the hypothesis that the drift rate,
# non decision time and boundary separation are affected by age.
# Here EMC2 deviates a bit from lme4, this is still an area of active development
# so bare with us if there's some bugs here and there. 
# To test that hypothesis, we set up a group-design, we'll additionally test
# the hypothesis that v_Sd:log_freq is affected by age. 
# let's first scale age, so our effect becomes more interpretable:

data_full$age <- data_full$age - mean(data_full$age)
data_full$age <- data_full$age/sd(data_full$age)
sampled_pars(alt_design)
# Here the `` are necessary for interaction effects
group_des <- group_design(list(v_Sd ~ age, `v_Sd:log_freq` ~ age, a ~age,
                               t0 ~age), data = data_full, 
             subject_design = alt_design)

summary(group_des)

# Update our prior
prior_age <- prior(alt_design, type = "standard", group_design = group_des,
                   update = prior_alt)
prior_age
DDM_age <- make_emc(data = data_full, design = alt_design, prior_list = prior_age,
                    group_design = group_des)

DDM_age <- fit(DDM_age, cores_per_chain = 3, fileName = "tmp.RData")
save(DDM_age, file = "samples/DDM_age.RData")
load("samples/DDM_age.RData")

# Unsurprisingly sv is the main culprit, but still looks super healthy
check(DDM_age)
summary(DDM_age)

# Initial assessment of our group-level parameters, here `mu`
# is the implied group-level means, `beta` is the untransformed parameters
credint(DDM_age, selection = "beta")

# Unfortunately on the main branch there's a bug in hypothesis
# but we can already see that the drift/caution goes up with age! 
hypothesis(DDM_age, parameter = "v_Sd_age", selection = "beta")

plot_pars(DDM_age, all_subjects = TRUE)


# Lastly we can also fit nested group-level models, which also include random effects
data$uni <- factor(data$uni)
table(data$uni)

# Here we have some participants from 2 different unis, we can fit a random intercept for each uni
group_des <- group_design(list(v_Sd ~ age, `v_Sd:log_freq` ~ age, a ~age,
                               t0 ~age*uni), data = data_full, 
                          subject_design = alt_design)

# AH previous fails for me, running main and dev, had to stop here.
# Error in if (any(level_counts > 1)) { : 
#   missing value where TRUE/FALSE needed

prior_uni <- prior(alt_design, type = "standard", group_design = group_des,
                   update = prior_age)
prior_uni

DDM_uni <- make_emc(data = data_full, design = alt_design, prior_list = prior_uni,
                    group_design = group_des)

DDM_uni <- fit(DDM_uni, cores_per_chain = 3, fileName = "tmp.RData")
save(DDM_uni, file = "samples/DDM_uni.RData")

credint(DDM_uni, selection = "beta")

# We would probably need the full set to make that type of inference!
# Note that we're working on nested random effects, so that will 
# soon be possible as well. 
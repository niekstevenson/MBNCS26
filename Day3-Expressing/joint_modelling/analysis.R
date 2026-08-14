rm(list = ls())
library(EMC2)
load("data_Eisenberg.RData")

data$task <- factor(data$task)
#
# In this joint modelling example we'll look at 7 decision-making tasks
# That all employ a cognitive conflict manipulation using congruent and
# incongruent stimuli. The data is stored as a list of data sets
plot_density(data, factors = c("task", "conflict_type"), layout = c(1,2))

# Split the data by task
data_list <- split(data, data$task)


# Let's set up our design for the DDM, for the first task we'll be using:
ADmat <- cbind(`_d` = c(-1/2,1/2))

# Since the design we can use across tasks is the same we can just create
# one design and use the same design across all tasks.
design_conflict <- design(data = data_list[[1]], contrasts = list(conflict_type = ADmat),
                          formula =list(v~conflict_type,a~1, t0~1),
                          model = DDM)

pars <- c(2, 1, log(1.5), log(.2))

plot(design_conflict, pars, factors = list(v = "conflict_type"),
     data = data_list[[1]])

# And then to make a joint design, simply replicate it across tasks
designs <- rep(list(design_conflict), 7)
names(designs) <- unique(data$task)
# NB: This also works with different designs for the different data sets

# Now let's set up a prior, here we use the same mean prior across tasks.
prior_blocked <- prior(design = designs,
                       mu_mean = rep(pars, 7))
prior_blocked
summary(prior_blocked)


# We're going to be estimating subsets of correlations, rather than all of them.
# Otherwise we would be estimating 378 correlations, that's more than we have
# subjects!
# To that end we set up parameter groups, all the correlations between
# parameters within a group are to be estimated.
par_names <- names(sampled_pars(designs))
v_pars <- par_names[grepl("|v", par_names, fixed = TRUE)]
v_pars <- v_pars[!grepl("conflict", v_pars, fixed = TRUE)]
d_pars <- par_names[grepl("|v_conflict_type_d", par_names, fixed = TRUE)]
a_pars <- par_names[grepl("|a", par_names, fixed = TRUE)]
t0_pars <- par_names[grepl("|t0", par_names, fixed = TRUE)]
par_groups <- list(v = v_pars, d = d_pars, a = a_pars, t0 = t0_pars)

# Now to complete our model, we give the joint design, data and joint_prior to the `make_emc`
# function and simply run it with `fit`. EMC2 automatically factorizes across joint designs
# making estimation very efficient compared to other samplers.
joint_blocked <- make_emc(data_list, designs,
                          prior_list = prior_blocked,
                          par_groups = par_groups)

joint_blocked <- fit(joint_blocked, cores_per_chain = 15, fileName = "joint_blocked.RData")
save(joint_blocked, file = "joint_blocked.RData")
load("joint_blocked.RData")

# Let's inspect our correlations!
# We make some nice parameter names for plotting
nice_names <- rep(c("v", "d", "a", "t0"),  7)
task_names <- rep(unique(data$task), each = 4)
nice_names <- paste0(task_names, "-", nice_names)

plot_relations(joint_blocked, only_cred = T, nice_names = nice_names)
plot_pars(joint_blocked, selection = "correlation", use_prior_lim = F)

# First let's plot the fit
pps <- predict(joint_blocked, n_cores = 10)

plot_cdf(data_list[[1]], pps[[1]], factors = "conflict_type")
plot_cdf(data_list[[2]], pps[[2]], factors = "conflict_type")

# Pretty poor fit, but the aim here is not model comparison of the subject-level model
# We can do some inference on the correlations, but it's a pretty convoluted
# story. But don't worry we'll get to that later!

# Instead we'll first show the effect of attenuation

# Attenuation -------------------------------------------------------------


# We first fit every subject individually.
# By specifying type is single, we are telling EMC2 not to estimate
# the fully hierarchical variance covariance matrix.
prior_single <- prior(design = designs, type = 'single',
                      pmean = rep(c(2, 1, log(1.5), log(.2)), 7))

joint_single <- make_emc(data_list, designs, prior_list = prior_single, type = "single")
joint_single <- fit(joint_single, cores_per_chain = 15, fileName = "joint_single.RData")
save(joint_single, file = "joint_single.RData")
load("joint_single.RData")
# In a typical application, the parameter estimates of each person would be
# correlated afterwards. This is how we would do that using EMC2.

# This is not typical for analyses (we'd advise you to use the proper)
# hierarchical model). Therefore, it does not have very inherent EMC2
# functionality.
alpha <- do.call(cbind, credint(joint_single, selection = "alpha", probs = .5))
alpha <- as.data.frame(alpha)
alpha$subjects <- rownames(alpha)

# Now we can run a group-level only (hyper) model on the subject-level medians.
two_step <- run_hyper(type = "standard", data = alpha, prior = prior_blocked, par_groups = joint_blocked[[1]]$par_group)
save(two_step, file = "two_step.RData")
load("two_step.RData")

#
# First let's take a look at the attenuation factor, by looking at the variances.
# For this we'll use the plot_pars function in a slightly different way.
# In green the hierarchical estimates, in black the two-step estimates
plot_pars(two_step, true_pars = joint_blocked, selection = "sigma2",
          use_prior_lim = FALSE, layout = c(2,2))

# We can see the consequences of this in the correlation
plot_pars(two_step, true_pars = joint_blocked, selection = "correlation",
          use_prior_lim = FALSE, use_par = "cue|a")



# Latent variable joint models --------------------------------------------


# We use blocked covariance matrices to help reduce the dimensionality
# of what we're estimating. Another, more sophisticated, way to reduce the
# dimensionality is using factor analysis. EMC2 allows you to run
# a factor analysis on top of the cognitive model in a fully hierarchical set up.
# This way you are left with unattenuated factor loadings.
#
# To set up factor analysis we construct it as a subclass of the SEM
# type model:
fa_settings <- make_sem_structure(
  design = designs,
  lambda_specs = list(
    v = v_pars,
    d = d_pars,
    a = a_pars,
    t0 = t0_pars)
)


prior_factor <- prior(design = designs, type = 'SEM',
                      mu_mean = rep(c(2, 1, log(1.5), log(.2)), 7),
                      sem_settings = fa_settings)

joint_factor <- make_emc(data_list, designs, type = "SEM", prior_list = prior_factor,
                         sem_settings = fa_settings)


make_SEM_diagram(joint_factor, par_names = nice_names,
                 width = 700, height = 700, layout = "nicely")

joint_factor <- fit(joint_factor, cores_per_chain = 15, fileName = "joint_factor.RData")
save(joint_factor, file = "joint_factor.RData")
load("joint_factor.RData")

make_SEM_diagram(joint_factor, par_names = nice_names,
                 width = 900, height = 900, layout = "nicely")

make_SEM_diagram(joint_factor, par_names = nice_names, cred_only = TRUE,
                 width = 900, height = 900, layout = "nicely")


# Some basic inference
plot(joint_factor, selection = "std_loadings")
credint(joint_factor, selection = "std_loadings")

par(mfrow = c(1,1))
plot_relations(joint_factor, selection = 'std_loadings', nice_names = nice_names)



# Full SEM analysis -------------------------------------------------------
# Sometimes you're also interested in more complex questions
# How are these factors related to each other?
# How are they related to IQ? or some personality scales?
# That's where the larger framework of SEMs is useful.
data$Age <- scale(data$Age)
sem_settings <- make_sem_structure(
  data = data,
  design = designs,
  covariate_cols = "Age",
  lambda_specs = list(
    v = v_pars,
    d = d_pars,
    a = a_pars,
    t0 = t0_pars),
  g_specs = list(
    v = "Age",
    d = "Age",
    a = "Age",
    t0 = "Age"
  )
)

prior_SEM <- prior(design = designs, type = 'SEM',
                   mu_mean = rep(c(2, 1, log(1.5), log(.2)), 7),
                   sem_settings = sem_settings)

joint_SEM <- make_emc(data_list, designs, type = "SEM", prior_list = prior_SEM,
                      sem_settings = sem_settings)

make_SEM_diagram(joint_SEM, par_names = nice_names)

joint_SEM <- fit(joint_SEM, cores_per_chain = 15, fileName = "joint_sem.RData")
# Ran in an hour. 190 subjects, 180.000 trials, more than 5500 parameters.

save(joint_SEM, file = "joint_sem.RData")
load("joint_sem.RData")

credint(joint_SEM, selection = "std_loadings")
plot(joint_SEM, selection = "factor_regressors")
credint(joint_SEM, selection = "factor_regressors")


# We can also make some nice diagrams of this
make_SEM_diagram(joint_SEM, par_names = nice_names, cred_only = TRUE,
                 width = 1500)



# Let's finish with Bayesian inference
credint(joint_SEM, selection = "factor_regressors")

# Evidence for v, d, and t0 with age.
hypothesis(joint_SEM, selection = "factor_regressors", parameter = "v.Age")
hypothesis(joint_SEM, selection = "factor_regressors", parameter = "d.Age")
hypothesis(joint_SEM, selection = "factor_regressors", parameter = "a.Age")
hypothesis(joint_SEM, selection = "factor_regressors", parameter = "t0.Age")






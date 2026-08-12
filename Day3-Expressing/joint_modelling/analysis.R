rm(list = ls())
library(EMC2)
load("data_Eisenberg.RData")

data$task <- factor(data$task)
#
# In this joint modelling example we'll look at 7 decision-making tasks
# That all employ a cognitive conflict manipulation using congruent and
# incongruent stimuli. The data is stored as a list of data sets
plot_density(data, factors = c("task", "conflict_type"), layout = c(2,2))

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
          use_prior_lim = FALSE)

# We can see the consequences of this in the correlation
plot_pars(two_step, true_pars = joint_blocked, selection = "correlation",
          use_prior_lim = FALSE, use_par = "cue|a")

rm(list = ls())
library(EMC2)
dat <- forstmann
hDDM <- get(load("samples/hDDM.RData"))
hLBA <- get(load("samples/hLBA.RData"))

# First we do model comparison between the DDM and the LBA
compare(list(LBA = hLBA, DDM = hDDM))

#### Model Fit ----

# We generate posterior predictions for the models to see how well it
# fits the data:

pp_LBA <- predict(hLBA, n_cores = 9)
pp_DDM <- predict(hDDM, n_cores = 9)

acc_fun <- function(data) factor(data$S == data$R)
plot_density(dat, post_predict = list(LBA = pp_LBA, DDM = pp_DDM),
         functions = list(correct = acc_fun),
         factors = "E",
         defective_factor = "correct",
         layout = c(1,3))


# As with single-subject, we can also compare arbitrary descriptives
# on the real data to the predictives of both models.
# For example differences in response time and
# error rates between the emphasis conditions:

drt <- function(data){
  all <- tapply(data$rt,data$E,mean)
  out <- c(all['neutral'] - all['speed'],all['accuracy'] - all['speed'])
  names(out) <- c("NTR-SPD", "ACC-SPD")
  return(out)
}

derr <- function(data){
  data$correct <- data$S == data$R
  all <- tapply(data$correct,data$E,mean)*100
  out <- c(all['neutral'] - all['speed'],all['accuracy'] - all['speed'])
  names(out) <- c("NTR-SPD", "ACC-SPD")
  return(out)
}

par(mfrow = c(1,2))
plot_stat(hLBA, list(LBA = pp_LBA), stat_fun = drt,
  xlab = "RT (s) difference", layout = NULL,legendpos = c("topleft","topright"))
plot_stat(hLBA, list(LBA = pp_LBA), stat_fun = derr,
  xlab = "Accuracy (%) difference", layout = NULL,legendpos = c("topleft","topright"))

# So accuracy is quite similar between Accuracy and Neutral,
# but response times are a little different, with people being slower in
# the accuracy condition. This cannot be captured by the model, since it
# forces these to be the same for these effects across all parameters.

# On our final model we can now also start to perform psychological inference
# So far we have mostly looked at 'which model is better'. But from a psychological
# point of view it might be more interesting to ask questions about the experimental
# manipulations.
# So as an example, which parameters are on-average affected by the manipulation.
sampled_pars(hLBA)
# Using hypothesis is different to compare here. Hypothesis uses a trick to compare
# the model to a null model for which the parameter is set to a constant value (H0).
# But this only works for strictly nested models. So here we compare the group-level
# mean of the LBA model, to alternative models for which there is no group-level difference from
# e.g. 0.

# In compare, we check if this model is preferred across subjects, not on average.
# It's like checking if there's a difference in alpha, or in mu.

hypothesis(hLBA, parameter = "v_E2nonspeed")
hypothesis(hLBA, parameter = "v_lMd:E2nonspeed")
hypothesis(hLBA, parameter = "B_E2nonspeed")

# To check the posteriors of inferential targets, we can use credint
credint(hLBA, map = T)
credint(hLBA, map = list(v = "E", B = "E"))








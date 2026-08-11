rm(list=ls())
library(EMC2)

# This script provides followup exercises for the 1-StandardHierarchical.R
# 1. Test whether a model that omits response bias altogether is preferred 

# Answers 1:
matchfun=function(data)data$S==data$lR
dat <- forstmann

ADmat <- cbind(d = c(-1/2,1/2))

E2 <- function(data) factor(data$E!="speed",labels=c("speed","nonspeed"))
# Now we use a custom contrast matrix
E_incr <- contr.increasing(3)
colnames(E_incr) <- c("nonspeed", "acc")

# Now with lM, we set v_Eacc:lMd to be a constant
design_LBABvE2_nolR <- design(data = dat,model=LBA,matchfun=matchfun,
                         formula=list(v~E*lM,sv~lM,B~E2,A~1,t0~1),
                         contrasts=list(lM = ADmat, E = E_incr),
                         constants=c(sv=log(1), `v_Eacc:lMd` = 0),
                         functions=list(E2=E2))

load("samples/LBABvE2.RData")
prior_LBABvE2_nolR <- prior(design_LBABvE2_nolR, update = get_prior(LBABvE2))

LBABvE2_nolR <- make_emc(dat, design = design_LBABvE2_nolR, prior_list = prior_LBABvE2_nolR)
LBABvE2_nolR <- fit(LBABvE2_nolR, cores_per_chain = 4)
save(LBABvE2_nolR, file = "samples/LBABvE2_nolR.RData")

# Now compare the two models
compare(list(lR = LBABvE2, nolR = LBABvE2_nolR), cores_per_prop = 4)

# The model with lR is supported. We generally find that it is important to 
# always include the possibility of response bias, even when the population mean
# is close to unbiased.

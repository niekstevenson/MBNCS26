library(EMC2)
load("complicatedLBA.RData")
fit(emc, cores_per_chain = 3, fileName="complicatedLBA.RData")
rm(list=ls())
library(EMC2)
print(load("data/Andrew_flanker.RData"))
names(dat)[1] <- "subjects"
dat <- dat[,c("subjects","S","CI","E","R","rt")]
dat <- dat[dat$rt>.2,]
ADmat <- cbind(d = c(-1/2,1/2))

# This script provides followup exercises for 3-Race.R.
# 1) Try dropping the st0 parameter from fitting.
# a) By what factor does it increase the speed of fitting?
# b) Use model selection to test the support for including st0.

# 1a) Fitting is took ~27sec secs vs.  5.8mins, a speedup of ~13x

designDDM1 <- design(model=DDM,data=dat,LT=.2,
                     formula=list(a~E,v~0+S/CI,Z~1,t0~1,sv~1,SZ~1)
)
pmean1 = c(a=log(0.17),a_Espeed=0,
           v_Sleft=-2.25,v_Sright=2.25,
           'v_Sleft:CIincongruent'=0,'v_Sright:CIincongruent'=0,
           t0=log(.45),Z=0,sv=log(1.2),SZ=qnorm(0.38))
psd1 = c(.7,.5,2.5,2.5,1,1,.4,.4,0.7,.75)
priorDDM1 <- prior(designDDM1,pmean=pmean1,psd=psd1,type="single")
sDDM1 <-  make_emc(dat,designDDM1,type="single",rt_resolution=.05,prior=priorDDM1)
sDDM1 <- fit(sDDM1)
save(sDDM1,file="samples/sDDM1.RData")

# 1b) The full DDM wins
load("samples/sDDM.RData")
compare(list(FullDDM=sDDM,DDMsvSZ=sDDM1))

# 2) The RDM did less well than the LBA, but had two disadvantages,
# i) it did not allow s to vary with lM (analogous to sv~lM in the LBA)
# ii) it did not estimate A
#
# a) Fit a model addressing (i), how long does it take relative to the original?
# b) Fit a model addressing (i) and (ii), how long does that take?
# c) Compare the three models and the LBA, which wins?
# d) Which is more beneficial appears more beneficial, adding s~lM or adding A~1
# as well?

# i) Time difference of 33.99654 secs
designRDM1 <- design(data=dat,model=RDM,LT=.2,
                     formula=list(B~E+lR,v~lM*CI,t0~1,s~lM),
                     matchfun=function(d)d$S==d$lR,
                     contrasts=list(lM=ADmat),constants=c(s=log(1)))
# Choose a prior assuming no difference in s between match and mismatch
pmean1 <- c(B=log(2),B_Espeed=log(1),B_lRright=log(1),v=log(2),
            v_lMd=log(2),v_CIincongruent=log(1),'v_lMd:CIincongruent'=log(1),t0=log(.2),
            s_lMd=log(1))

# With the aid of plot_prior we choose standard deviations so that the prior
# distributions cover a reasonable range.
psd1 <-  c(1,.5,.5,1,.5,.5,.5,.5,.5)
priorRDM1 <- prior(designRDM1,pmean=pmean1,psd=psd1,type="single")
plot(priorRDM1,designRDM1,layout=c(2,6))
sRDM1 <-  make_emc(dat,designRDM1,type="single",rt_resolution=.05,prior=priorRDM1)

sRDM1 <- fit(sRDM1)
# save(sRDM1,file="samples/sRDM1.RData")

designRDM2 <- design(data=dat,model=RDM,LT=.2,
                     formula=list(B~E+lR,v~lM*CI,t0~1,s~lM,A~1),
                     matchfun=function(d)d$S==d$lR,
                     contrasts=list(lM=ADmat),constants=c(s=log(1)))
# Choose a prior assuming no difference in s between match and mismatch
pmean2 <- c(B=log(2),B_Espeed=log(1),B_lRright=log(1),v=log(2),
            v_lMd=log(2),v_CIincongruent=log(1),'v_lMd:CIincongruent'=log(1),t0=log(.2),
            s_lMd=log(1),A=log(.5))

# With the aid of plot_prior we choose standard deviations so that the prior
# distributions cover a reasonable range.
psd2 <-  c(1,.5,.5,1,.5,.5,.5,.5,.5,.5)
priorRDM2 <- prior(designRDM2,pmean=pmean2,psd=psd2,type="single")
plot(priorRDM2,designRDM2,layout=c(2,6))
sRDM2 <-  make_emc(dat,designRDM2,type="single",rt_resolution=.05,prior=priorRDM2)

sRDM2 <- fit(sRDM2)
# save(sRDM2,file="samples/sRDM2.RData")

# 2a) ~1/3 slow down: 10.4/7.8

# 2b) ~1/3 slow down: 10.1/7.8
# 2c/d) The LBA still wins. Adding s~lM garners most improvement.
load("samples/sRDM.RData")
load("samples/sLBA.RData")
compare(list(LBA=sLBA,RDM=sRDM,RDM1=sRDM1,RDM2=sRDM2))


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

# 3) Note that we simply removed the contrast, so sv reverts to treatment
#    coding, equally we could have used contrasts=list(sv=ADmat)

designRDM0 <- design(data=dat,model=RDM,LT=.2,
                     formula=list(B~E+lR,v~lM*CI,t0~1,s~lM),
                     matchfun=function(d)d$S==d$lR,
                     constants=c(s=log(1),v=log(1),v_CIincongruent=log(1)))
# Choose a prior assming no difference in s between match and mismatch
pmean0 <- c(B=log(2),B_Espeed=log(1),B_lRright=log(1),
            v_lMd=log(2),'v_lMd:CIincongruent'=log(1),t0=log(.2),
            s_lMd=log(1))

# With the aid of plot_prior we choose standard deviations so that the prior
# distributions cover a reasonable range.
psd0 <-  c(1,.5,.5,.5,.5,.5,.5)
priorRDM0 <- prior(designRDM0,pmean=pmean0,psd=psd0,type="single")
plot(priorRDM0,designRDM0,layout=c(2,4))
sRDM0 <-  make_emc(dat,designRDM0,type="single",rt_resolution=.05,prior=priorRDM0)
sRDM0 <- fit(sRDM0)
# save(sRDM0,file="samples/sRDM0.RData")

ppRDM0 <- predict(sRDM0)

# a) The old model wins
compare(list(RDM=sRDM,RDM0=sRDM0))

#  b) Both minD and meanD indicate slighlty better fit for the old model.
#     Graphically there is little to choose between them.
par(mfrow=c(2,8))
ppRDM <- predict(sRDM)
plot_density(dat,ppRDM0,factors=c("S","CI"))
plot_density(dat,ppRDM,factors=c("S","CI"))

# c) However, for the new model posterior correlations are less extreme, and
# updating is improved.
pairs_posterior(sRDM0)
pairs_posterior(sRDM1)

plot_pars(sRDM0,layout=c(2,5))
plot_pars(sRDM1,layout=c(2,5))


designLBA0 <- design(data=dat,model=LBA,LT=.2,
                     formula=list(B~E+lR,v~lM*CI,t0~1,sv~lM,A~1),
                     matchfun=function(d)d$S==d$lR,
                     constants=c(sv=log(1),v=log(1),v_CIincongruent=log(1)))
# Choose a prior assuming no difference in s between match and mismatch
pmean0 <- c(B=log(2),B_Espeed=log(1),B_lRright=log(1),
            v_lMd=2,'v_lMd:CIincongruent'=0,t0=log(.2),
            s_lMd=log(1),A=log(.5))

# With the aid of plot_prior we choose standard deviations so that the prior
# distributions cover a reasonable range.
psd0 <-  c(1,.5,.5,.5,.5,.5,.5,.5)
priorLBA0 <- prior(designLBA0,pmean=pmean0,psd=psd0,type="single")
plot(priorLBA0,designLBA0,layout=c(2,5))
sLBA0 <-  make_emc(dat,designLBA0,type="single",rt_resolution=.05,prior=priorLBA0)
sLBA0 <- fit(sLBA0)
# save(sLBA0,file="samples/sLBA0.RData")

ppLBA0 <- predict(sLBA0)

# d) The old model wins.
compare(list(LBA=sLBA,LBA0=sLBA0))

#  e) Both minD and Dmean indicate slightly better fit for the old model.
#     Graphically there is little to choose between them.
ppLBA <- predict(sLBA)
par(mfrow=c(2,8))
plot_density(dat,ppLBA0,factors=c("S","CI"))
plot_density(dat,ppLBA,factors=c("S","CI"))

# f) As for the RDM, posterior correlations are less extreme, and updating is
# improved. This pattern gives greater confidence when interpreting parameter
# estimates.
pairs_posterior(sLBA0)
pairs_posterior(sLBA)

plot_pars(sLBA0,layout=c(2,5))
plot_pars(sLBA,layout=c(2,5))

# 4) Compare the mapped parameters for the new models that are analogous.
# a) Non-decision time
# b) Rate variability
# c) Thresholds
# d) Mean rates
# e) How does the picture of performance provided by each model differ?

# 4) Note that t0 can be quantitatively compared between models as it is on the
#    same scale. For the other parameters relative sizes of effects can be
#    compared.
credint(sRDM0,map=TRUE)
credint(sLBA0,map=TRUE)

# a) t0 is ~ .02s quicker in the RDM than LBA
# b) In both models rate variability (s/sv) is smaller for the matching
#    accumulator, but to a greater degree (in a ratio sense) for the LBA
# c) There is little evidence of response bias and thresholds are lower in speed
#    than accuracy, more so for the RDM (in a ratio sense)
# d) Match rates are lower for incongruent than congruent, slightly more so for
#    the LBA (in  ratio sense)
# e) Overall the two models present a quite consistent qualitative picture, with
#    some quantitative differences for t0 and sv


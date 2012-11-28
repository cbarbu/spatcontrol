source("spatcontrol.R",chdir=TRUE)
graphics.off()
sampled<-get.sampled()
traces(sampled[-1,])
estimates<-posteriors(sampled[-1,]

cofs<-get.cofactors()
traces(cofs[-1,])
estimates<-c(estimates,posteriors(cofs[-1,]))

betas<-get.betas()
traces(betas[-1,])
estimates<-c(estimates,posteriors(betas[-1,]))
# tb<-summary.spatcontrol(estimates=estimates)

# estimate post-spraying
dbFitted<-read.csv("dbFitted.csv")
ypsamples<-get.field("ypsamples.txt")
ProbaInfByHouse<-apply(ypsamples,2,mean)
ProbaPostSprayBH<-ProbaInfByHouse*(1-dbFitted$observed)
totalInfPostSpray<-sum(ProbaPostSprayBH)


## estimate post spraying:
naivePostSpray<-count(dbFitted$observed==0)*count(dbFitted$positive==1)/count(dbFitted$observed==1)

cat("Naive:",naivePostSpray,"houses infested post spraying\n")
cat("GMRFpino Model:",
totalInfPostSpray,"houses infested post spraying")
cat("Or",naivePostSpray/totalInfPostSpray,"times less\n")





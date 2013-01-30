source("spatcontrol/spatcontrol.R",chdir=TRUE)
graphics.off()
sampled<-get.sampled()
traces(sampled[-1,])
estimates<-posteriors(sampled[-1,])

cofs<-get.cofactors()
if(!is.null(cofs)){
	traces(cofs[-1,])
	estimates<-c(estimates,group.posteriors(cofs[-1,]))
}

betas<-get.betas()
if(!is.null(betas)){
traces(betas[-1,])
estimates<-c(estimates,group.posteriors(betas[-1,]))
}
tb<-summary.spatcontrol(estimates=estimates)

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

# convergence<-cb.diag(betas)

# field analysis
dev.new()
par(mfrow=c(2,2))
plot_reel(dbFitted$X,dbFitted$Y,2*dbFitted$positive+(1-dbFitted$observed),base=0,top=2,main="Data")
us<-getField("usamples.txt")
if(!is.null(us)){
umean<-apply(us,2,mean)
plot_reel(dbFitted$X,dbFitted$Y,umean,base=-2,top=2,main="u mean")
}

yps<-getField("ypsamples.txt")
if(!is.null(yps)){
ypmean<-apply(yps,2,mean)
plot_reel(dbFitted$X,dbFitted$Y,ypmean,base=0,top=1,main="yp mean")
}

os<-getField("osamples.txt")
if(!is.null(os)){
omean<-apply(os,2,mean)
plot_reel(dbFitted$X,dbFitted$Y,omean,base=-2,top=2,main="o mean")
}





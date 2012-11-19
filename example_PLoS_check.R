source("spatcontrol.R")
graphics.off()
samples<-trace.mcmc()
estimates<-posteriors.mcmc(samples=samples)

tb<-summary.spatcontrol(estimates=estimates)


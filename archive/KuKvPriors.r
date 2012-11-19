library("pscl")
KInfLim<-0.1
KSupLim<-20
Kshape<-0.001# 0.5,0.75,1.1
Kscale<-1000 # 100,10,2
aK<-Kshape
bK<-1/Kscale # bK is the rate parameter = 1/scale parameter for the gamma

par(mfrow=c(2,2))
xabs<-seq(sqrt(1/KSupLim)/5,sqrt(1/KInfLim)*5,0.01)
# xabs<-seq(0.01,100,0.01)
plot(xabs,pigamma(xabs^2,aK,bK),col="blue",type="l",xlab="sd",ylab="CDF",ylim=c(0,1))
abline(h=c(0.025,0.975));abline(v=c(sqrt(1/KSupLim),sqrt(1/KInfLim)))
abline(v=c(1,10))
plot(xabs,densigamma(xabs^2,aK,bK),col="blue",type="l",xlab="sd",ylab="PDF")
abline(v=c(1,10))

xabs<-seq(KInfLim/5,KSupLim*5,0.01)
plot(xabs,pgamma(xabs,shape=Kshape,scale=Kscale),col="blue",type="l",xlab="K",ylab="CDF",ylim=c(0,1))
abline(h=c(0.025,0.975))
abline(h=c(0.025,0.975));abline(v=c(KSupLim,KInfLim))
plot(xabs,dgamma(xabs,shape=Kshape,scale=Kscale),col="blue",type="l",xlab="K",ylab="PDF")


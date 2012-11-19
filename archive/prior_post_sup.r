library(locfit)
## T
par(mfcol=c(2,3))
if(use.streets){
	#pdf("Tfit.pdf",width=4,height=4)
	#par(mar=c(4,4,0.5,0.5))
Tfit<-locfit(~sampled[(beginEstimate):(nbsimul),1],xlim=c(0,15),alpha=0.5)
plot(Tfit,main="autocorrelation ratio accross streets\n versus inside block",xlab="T")
xabs<-seq(0.01,15,0.01)
lines(xabs,lik.T(xabs,mT,sdlT,log=FALSE),type="l",lty=2)
# dev.off()
}

	
#pdf("ffit.pdf",width=4,height=4)
#	par(mar=c(4,4,0.5,0.5))
ffit<-locfit(~sampled[(beginEstimate):(nbsimul),3],xlim=c(1,150),alpha=1)
plot(ffit,main="kernel slope",xlab="f")
xabs<-seq(0.01,500,0.01)
lines(xabs,lik.f(xabs,mf,sdlf,log=FALSE),type="l",lty=2)
# dev.off()

Kufit<-locfit(~sampled[(beginEstimate):(nbsimul),5],xlim=c(0,150),alpha=1)
plot(Kufit,main="Spatial precision",xlab="Ku",xlim=c(0,10))
xabs<-seq(0.01,500,0.01)
lines(xabs,dgamma(xabs,shape=Kushape,scale=Kuscale),type="l",lty=2)
# lines(xabs,dgamma(xabs,shape=1.1,scale=0.5),type="l",lty=2)

Kvfit<-locfit(~sampled[(beginEstimate):(nbsimul),10],xlim=c(0,150),alpha=1)
plot(Kvfit,main="Non-Spatial precision",xlab="Kv",xlim=c(0,10))
xabs<-seq(0.01,500,0.01)
lines(xabs,dgamma(xabs,shape=Kvshape,scale=Kvscale),type="l",lty=2)
# lines(xabs,dgamma(xabs,shape=1.1,scale=0.5),type="l",lty=2)

if(use.insp){
	nbbetas<-length(betas[,1])

	plot(c(0,1),c(0,15),main="inspectors detection rates",xlab="beta",type="n")
	for(insp in 1:length(betas[1,])){
		bfit<-locfit(~betas[beginEstimate:nbbetas,insp],xlim=c(0,1),alpha=0.5)
		lines(bfit)
	}
	xabs<-seq(0,1,0.01)
	lines(xabs,dbeta(xabs,abeta,bbeta),type="l",lty=2,col=3)
}

if(use.cofactors){
	xabs<-seq(-2,2,0.01)
	plot(xabs,dnorm(xabs,0,1/sqrt(Kc)),xlab="c.val",ylab="Density",type="l",lty=2,ylim=c(0,6),main="cofactors")
	for(cof in 1:length(c.vals[1,])){
		cfit<-locfit(~c.vals[beginEstimate:nbsimul,cof],xlim=c(-15,15),alpha=1)
		lines(cfit)
	}
}
printdev(device=pdf,"prior_post_sup.pdf")




par(mfrow=c(4,2))
minx<-min(c(u,v,w))
maxx<-max(c(u,v,w))

#### all splited
## basic histogram
hist(v,xlim=c(minx,maxx),col="grey",main="risks induced for one simulation")
hist(u,add=TRUE,col="black",density=15,angle=-45)
hist(c.comp,add=TRUE,col="black",density=50,angle=-45)

hist(est.v,col="grey",xlim=c(minx,maxx),main="estimated risks induced")
hist(est.u,add=TRUE,col="black",density=15,angle=-45)
hist(c.map%*%est.c.val,add=TRUE,col="black",density=50,angle=-45)

library(locfit)

## locfit
vfit<-locfit(~v,xlim=c(minx,maxx),alpha=0.5)
plot(vfit,main="comparison of u,v,c for specific values",xlab="probit risk",xlim=c(-5,5),lty=2,ylim=c(0,1))

ufit<-locfit(~u,xlim=c(minx,maxx),alpha=0.5)
lines(ufit,lty=1)

cfit<-locfit(~c.comp,xlim=c(minx,maxx),alpha=0.5)
lines(cfit,lty=3)

vfit<-locfit(~est.v,xlim=c(minx,maxx),alpha=0.5)
plot(vfit,main="comparison of u,v,c for specific values",xlab="probit risk",xlim=c(-5,5),lty=2,ylim=c(0,1))

ufit<-locfit(~est.u,xlim=c(minx,maxx),alpha=0.5)
lines(ufit,lty=1)

cfit<-locfit(~(c.map%*%est.c.val),xlim=c(minx,maxx),alpha=0.5)
lines(cfit,lty=3)


#### spatial versus non-spatial
## basic histogram
hist(v+c.comp,xlim=c(minx,maxx),col="grey",main="risks induced for one simulation")
hist(u,add=TRUE,col="black",density=15,angle=-45)

hist(est.v+c.map%*%est.c.val,col="grey",xlim=c(minx,maxx),main="estimated risks induced")
hist(est.u,add=TRUE,col="black",density=15,angle=-45)

## locfit
totnonspat<-(v+c.comp)
vfit<-locfit(~totnonspat,xlim=c(minx,maxx),alpha=0.5)
plot(vfit,main="comparison of u,v a specific value",xlab="probit risk",xlim=c(minx,maxx),lty=2,ylim=c(0,1))

ufit<-locfit(~u,xlim=c(minx,maxx),alpha=0.5)
lines(ufit,lty=1)

totnonspat<-(est.v+c.map%*%est.c.val)
vfit<-locfit(~totnonspat,xlim=c(minx,maxx),alpha=0.5)
plot(vfit,main="comparison of u,v estimated",xlab="probit risk",xlim=c(minx,maxx),lty=2)

ufit<-locfit(~est.u,xlim=c(minx,maxx),alpha=0.5)
lines(ufit,lty=1)

printdev(device=pdf,"comparison_risks.pdf")

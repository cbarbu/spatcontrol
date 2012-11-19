
samplepriorKu <- function(dim,Q,u,Kushape,Kuscale) {
	pos.shape <- (Kushape);
	pos.scale <- (Kuscale^(-1))^(-1);
	Ku <- rgamma(n=1, shape=pos.shape, scale=pos.scale);
	return(Ku);
}

check_one<-function(Q,u,Kushape,Kuscale){
	nbsamp<-10000
	dimension<-length(u)
	xabs<-seq(0.0001,1,0.01)
	plot(xabs,dgamma(xabs,shape=Kushape,scale=Kuscale),main="Ku prior",xlab="Ku",ylab="Density",type="l")

	sampledKu<-rep(0,nbsamp)
	for(i in 1:nbsamp){
		sampledKu[i]<-samplepriorKu(dimension,Q,u,Kushape,Kuscale);
	}
	hist(sampledKu)
}
par(mfrow=c(2,3))
hist(u)

Kushape<-0.75
Kuscale<-0.1
check_one(Q,u,Kushape,Kuscale);

Kushape<-0.75
Kuscale<-0.01
check_one(Q,u,Kushape,Kuscale);




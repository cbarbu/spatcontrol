## DIC calculus (see Spiegelhalter2002)
source("DeltaSampling.r")
if(!INTERMEDIARY){
try(source("estimated.txt"),silent=TRUE)
}
source("make_mean_field.r")
if(INTERMEDIARY){
	# if(!exists("betas")){
		betas<-read.table("betasamples.txt");
	# }
	est.beta<-apply(betas,2,mean)
}

# we are in fact only interested by the influence of f/T but filtering out everything else
# may introduce a bad bias so we calculate with everything
Q<-QfromfT(Dmat,AS,SB,meanf,meanT);
meanKu<-mean(smallsampled[,5]);
# Dmean<-llh.ugivQ(dimension,est.u,Q,meanKu)+llh.zgivw(est.w,zpos,zneg,inspector%*%est.beta)+sum(dnorm(v,0,Kv,log=TRUE))
# meanD<-mean(sampled[,4]+sampled[,7]+sampled[,13]+sampled[,14]+sampled[,15]);
# ptheta<-meanD-Dmean;
# DIC <- ptheta+meanD;

# as the variables of interest here are the variables describing the spatial relationship, 
# we use a partial likelihood in the calculus of the DIC, restricted to the description of the data given the mean spatial parameters: f,T,Ku,u
# see 
# Dmean<- -2*(llh.ugivQ(dimension,est.u,Q,meanKu)+llh.zgivw(est.w,zpos,zneg,inspector%*%est.beta))
# meanD<-mean(-2*(smallsampled[,4]+smallsampled[,7]));
# if we condiser that D(theta)=-2*log(y|theta)+c, c being a constant that disappear in the comparisons
# and that y depends directly on w and beta then we can simplify the expression of D:
Dmean<- -2*(llh.zgivw(est.w,zpos,zneg,inspector%*%est.beta))
meanD<-mean(-2*(smallsampled[,7]));
ptheta<-meanD-Dmean;
DIC <- ptheta+meanD;
cat("\nDIC (partial):",DIC,"ptheta:",ptheta,"meanD:",meanD,"Dmean:",Dmean,"\n");



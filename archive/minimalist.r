## parameters
max.simul<-100000

## import data
data <- read.csv("DB_mm_blocks_05May2011.csv",header=TRUE);

## libraries
library(spam);

## sampling functions
sample_u <- function(dimension,Q,Ku,y,cholQ=NULL){
	R <- Ku*Q+diag.spam(1,dimension);
	center <- y;
	u <- rmvnorm.canonical(n=1, b=center, Q=R,Rstruct=cholQ);
	
	return(drop(u));
}
QfromfT<-function(Dmat,AS,SB,f=f.r,T=T.r,addeps=epsilon,origin=NULL,kern=Kernel,dimension=dim(AS)[1]){
	Q<- SpecificMultiply.spam(1/T,-Kernel(T,Dmat,f),SB);#the precision matrix of u 
	rsumQ <- -Q %*% rep(1,dimension); # the precision for each yi
	rsumQ<-rsumQ+addeps # NB we want addeps to be in the normalisation so that we don't have too heavy outliers, addeps is the variance of isolated houses and locically is included in the normalization of the isolation component
	# norm_fact<- mean(1/rsumQ) # the global variance of y
	# cat("norm_fact (mean)",norm_fact,"\n");
	diag.spam(Q) <- rsumQ;

	# norm_fact<- median(1/rsumQ) # the global variance of y
	# diag.spam(Q) <- rsumQ+addeps/norm_fact;
	# cat("norm_fact",norm_fact,"\n");

	# Q<-Q*norm_fact;

	return(Q);
}

## init
num.simul<-1

## main loop
while (num.simul < max.simul) {
	# sample fields/parameters
	u<-sample_u(dimension,Q,K,y,cholQ); # sample the spatial component
	Ku<-sampleKu(dimension,Q,u,Kushape,Kuscale); 
	w<-u # only the spatial component in the risk
	y<-sample_y_direct(w,zpos,zneg,zNA,bivect); # sample the infestation (continuous)
	yprime <- (y>0); # transform the continous infestation into a binary field

	# get likelihood
	LLHz<-llh.zgivw(w,zpos,zneg,bivect);
	LLHu<-llh.ugivQ(dimension,u,Q,Ku,cholQ=cholQ)
	llhKu<-dgamma(Ku,shape=Kushape,scale=Kuscale,log=TRUE)
	LLHTotal<-LLHz+LLHu+llhKu;
	cat("LLHTotal:",LLHTotal,"LLHu",LLHu,"LLHz",LLHz,"llhKu",llhKu,"\n")
	num.simul<-num.simul+1
}


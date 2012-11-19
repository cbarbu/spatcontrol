# this is the gibbs/MH hasting sampler, on both generated or true data:
graphics.off()

### general settings
name<-"trueDataMediumChangedPriorKuKvNoNANoInsp"
use.generated <- FALSE	# TRUE use generated data  FALSE use real data
use.yprime<-FALSE;
use.v<-TRUE;
use.insp<-FALSE;
use.NA<-FALSE;
use.cofactors<-TRUE;
mh.cof<-TRUE;
# for the parameters of the generated data refer to pseudo_data_generation.r
nbsimul <- 200000;
freqsave=10;

set.seed(500)
## set the name
if(use.generated){racgen<-"gen"}else{racgen<-"true"}
if(use.yprime){racy<-"yprime"}else{racy<-"ydirect"}
if(use.v){racv<-"v"}else{racv<-"uv"}
if(use.insp){racinsp<-"Insp"}else{racinsp<-"NoInsp"}

source("pseudo_data_generation.r")
name<-paste(name,racgen,racy,racv,racinsp,"tr",threshold,nbsimul,sep="_")

library("spam")
spam.options(cholsymmetrycheck=FALSE, safemode=c(FALSE,FALSE,FALSE)) ## to speed up things once everything is ok, check same results
# powerboost() ## not sure it is usefull after spam.options
## prep data 
if(use.generated){
	zpos <- which(z.r==1);
	zneg <- which(z.r==0);
	zNA <- which(z.r==9);
}else{
	z.r<-data$status
	graphics.off()
	## initialisation specific to true data
	zpos <- which(data$status==1);
	zneg <- which(data$status==0);
	zNA <- which(data$status==9);
	dev.new()
	par(mfrow=c(1,4));
	select<-c(zneg,zpos)
	plot_reel(data$easting[select],data$northing[select],2*data$status[select]-1,main="data")
}

dev.new()

# priors and sampling settings
nparam<-3
if(use.v){
	nparam=nparam+1;
}
if(use.cofactors){
	nparam=nparam+1;
}
if(use.insp){
	nparam=nparam+1;
}
par(mfrow=c(2,nparam))
# # f with gamma law:
# mf<-0.05
# vf<- 0.01
# shf<-(2+mf^2/vf)+sqrt(mf^2/vf * (4+mf^2/vf))
# scf<-mf/(shf-1)
# xabs<-seq(0,0.2,0.0001)
# plot(xabs,dgamma(xabs,shape=shf,scale=scf),main="f prior",xlab="f",ylab="density")
# f with log normal prior
mf<-1+log(20);
sdlf<-1 # 2
lik.f<-function(f,mf,sdlf,log=TRUE){
	LLH<-dlnorm(f,mf,sdlf,log=log)
	return(LLH);
}
xabs<-seq(-1,7,0.1)
# plot(xabs,lik.f(exp(xabs),mf,sdlf),main="f prior",xlab="f",ylab="density")
plot(xabs,lik.f(exp(xabs),mf,sdlf,log=FALSE),main="f prior",xlab="f (log)",ylab="density",type="l")
logsdfprop<-0.1

# T
mT<-1;
sdlT<-1;
# lik.T<-function(x,...){return(rep(1,length(x)))}; # the uniform prior
lik.T<-lik.f # the log normal prior of f
xabs<-seq(-5,5,0.1)
# plot(xabs,lik.T(exp(xabs),mT,sdlT),main="T prior",xlab="T (log)",ylab="density (log)")
plot(xabs,lik.T(exp(xabs),mT,sdlT,log=FALSE),main="T prior",xlab="T (log)",ylab="density",type="l")
logsdTprop<-0.1

# c
sdc.val<-0.1/nbfact.gen;

# Ku
Kushape <- 2
Kuscale <- 1;
xabs<-seq(0.0001,10,0.1)
plot(xabs,dgamma(xabs,shape=Kushape,scale=Kuscale),main="Ku prior",xlab="Ku",ylab="Density",type="l")

if(use.v){
# Kv
Kvshape<-2
Kvscale<-1
xabs<-seq(0.0001,10,0.1)
plot(xabs,dgamma(xabs,shape=Kvshape,scale=Kvscale),main="Kv prior",xlab="Kv",ylab="Density",type="l")
K.hyper <- c(Kushape,Kuscale,Kvshape,Kvscale); 
}

if(use.cofactors){
# Kc
Kcshape<-2
Kcscale<-1
xabs<-seq(0.0001,10,0.1)
plot(xabs,dgamma(xabs,shape=Kcshape,scale=Kcscale),main="Kc prior",xlab="Kc",ylab="Density",type="l")
}

# b inspectors found rate when infested
## high prior
# abeta <- 35;
# bbeta <- 4;
## uniform prior
abeta <- 18; ## allow to have the mean at 0.9
bbeta <- 2; ## allow to have a "direct" slot to 0 in 1
if(use.insp){
	xabs<-seq(0,1,0.01)
	plot(xabs,dbeta(xabs,abeta,bbeta),main="Betas prior",xlab="Beta",ylab="Density",type="l")
}

# likelihood functions
# general function for likelihood of u given Q, used for Delta and f MH
lik.ugivQ<-function(dimension,u,Q,Ku,cholQ=NULL){
	# cholQ<-chol(Q);
	if(is.null(cholQ)){
		cholQ<-chol.spam(Q);
	}else{
		cholQ<-update.spam.chol.NgPeyton(cholQ,Q);
	}
	exp_part<- Ku*(t(u)%*%Q%*%u); 
	det_part<- dimension*log(Ku)+2*determinant(cholQ)$modulus
	logpi<-dimension*log(2*pi);
	LLH= -1/2*(exp_part-det_part+logpi); 
	cat("LLH",LLH,"exppart:",exp_part,"det_part",det_part,"logpi",logpi,"\n");
	return(LLH);
}
lik.ygivw<-function(y,w){
	## return the loglikelihood of y given w
	LLH<- -1/2 * t(y-w)%*%(y-w);
	return(LLH);
}
llh.zgivy<-function(y,zpos,zneg,bivect){
	probypos<-(y>0)
	LLHpos<-sum(log(probypos[zpos]*bivect[zpos]));
	LLHneg<-sum(log(probypos[zneg]*(1-bivect[zneg])+(1-probypos[zneg])));
	LLH<-LLHpos+LLHneg;
	return(LLH);
}
llh.zgivw<-function(w,zpos,zneg,bivect){
	probypos<-1-pnorm(0,w,1)
	LLHpos<-sum(log(probypos[zpos]*bivect[zpos]));
	LLHneg<-sum(log(probypos[zneg]*(1-bivect[zneg])+(1-probypos[zneg])));
	LLH<-LLHpos+LLHneg;
	return(LLH);
}
general.likelihood<-function(w,u,Ku,y,v=NULL,c.compprop=NUL,c.valprop.c=NULL,Kv=1,Kc=1){
	 LLH=llh.ugivQKu+llh.KugivpKu+llh.TgivpT+llh.fgivpf;
	 if(use.insp){
		 LLH=llh.zgivy+llh.ygivw+llh.ugivQKu+llh.KugivpKu+llh.TgivpT+llh.fgivpf;
		 print("not yet implemented\n")
		 break()
	 }else{
		 LLH=LLH+llh.zgivw(w,zpos,zneg,bivect);
	 }
	 if(use.v){
		 LLH=LLH+llh.vgivKv+llh.KvgivpKv
	 }
	 if(use.cofactors){
		 LLH=LLH+llh.c(c.valprop,c.compprop,Kc,w-c.compprop,y); 
	 }

	return(LLH);
}
# 	# but if we want to do model selection, namely on T, the kernel and cofactors
# 	# we may not want to introduce the probability of the model itself but only
# 	# the probability of the data and intermediate state given the model
# 	# or even say that u and v are parameters of the model, and then we only
# 	# look at the likelihood of z|w
# 	return(LLH);
# }


## sampling functions
acceptDelta<-{};
sample_Delta <- function(u,K,Delta,sdDprop,aD,bD,f,Q,LLHu,cholQ=NULL,Dmat=NULL){

	# Q<-QfromDelta(Delta,dist_mat,AS,f=f);

	Delta_prop=rtnorm(1,mean=Delta,sd=sdDprop,lower=0,upper=threshold);
	out<-QfromDelta(Delta_prop,dist_mat,AS,f=f);
	Qprop<-out[[1]]
	Dmatprop<-out[[2]]

	hasting_term=dtnorm(Delta,mean=Delta_prop,sd=sdDprop,lower=0,upper=threshold,log=TRUE)-dtnorm(Delta_prop,mean=Delta,sd=sdDprop,lower=0,upper=threshold,log=TRUE);
	LLHuprop<-lik.ugivQ(dimension,u,Qprop,K[1])
	LLHproposal <- LLHuprop+lik.Delta(Delta_prop,aD,bD); 
	LLH <- 	       LLHu+lik.Delta(Delta,aD,bD);
	lnr <- LLHproposal-LLH+hasting_term;
	cat("Delta:",Delta,"LLH:",LLH,"Delta prop:",Delta_prop,"LLHproposal:",LLHproposal,"h_t",hasting_term,"lnr",lnr);

	if(lnr>=log(runif(1))) {
		Delta <- Delta_prop;
		Q<-Qprop;
		LLHu<-LLHuprop;
		Dmat<-Dmatprop;
		acceptDelta<<-c(acceptDelta,1);
		cat(" accept 1\n");
	}else{
		acceptDelta<<-c(acceptDelta,0);
		cat("accept 0\n");
	}
	return(list(Delta,Q,LLHu,cholQ,Dmat));
}
library("stats")
acceptf<-{};
sample_f <- function(u,K,Delta,logsdfprop,f,mf,sdlf,Q,LLHu,AS,SB,cholQ=NULL,Dmat=NULL){
	# Q<-QfromDelta(Delta,dist_mat,AS,f=f);

	f_prop<-rlnorm(1,log(f),logsdfprop);

	Qprop<-QfromfT(Dmat,AS,SB,f=f_prop,T=T);
	# Qprop<-Qfromf(Dmat,f=f_prop);

	hasting_term=dlnorm(f,log(f_prop),logsdfprop,log=TRUE)-dlnorm(f_prop,log(f),logsdfprop,log=TRUE);

	LLHuprop<-lik.ugivQ(dimension,u,Qprop,K[1],cholQ=cholQ)

	llhfprop<-lik.f(f_prop,mf,sdlf);
	llhf<-lik.f(f,mf,sdlf);

	LLHproposal <- LLHuprop+llhfprop; 
	LLH <- 	       LLHu+llhf;
	lnr <- LLHproposal-LLH+hasting_term;

	cat("f:",f," LLH:",LLH,"(",LLHu,"+",llhf,") f prop:",f_prop," LLHproposal:",LLHproposal,"(",LLHuprop,"+",llhfprop,") h_t:",hasting_term," lnr",lnr,sep="");

	if(lnr>=log(runif(1))) {
		f <- f_prop;
		Q<-Qprop;
		LLHu<-LLHuprop;
		acceptf<<-c(acceptf,1);
		cat(" accept 1\n");
	}else{
		acceptf<<-c(acceptf,0);
		cat("accept 0\n");
	}
	return(list(f,Q,LLHu,cholQ));
}
acceptT<-{};
sample_T <- function(u,K,f,T,logsdTprop,mT,sdT,Q,LLHu,AS,SB,cholQ=NULL,Dmat=NULL){
	T_prop<-rlnorm(1,log(T),logsdTprop);

	Qprop<-QfromfT(Dmat,AS,SB,f=f,T=T_prop);

	hasting_term=dlnorm(T,log(T_prop),logsdTprop,log=TRUE)-dlnorm(T_prop,log(T),logsdTprop,log=TRUE);

	LLHuprop<-lik.ugivQ(dimension,u,Qprop,K[1],cholQ=cholQ);
	llhTprop<-lik.T(T_prop,mT,sdlT);
	llhT<-lik.T(T,mT,sdlT);
	LLHproposal <- LLHuprop+llhTprop; 
	LLH <- 	       LLHu+llhT;
	lnr <- LLHproposal-LLH+hasting_term;

	cat("T:",T," LLH:",LLH,"(",LLHu,"+",llhT,") T prop:",T_prop," LLHproposal:",LLHproposal,"(",LLHuprop,"+",llhTprop,") h_t:",hasting_term," lnr",lnr,sep="");

	if(lnr>=log(runif(1))) {
		T <- T_prop;
		Q<-Qprop;
		LLHu<-LLHuprop;
		acceptT<<-c(acceptT,1);
		cat(" accept 1\n");
	}else{
		acceptT<<-c(acceptT,0);
		cat("accept 0\n");
	}
	return(list(T,Q,LLHu,cholQ));
}
sampleDeltaf<-function(dimension,u,K,Delta,logsdfprop,sdDprop,f,mf,sdlf,aD,bD,Q,LLHu,cholQ=NULL,Dmat=NULL){
	Delta_prop=rtnorm(1,mean=Delta,sd=sdDprop,lower=0,upper=threshold);
	out<-QfromDelta(Delta_prop,dist_mat,AS,f=f);
	Qprop<-out[[1]]
	Dmatprop<-out[[2]]

	hasting_term=dtnorm(Delta,mean=Delta_prop,sd=sdDprop,lower=0,upper=threshold,log=TRUE)-dtnorm(Delta_prop,mean=Delta,sd=sdDprop,lower=0,upper=threshold,log=TRUE);
	LLHuprop<-lik.ugivQ(dimension,u,Qprop,K[1])
	LLHproposal <- LLHuprop+lik.Delta(Delta_prop,aD,bD); 
	LLH <- 	       LLHu+lik.Delta(Delta,aD,bD);
	lnr <- LLHproposal-LLH+hasting_term;
	cat("Delta:",Delta,"LLH:",LLH,"Delta prop:",Delta_prop,"LLHproposal:",LLHproposal,"h_t",hasting_term,"lnr",lnr);

	if(lnr>=log(runif(1))) {
		Delta <- Delta_prop;
		Q<-Qprop;
		LLHu<-LLHuprop;
		Dmat<-Dmatprop;
		acceptDelta<<-c(acceptDelta,1);
		cat(" accept 1\n");
	}else{
		acceptDelta<<-c(acceptDelta,0);
		cat("accept 0\n");
	}
	return(list(Delta,Q,LLHu,cholQ,Dmat));

}

sample_u <- function(dimension,Q,K,y,cholQ=NULL){
	R <- K[1]*Q+diag.spam(1,dimension);
	center <- y;
	u <- rmvnorm.canonical(n=1, b=center, Q=R,Rstruct=cholQ);
	
	return(drop(u));
}

source("sample_y_direct.r")

sampley <- function(dim,u,yprime) {
  center <- u;
  lwbd <- rep(0,dim);
  lwbd[(yprime==0)] <- -Inf;
  upbd <- rep(0,dim);
  upbd[(yprime==1)] <- Inf;
  y <- rtnorm(n=dim, mean=center, sd=1, lower=lwbd, upper=upbd);
  return(y);
}
sampley2 <- function(dim,u,yprime) {
  center <- u;
  lwbd <- rep(0,dim);
  lwbd[(yprime==0)] <- -Inf;
  upbd <- rep(0,dim);
  upbd[(yprime==1)] <- Inf;
  y[(yprime!=-1)] <- rtnorm(n=dim, mean=center, sd=1, lower=lwbd, upper=upbd);
  y[(yprime==-1)]<-u[(yprime==-1)]; # to be changed to a mean of the other y but for now a good approximation
  return(y);
}

sampleyprime <- function(dim,w,betaprime,zpos,zneg,zNA) {
  yprime <- rep(0,dim);
  yprime[zpos] <- 1;
  p <- (1-betaprime[zneg])*(1-pnorm(q=0, mean=(w[zneg]), sd=1));
  q <- (pnorm(q=0, mean=(w[zneg]), sd=1));
  area<-p+q
  area[area==0]<-1; # avoid NA and associated warnings
  Pzneg <- p/(area);
  PzNA <- 1-pnorm(q=0, mean=(w[zNA]), sd=1);
  yprime[zneg] <- rbinom(length(zneg),1,Pzneg);
  yprime[zNA] <-  rbinom(length(zNA),1,PzNA); # w[zNA]  
  return(yprime);
}
sampleyprime2 <- function(dim,w,betaprime,zpos,zneg,zNA) {
  yprime <- rep(0,dim);

  yprime[zpos] <- 1;

  ## identify which zneg will be considered NA
  ## we should not consider the probability to be positive but only the probability to be bad observed
  ## for the badly observed we put -1 as for the NA and for the well observed we put 0
  ## P(obs-)=(1-beta), P(obs)=beta
  goodobs<-rbinom(length(zneg),1,betaprime);
  goodobs<-goodobs-1; # the good observations are 0 in terms of yprime, the other are -1
  yprime[zneg]<- goodobs;
  yprime[zNA] <- -1;

  return(yprime);
}
sampleK <- function(dim,Q,x,K.hyper) {
  Ku.a <- K.hyper[1];
  Ku.b <- K.hyper[2];
  Kv.a <- K.hyper[3];
  Kv.b <- K.hyper[4];
  u <- x[1:dim];
  v <- x[(dim + (1:dim))] - u;
  u.pshape <- (0.5*(dim-1) + Ku.a);
  u.pscale <- (0.5*as.numeric(u %*% (Q %*% u)) + Ku.b^(-1))^(-1);
  v.pshape <- (0.5*dim + Kv.a);
  v.pscale <- (0.5*as.numeric(v %*% v) + Kv.b^(-1))^(-1);
  Ku <- rgamma(n=1, shape=u.pshape, scale=u.pscale);
  Kv <- rgamma(n=1, shape=v.pshape, scale=v.pscale);
  K <- c(Ku,Kv);
  return(K);
}

## ok for both Kv and Kc (and could be used for Ku with dim<-dimension-1)
sampleKg <- function(dim,Q,field,K.shape,K.scale,n=1) {
  pshape <- (0.5*dim + K.shape);
  pscale <- (0.5*as.numeric(t(field) %*%Q%*%field) + K.scale^(-1))^(-1);
  K <- rgamma(n=n, shape=pshape, scale=pscale);
  return(K);
}
sampleKu <- function(dim,Q,u,Kushape,Kuscale) {
	pos.shape <- (0.5*(dim-1) + Kushape);
	pos.scale <- (0.5*as.numeric(u %*% Q %*% u) + Kuscale^(-1))^(-1);
	Ku <- rgamma(n=1, shape=pos.shape, scale=pos.scale);
	return(Ku);
}

makeR <- function(dim,Q,K) {
  Ku <- K[1];
  Kv <- K[2];
  R <- Ku*Q;
  diag.spam(R) <- diag.spam(R) + Kv;
  R <- cbind(R, diag.spam(-1*Kv,dim,dim));
  R <- rbind(R, cbind(diag.spam(-1*Kv,dim,dim), diag.spam((Kv+1), dim, dim)));
  return(R);
}
# makeRc <- function(Q,K,Qc){
#   Ku <- K[1];
#   Kv <- K[2];
#   Kc <- K[3];
#   KQc<-Kc*Qc;
#   R <- Ku*Q+KQc;
#   R <- cbind.spam(R, -KQc);
#   # dimension<-nrow(Qc)
#   # zero<-spam(0,dimension,dimension)
#   R <- rbind.spam(R, cbind.spam(-KQc, KQc+diag.spam(1,nrow(Qc))));
#   return(R);
# }
# makeRc<- function(Q,K,Qc){
# 	# make Rc using kronecker product
# 	# substantively faster
# 	Ku <- K[1];
# 	Kv <- K[2];
# 	Kc <- K[3];
# 	KQc<-Kc*Qc;
# 	R <- Ku*Q;
# 	dimension<-nrow(Qc)
# 	zero<-spam(0,dimension,dimension)
# 	R <- cbind.spam(R, zero);
# 	R <- rbind.spam(R, cbind.spam(zero, diag.spam(1,nrow(Qc))));
# 
# 	K<-kronecker.spam(matrix(c(1,-1,-1,1),2),KQc)
# 	R <-R+K;
# 	return(R);
# }
makeRc<- function(Q,K,Qc){
	# make Rc using kronecker product
	# substantially faster
	Ku <- K[1];
	Kv <- K[2];
	Kc <- K[3];
	R <- Ku*Q;
	dimension<-nrow(Qc)
	zero<-spam(0,dimension,dimension)
	R <- cbind.spam(R, zero);
	R <- rbind.spam(R, cbind.spam(zero, diag.spam(1,nrow(Qc))));

	R <-R+Kc*kQc;
	return(R);
}
baseRc<- function(Q,K,Qc){
	Ku <- K[1];
	Kv <- K[2];
	Kc <- K[3];
	R <- Ku*Q;
	dimension<-nrow(Qc)
	zero<-spam(0,dimension,dimension)
	R <- cbind.spam(R, zero);
	R <- rbind.spam(R, cbind.spam(zero, zero));

	return(R);
}
basediag1<-function(dimension){
	zero<-spam(0,dimension,dimension)
	Rdiag1 <- cbind.spam(zero, zero);
	Rdiag1 <- rbind.spam(Rdiag1, cbind.spam(zero, diag.spam(1,nrow(Qc))));
	return(Rdiag1);
}

updateRc<- function(Q,RQc,QentriesInR,diag1entriesInR,K){
	# make Rc using kronecker product
	# substantially faster
	Ku <- K[1];
	Kc <- K[3];
	RKQc <- Kc*RQc;

	RKQc@entries[QentriesInR]<-Ku*Q@entries+RKQc@entries[QentriesInR];
	RKQc@entries[diag1entriesInR]<-RKQc@entries[diag1entriesInR]+1;
	return(RKQc);
}

# # for(i in 1:1000){
# system.time(R <- makeRc(Q,K,Qc))
# # }

samplexcof<-function(dimension,Q,K,y,Qc,cholR){
	x <- rnorm(n=(2*dimension), mean=0, sd=1);
	center <- c(rep(0,dimension), y);
	# R <- makeRc(Q,K,Qc);
	R <-updateRc(Q,RQc,QentriesInR,diag1entriesInR,K) ## consistently faster than makeRc
	if(is.null(cholR)){
		cholR <- chol.spam(R, memory=list(nnzcolindices=4e6),);
	}else{
		cholR <- update.spam.chol.NgPeyton(cholR,R);
	}

	center <- backsolve(cholR, forwardsolve(cholR, center));
	x <- backsolve(cholR,x);
	x <- x + center;
	cholR<<-cholR;
	return(x);
}
##  as we quickly run into memory troubles with the gibbs sampler
##  a metropolis hasting version is welcome
# pi(c | y,u,v)  \propto pi(y | c.val,u,v)pi(c.val)
#                \propto N(c.comp+u+v,1) N(c.val,Kc)
llh.c<-	function(c.val,c.comp,Kc,wnoc,y){
	w<-wnoc+c.comp;
	LLH<-lik.ygivw(y,w)# +sum(dnorm(c.val,mean=0,sd=Kc,log=TRUE));
	return(LLH);
}
acceptc.val<-{}
mhsamplec<-function(c.val,c.comp,sdc.val,Kc,wnoc,y){
	nbc<-length(c.val);
	c.valprop<-rnorm(nbc,c.val,sdc.val);
	# hasting_term=dnorm(c.val,log(c.valprop),logsdc.val,log=TRUE)-dnorm(c.valprop,log(c.val),logsdc.val,log=TRUE);
	c.compprop<-drop(c.map%*%c.valprop);
	LLHproposal <- llh.c(c.valprop,c.compprop,Kc,wnoc,y); 
	LLH <-llh.c(c.val,c.comp,Kc,wnoc,y);
	lnr <- LLHproposal-LLH# +hasting_term;
	cat("c.val:",c.val,"LLH:",LLH,"c.valprop:",c.valprop,"LLH proposal:",LLHproposal,"lnr",lnr);

	if(lnr>=log(runif(1))) {
		c.val <- c.valprop;
		c.comp<-c.compprop;
		acceptc.val<<-c(acceptc.val,1);
		cat(" accept 1\n");
	}else{
		acceptc.val<<-c(acceptc.val,0);
		cat("accept 0\n");
	}

	return(list(c.val,c.comp));
}
# the problem is more the sampling of u and v through the gibbs sampler
# I think that we can simply use y-c.comp as the center
# using x <- samplex(dimension,Q,K,y-c.comp,cholR);

samplex <- function(dim,Q,K,y,cholR=NULL) {
  x <- rnorm(n=(2*dim), mean=0, sd=1);
  center <- c(rep(0,dim), y);
  R <- makeR(dim,Q,K);
  if(is.null(cholR)){
	  cholR <- chol.spam(R, memory=list(nnzcolindices=4e6),);
  }else{
	  cholR <- update.spam.chol.NgPeyton(cholR,R);
  }

  center <- backsolve(cholR, forwardsolve(cholR, center));
  x <- backsolve(cholR,x);
  x <- x + center;
  cholR<<-cholR;
  return(x);
}

samplebeta <- function(zpos,zneg,matrix,yprime,a,b) {
	yp.positive <- yprime;
	yp.negative <- yprime; 
	yp.positive[-zpos] <- 0; # z- turn y+ into 0 a.pos being then the sum of the y+ well observed
	yp.negative[-zneg] <- 0; # z+ turn y+ into 0 a.neg being then the sum of the y+ badly observed 
	# the NA are not apparent here
	a.pos <- (as.vector(t(matrix) %*% yp.positive) + a); 
	b.pos <- (as.vector(t(matrix) %*% yp.negative) + b); 
	beta <- rbeta(n=ncol(matrix), shape1=a.pos, shape2=b.pos); 
	return(beta);
}

samplebeta2 <- function(zpos,zneg,matrix,yprime,a,b) {
	yprime<-abs(yprime) # when inspected, all the NA are y+
	yp.positive <- yprime;
	yp.negative <- yprime; 
	yp.positive[-zpos] <- 0; # z- turn y+ into 0 a.pos being then the sum of the y+ well observed
	yp.negative[-zneg] <- 0; # z+ turn y+ into 0 a.neg beint then the sum of the y+ badly observed 
	# the NA are not apparent here
	a.pos <- (as.vector(t(matrix) %*% yp.positive) + a); 
	b.pos <- (as.vector(t(matrix) %*% yp.negative) + b); 
	beta <- rbeta(n=ncol(matrix), shape1=a.pos, shape2=b.pos); 
	return(beta);
}

## general
# starting values
Delta=0;
f=20;
T<-1;
Ku<-1;
Kv<-1;
Kc<-1;

K<-c(Ku,Kv,Kc);
Q<-QfromfT(dist_mat,AS,SB,f,T);
# out<-QfromDelta(Delta,dist_mat,AS,f);
# Q<-out[[1]]
# Dmat<-out[[2]]
cholQ<-chol(Q);
LLHTu<-0
LLHfu<-0

u<-rep(0,dimension);
y<-rep(0,dimension);
yprime <- (y>0);
w <- rnorm(dimension,u,sqrt((K[2])^(-1)));
x <- c(u, w);
beta<-rep(1,inspector@dimension[2])
###### start on true values 
Delta=0;
f=f.r;
T<-T.r;
Ku<-Ku.r;
Kv<-Kv.r;
cholR<-NULL;
c.val<-rep(0,nbfact.gen)
c.comp<-c.map%*%c.val

## initialisation for updateRc
if(!mh.cof && use.cofactors){
kQc<-kronecker.spam(matrix(c(1,-1,-1,1),2),Qc)
R<-makeRc(Q,K,Qc);
R0<-makeRc(Q,K,Qc);
R0@entries<-rep(0,length(R0@entries));
RQ<-baseRc(Q,K,Qc);
RQ0<-R0+RQ;
QentriesInR<-which(RQ0@entries!=0);
RQc<-R0+kQc;
Rdiag1<-basediag1(dimension);
Rdiag1<-R0+Rdiag1;
diag1entriesInR<-which(Rdiag1@entries==1);
rm(R0,RQ,kQc,RQ0,Rdiag1);
system.time(R<-updateRc(Q,RQc,QentriesInR,diag1entriesInR,K)) ## only to test beforehand
}

u<-u.r;
y<-y.r;
yprime <- (y>0);
w <- w.r;
x <- c(u, w);
beta<-beta.r
######
cholR<-NULL;

write.table(t(u), "usamples.txt", sep="\t",col.names=FALSE,row.names=FALSE)

## sampler
nbtraced=12;
sampled<-as.matrix(mat.or.vec(nbsimul+1,nbtraced));
sampled[1,1]<-T;
sampled[1,2]<-0;
sampled[1,3]<-f;
sampled[1,4]<-0;
sampled[1,5]<-Ku;
sampled[1,6]<-llh.zgivy(y,zpos,zneg,bivect);
sampled[1,7]<-llh.zgivw(u,zpos,zneg,bivect);
sampled[1,8]<-lik.ygivw(y,u);
sampled[1,9]<-0;
sampled[1,10]<-Kv;
sampled[1,11]<-mean(u);
sampled[1,12]<-Kc;
c.val<-rep(0,nbfact.gen)

write.table(t(sampled[1,]), "sampled.txt", sep="\t",col.names=FALSE,row.names=FALSE)
write.table(t(u), "usamples.txt", sep="\t",col.names=FALSE,row.names=FALSE)
write.table(t(w), "wsamples.txt", sep="\t",col.names=FALSE,row.names=FALSE)
if(use.cofactors){
	write.table(t(c.val), "cofactors.txt", sep="\t",col.names=FALSE,row.names=FALSE)
}
if(use.insp){
	write.table(t(beta), "betasamples.txt", sep="\t",col.names=FALSE,row.names=FALSE)
}

starter<-1

dev.new()
nbploted<-2
if(use.cofactors){
	nbploted<-nbploted+1
}
if(use.v){
	nbploted<-nbploted+1
}
if(use.insp){
	nbploted<-nbploted+1
}
par(mfcol=c(nbploted,5))
# Rprof();
for (i in starter:nbsimul) {
	cat("\n",i,", out of:",nbsimul," ");
	
	# K <- sampleK(dimension,Q,x,K.hyper);
	# x <- samplex(dimension,Q,K,y,cholR);
	# y <- sampley(dimension,x,yprime);
	# yprime <- sampleyprime(dimension,x,betaprime,zpos,zneg,zNA);
	# beta <- samplebeta(zpos,zneg,inspector,yprime,a,b);

	if(use.yprime){
		yprime<-sampleyprime(dimension,w,bivect,zpos,zneg,zNA);
		y<-sampley(dimension,w,yprime);
	}else{
		y<-sample_y_direct(w,zpos,zneg,zNA,bivect);
		yprime <- (y>0);
		# cat("likugivQ after y sampling",lik.ugivQ(dimension,u,Q,K[1]),"\n");
		# y<-y.r;

	}
	if(use.cofactors && use.v){
		if(mh.cof){
			x <- samplex(dimension,Q,K,y-c.comp,cholR);
			u<-x[1:dimension];
			v<-x[dimension+(1:dimension)]-u;
			c.all<-mhsamplec(c.val,c.comp,sdc.val,Kc,x[dimension+(1:dimension)],y);
			c.val<-c.all[[1]]
			c.comp<-c.all[[2]]
			# Kc <- sampleKg(dimension,Qc,c.comp,Kcshape,Kcscale);
			K <- sampleK(dimension,Q,x,K.hyper);
			Ku<-K[[1]]
			Kv<-K[[2]]
			w<-x[dimension+(1:dimension)]+c.comp;
		}else{
			print("not implemented yet")
			break()
		}

	}else if(use.cofactors && !use.v){
		x <- samplexcof(dimension,Q,K,y,Qc,cholR);
		u<-x[1:dimension];
		w<-x[dimension+(1:dimension)];
		Kc <- sampleKg(dimension,Qc,w-u,Kcshape,Kcscale);
		Ku <- sampleKu(dimension,Q,u,Kushape,Kuscale); 
		K<-c(Ku,1,Kc);
	}else if(!use.cofactors && use.v){
		x <- samplex(dimension,Q,K,y,cholR);
		u<-x[1:dimension];
		w<-x[dimension+(1:dimension)]
		K <- sampleK(dimension,Q,x,K.hyper);
		cat(" Ku:",K[1]," Kv:",K[2]);
	}else{
		u<-sample_u(dimension,Q,K,y,cholQ);
		w<-u
		# cat("likugivQ after u sampling",lik.ugivQ(dimension,u,Q,K[1]),"\n");
		##

		Ku<-sampleKu(dimension,Q,u,Kushape,Kuscale); 
		K[1]<-Ku;
		cat(" Ku:",Ku);
	}
	cat("\nKu",Ku,"Kv",Kv,"Kc",Kc,"\n");
	LLHKuu<-lik.ugivQ(dimension,u,Q,K[1],cholQ=cholQ)

	out<-sample_f(u,K,Delta,logsdfprop,f,mf,sdlf,Q,LLHKuu,AS,SB,cholQ=cholQ,Dmat=Dmat); 
	f<-out[[1]];
	Q<-out[[2]];
	LLHfu<-out[[3]];
	cholQ<-out[[4]];

	out<-sample_T(u,K,f,T,logsdTprop,mT,sdT,Q,LLHfu,AS,SB,cholQ=cholQ,Dmat=Dmat);
	T<-out[[1]];
	Q<-out[[2]];
	LLHTu<-out[[3]];
	cholQ<-out[[4]];
	
	if(use.insp){
		beta <- samplebeta(zpos,zneg,inspector,yprime,abeta,bbeta);
		bivect <- as.vector(inspector %*% beta);
	cat("beta (",mean(beta),",",sd(beta),")");
	}

	# if((i)%%20==0){
	# 	# adapt sampling of Delta
	# 	rateaccept<-mean(tail(acceptT,20))
	# 	if(rateaccept<0.333){
	# 		logsdTprop<-0.9*logsdTprop;
	# 		cat("update of logsdTprop to:",logsdTprop);
	# 	}else if(rateaccept>0.666){
	# 		logsdTprop<-1.1*logsdTprop;
	# 		cat("update of logsdTprop to:",logsdTprop);
	# 	}
	# 	# adapt sampling of f
	# 	rateaccept<-mean(tail(acceptf,20))
	# 	if(rateaccept<0.333){
	# 		logsdfprop<-0.9*logsdfprop;
	# 		cat("update of logsdfprop to:",logsdfprop);
	# 	}else if(rateaccept>0.666){
	# 		logsdfprop<-1.1*logsdfprop;
	# 		cat("update of logsdfprop to:",logsdfprop);
	# 	}
	# }
	LLHyw<-lik.ygivw(y,w);
	LLHy<-llh.zgivy(y,zpos,zneg,bivect);
	LLH<-llh.zgivw(w,zpos,zneg,bivect);
	cat("LLHyw:",LLHyw,"LLHy:",LLHy,"LLH:",LLH,"mu:",mean(u),"sdu:",sd(u));

	sampled[i+1,1]<-T;
	sampled[i+1,2]<-LLHTu;
	sampled[i+1,3]<-f;
	sampled[i+1,4]<-LLHfu;
	sampled[i+1,5]<-K[1];
	sampled[i+1,6]<-LLHy;
	sampled[i+1,7]<-LLH;
	sampled[i+1,8]<-LLHyw;
	sampled[i+1,9]<-i;
	sampled[i+1,10]<-K[2];
	sampled[i+1,11]<-mean(u);
	sampled[i+1,12]<-Kc;

	if(i%%freqsave==0 || i==(nbsimul)){
		write.table(t(u), "usamples.txt", sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
		write.table(t(w), "wsamples.txt", sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
		
		write.table(sampled[(i+1-(freqsave-1)):(i+1),], "sampled.txt",append=TRUE, sep="\t",col.names=FALSE,row.names=FALSE)
		plot_reel(data$easting,data$northing,y,main="y")
		plot_reel(data$easting,data$northing,u,main="u")
		if(use.cofactors){
			# plot_reel(data$easting,data$northing,w-u,main="cofactors")
			c.val<-c.map.plus%*%(w-u)
			if(use.generated){
				plot(c.val~c.val.r)
			}else{
				plot(c.val)
			}
			write.table(t(c.val), "cofactors.txt", sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
			cat("\n c.val:",c.val,"(",c.val.r,")")
		}
		if(use.v){
			v<-x[dimension+(1:dimension)]-u;
			plot_reel(data$easting,data$northing,v,main="v")
		}
		if(use.insp){
			write.table(t(beta), "betasamples.txt", sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
			plot_reel(data$easting,data$northing,bivect*2-1,main="insp")
		}
	}
}
dev.print(device=png,file="finalmaps.pdf",width=1000,height=500)
nbsimul<-(i-1)
simulname<-paste(name,"_tr",threshold,"_T",T.r,"_f",f.r,"Ku",Ku.r,"_",nbsimul,sep="")
dump("sampled",file=paste(simulname,".txt",sep=""))
dump("u",file=paste(simulname,"_finalu.txt",sep=""))
dump("y",file=paste(simulname,"_finaly.txt",sep=""))
sampled<-sampled[1:nbsimul,];
# Rprof(NULL);
cat("\n")
### traces
dev.new()
par(mfcol=c(2,nparam+1))
plot(c(1,nbsimul),log10(c(min(sampled[-1,1],T.r),max(sampled[-1,1],T.r))),xlab="T",ylab="log10(T)",type="n")
lines(log10(sampled[-1,1]),type="l")
abline(h=log10(T.r))
plot(sampled[-1,2],xlab="LLHT",type="l")
plot(c(1,nbsimul),log10(c(min(sampled[-1,3],f.r),max(sampled[-1,3],f.r))),xlab="f",ylab="log10(f)",type="n")
lines(log10(sampled[-1,3]),xlab="f",type="l")
if(use.generated){
abline(h=log10(f.r))
}
# plot(sampled[-1,4],xlab="LLHf",type="l")
plot(sampled[-1,5],xlab="Ku",type="l")
if(use.generated){
abline(h=Ku.r)
}
plot(sampled[-1,7],xlab="LLH z given w",type="l")
if(use.v){
	plot(sampled[-1,10],xlab="Kv",type="l")
	if(use.generated){
	abline(h=Kv.r)
	}
}
if(use.cofactors){
	c.vals<-read.table("cofactors.txt")
	if(use.generated){
		plot(c(1,length(c.vals[,1])),c(min(c.vals,c.val.r),max(c.vals,c.val.r)),type="n")
		for(i in 1:nbfact.gen){
			abline(h=c.val.r[i],col=i);	
		}
	}else{
	plot(c(1,nbsimul),c(min(c.vals),max(c.vals)),type="n")
	}
	for(i in 1:nbfact.gen){
		lines(c.vals[,i],col=i)
	}

	# plot(sampled[-1,12],xlab="Kc",type="l")
	# if(use.generated){
	# abline(h=Kc.r)
	# }
}
	plot(sampled[-1,6],xlab="LLH z given y",type="l")
plot(sampled[-1,11],xlab="mean(u)",type="l")
abline(h=mu.r)
plot(sampled[-1,8],xlab="LLH y given w",type="l")

dev.print(device=png,paste("traces.png",sep=""),width=800,height=400)

### posteriors
dev.set(3)
par(mfrow=c(2,nparam))
# xabs<-seq(-1,3,0.05)
# plot(xabs,lik.f(10^(xabs),mf,sdlf,log=FALSE),main="f prior",xlab="f (log10)",ylab="density",type="l")
xabs<-seq(1,500,1)
plot(xabs,lik.f(xabs,mf,sdlf,log=FALSE),main="f prior",xlab="f (log10)",ylab="density",type="l")
xabs<-seq(-2,2,0.1)
plot(xabs,lik.T(10^(xabs),mT,sdlT,log=FALSE),main="T prior",xlab="T (log10)",ylab="density",type="l")
xabs<-seq(0,30,0.1)
plot(xabs,dgamma(xabs,shape=Kushape,scale=Kuscale),main="Ku prior",xlab="Ku",ylab="Density",type="l")
if(use.cofactors){
	plot(xabs,dgamma(xabs,shape=Kcshape,scale=Kcscale),main="Kc prior",xlab="Kc",ylab="Density",type="l")
}
if(use.v){
plot(xabs,dgamma(xabs,shape=Kvshape,scale=Kvscale),main="Kv prior",xlab="Kv",ylab="Density",type="l")
}
if(use.insp){
	xabs<-seq(0,1,0.01)
	plot(xabs,dbeta(xabs,abeta,bbeta),main="Betas prior",xlab="Beta",ylab="Density",type="l")
}

hist(sampled[(nbsimul/2):nbsimul,3],main=paste("f posterior (mean=",signif(mean(sampled[(nbsimul/2):(nbsimul),3]),4),")",sep=""))
abline(v=mean(sampled[(nbsimul/2):(nbsimul),3]))
if(use.generated){
abline(v=f.r,col=4)
}
hist(sampled[(nbsimul/2):(nbsimul),1],main=paste("T posterior (mean=",signif(mean(sampled[(nbsimul/2):(nbsimul),1]),4),")",sep=""))
abline(v=mean(sampled[(nbsimul/2):(nbsimul),1]))
if(use.generated){
abline(v=T.r,col=4)
}
hist(sampled[(nbsimul/2):(nbsimul),5],main=paste("Ku posterior (mean=",signif(mean(sampled[(nbsimul/2):(nbsimul),5]),4),")",sep=""))
abline(v=mean(sampled[(nbsimul/2):(nbsimul),5]))
if(use.generated){
abline(v=Ku.r,col=4)
}
if(use.v){
	hist(sampled[(nbsimul/2):(nbsimul),10],main=paste("Kv posterior (mean=",signif(mean(sampled[(nbsimul/2):(nbsimul),10]),4),")",sep=""))
	abline(v=mean(sampled[(nbsimul/2):(nbsimul),10]))
	if(use.generated){
		abline(v=Kv.r,col=4)
	}
}
if(use.insp){
	betas<-read.table("betasamples.txt");
	nbbetas<-length(betas[,1])
	meanbetas<-apply(betas[(nbbetas/2):(nbbetas),],2,mean)
	hist(as.vector(as.matrix(betas[(nbbetas/2):nbbetas,])),main=paste("Betas posterior (mean=",signif(mean(meanbetas),4),")",sep=""),xlim=c(0,1))
	abline(v=mean(meanbetas))
	if(use.generated){
		abline(v=beta.r,col=4)
	}
}
dev.print(device=pdf,paste("prior_post.pdf",sep=""))

## final fields
dev.new()
ncolvisu=3;
if(use.insp){
	ncolvisu=ncolvisu+1;
}
par(mfrow=c(2,ncolvisu))
visudata<-z.r
visudata[visudata==9]<-0
plot_reel(data$easting,data$northing,2*visudata-1,main="data")
if(use.NA){
dataPaleNA<-z.r
dataPaleNA[zNA]<-0.5
plot_reel(data$easting,data$northing,2*dataPaleNA-1,main="data with pale NA")
}
visupseudodata<-2*generate_z(y,bivect)-1;
visupseudodata[data$status==9 &visupseudodata==1]<-0
plot_reel(data$easting,data$northing,visupseudodata,main="generated z final")
if(use.insp){
plot_reel(data$easting,data$northing,bivect*2-1,main="beta final")
}
plot_reel(data$easting,data$northing,y,main="y final")
plot_reel(data$easting,data$northing,w,main="w final")
plot_reel(data$easting,data$northing,w-u,main="c final")
plot_reel(data$easting,data$northing,u,main="u final")

dev.print(device=pdf,paste("map.pdf",sep=""))

## hist LLH
dev.set(4)
hist(sampled[(nbsimul/2):(nbsimul),7],main="LLH z given u")
abline(v=mean(sampled[(nbsimul/2):(nbsimul),7]))
hist(sampled[(nbsimul/2):(nbsimul),4],main="LLHf (u given Q)")

### correlations
dev.new()
smallsampled<-sampled[(nbsimul/2):nbsimul,]
smallsampled<-smallsampled[seq(1,length(smallsampled[,1]),freqsave),]
smallsampled<-as.data.frame(smallsampled)
namessmallsampled<-c("T","LLHTu","f","LLHfu","Ku","LLHy","LLH","LLHyw","Kv","i","mu","Kc")
par(mfcol=c(3,nparam))

colnames(smallsampled)<-namessmallsampled

plot(smallsampled$LLH~smallsampled$f)
if(use.generated){
abline(v=f.r,col=4)
}
plot(smallsampled$LLH~smallsampled$T)
if(use.generated){
abline(v=Delta.r,col=4)
}
plot(smallsampled$LLH~smallsampled$Ku)
if(use.generated){
abline(v=Ku.r,col=4)
}

plot(smallsampled$LLHfu~smallsampled$f)
if(use.generated){
abline(v=f.r,col=4)
}
plot(smallsampled$LLHTu~smallsampled$T)
if(use.generated){
abline(v=Delta.r,col=4)
}
plot(smallsampled$LLH~smallsampled$Ku)
if(use.generated){
abline(v=Ku.r,col=4)
}

plot(smallsampled$T~smallsampled$f)
if(use.generated){
abline(v=f.r,col=4)
abline(h=T.r,col=4)
}
plot(smallsampled$T~smallsampled$Ku)
if(use.generated){
abline(h=T.r,col=4)
abline(v=Ku.r,col=4)
}
plot(smallsampled$f~smallsampled$Ku)
if(use.generated){
abline(h=f.r,col=4)
abline(v=Ku.r,col=4)
}
if(use.cofactors){
	plot(smallsampled$T~smallsampled$Kc)
	if(use.generated){
		abline(v=Kc.r,col=4)
		abline(h=T.r,col=4)
	}
	plot(smallsampled$f~smallsampled$Kc)
	if(use.generated){
		abline(h=f.r,col=4)
		abline(v=Kc.r,col=4)
	}
	plot(smallsampled$Ku~smallsampled$Kc)
	if(use.generated){
		abline(h=Ku.r,col=4)
		abline(v=Kc.r,col=4)
	}
}
if(use.v){
plot(smallsampled$T~smallsampled$Kv)
if(use.generated){
abline(v=Kv.r,col=4)
abline(h=T.r,col=4)
}
plot(smallsampled$f~smallsampled$Kv)
if(use.generated){
abline(h=f.r,col=4)
abline(v=Kv.r,col=4)
}
plot(smallsampled$Ku~smallsampled$Kv)
if(use.generated){
abline(h=Ku.r,col=4)
abline(v=Kv.r,col=4)
}
}
dev.print(device=png,paste("LLH_and_cov.png",sep=""),width=600,height=600)

## mean kernel
dev.new()
meanT<-mean(smallsampled$T)
meanf<-mean(smallsampled$f)
xabs<-seq(0,threshold);
plot(xabs,exp(-xabs/meanf),type="l",ylim=c(0,1),col=4,main=paste("mean kernel (T: ",meanT,", f: ",meanf,") \n blue:intrablocks ; red: interblocks", seq=""))
lines(xabs,meanT*exp(-xabs/meanf),col=2)
dev.print(device=pdf,"mean_kernel.pdf")

## print main results
cat(file="in_brief.txt","mean T:",meanT,"\n");
cat(file="in_brief.txt","mean f:",meanf,"\n",append=TRUE);
cat(file="in_brief.txt","mean Ku:",mean(smallsampled$Ku),"\n",append=TRUE);
if(use.v){
cat(file="in_brief.txt","mean Kv:",mean(smallsampled$Kv),"\n",append=TRUE);
}
cat(file="in_brief.txt","mean final u:",mean(u),"\n",append=TRUE);
if(use.v){
cat(file="in_brief.txt","mean final v:",mean(v),"\n",append=TRUE);
}
cat(file="in_brief.txt","mean final y:",mean(y),"\n",append=TRUE);


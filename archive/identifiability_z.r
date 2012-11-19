# this is the new sampler, on generated data with:
# sampler on Delta and f separatly by MH according to P(u|.)P(.) 
graphics.off()

source("pseudo_data_generation.r")
spam.options(cholsymmetrycheck=FALSE, safemode=c(FALSE,FALSE,FALSE)) ## to speed up things once everything is ok, check same results
dev.new()
source("/home/cbarbu/Documents/recherche/Penn/anal_bug_spread_cross_sectionnal/junction_CB_AH/test_inspectors/sample_y_direct.r")

### general settings
name<-"cosamp_nomedcorrec_stablesdsamp_nodetLLHDf"
nbsimul <- 200000;
freqsave=10

# priors and sampling settings
par(mfrow=c(2,4))
# # f with gamma law:
# mf<-0.05
# vf<- 0.01
# shf<-(2+mf^2/vf)+sqrt(mf^2/vf * (4+mf^2/vf))
# scf<-mf/(shf-1)
# xabs<-seq(0,0.2,0.0001)
# plot(xabs,dgamma(xabs,shape=shf,scale=scf),main="f prior",xlab="f",ylab="density")
# f with log normal prior
mf<-0.05;
sdlf<-100
lik.f<-function(f,mf,sdlf){
	LLH<-dlnorm(f,log(mf),sdlf)
	return(LLH);
}
xabs<-seq(0,0.2,0.0001)
plot(xabs,lik.f(xabs,mf,sdlf),main="f prior",xlab="f",ylab="density")
logsdfprop<-0.07

# Delta
aD<-1   # shape of the gamma
bD<-1 # scale of the gamma
sdDprop <- 5 # standard deviation of the gaussian function to propose Delta in MH process
lik.Delta<-function(Delta,aD,bD,tr=threshold){
	LLH<-dbeta(Delta/threshold,aD,bD,log=TRUE)-log(threshold);
	return(LLH);
}
xabs<-seq(0,threshold,1)
plot(xabs,lik.Delta(xabs,aD,bD,threshold),main="Delta prior",xlab="Delta",ylab="Density")

# Ku
Kushape <- 1
Kuscale <- 1;
xabs<-seq(0.0001,10,0.1)
plot(xabs,dgamma(xabs,shape=Kushape,scale=Kuscale),main="Ku prior",xlab="Ku",ylab="Density")

# b inspectors found rate when infested
## high prior
# abeta <- 35;
# bbeta <- 4;
## uniform prior
abeta <- 1;
bbeta <- 1;
xabs<-seq(0,1,0.01)
plot(xabs,dbeta(xabs,abeta,bbeta),main="Betas prior",xlab="Beta",ylab="Density")

# likelihood functions
# general function for likelihood of u given Q, used for Delta and f MH
lik.ugivQ<-function(dimension,u,Q,Ku,cholQ=NULL){
	cholQ<-chol(Q,Rstruct=cholQ);
	if(is.null(cholQ)){
		cholQ<-chol.spam(Q);
	}else{
		cholQ<-update.spam.chol.NgPeyton(cholQ,Q);
	}
	exp_part<- t(u)%*%Q%*%u; 
	det_part<- 2*determinant(cholQ)$modulus
	logpi<-dimension*log(2*pi);
	LLH= -1/2*(Ku*exp_part-(dimension*log(Ku)+det_part)+logpi); 
	# cat("LLH",LLH,"exppart:",exp_part,"det_part",det_part,"logpi",logpi,"\n");
	return(list(LLH,exp_part,det_part,logpi,cholQ));
}
lik.ygivw<-function(y,w){
	LLH<- -1/2 * t(y-w)%*%(y-w);
	return(LLH);
}
lik.zgivy<-function(y,zpos,zneg,bivect){
	probypos<-(y>0)
	LLHpos<-sum(log(probypos[zpos]*bivect[zpos]));
	LLHneg<-sum(log(probypos[zneg]*(1-bivect[zneg])+(1-probypos[zneg])));
	LLH<-LLHpos+LLHneg;
	return(LLH);
}
lik.zgivw<-function(w,zpos,zneg,bivect){
	probypos<-1-pnorm(0,w,1)
	LLHpos<-sum(log(probypos[zpos]*bivect[zpos]));
	LLHneg<-sum(log(probypos[zneg]*(1-bivect[zneg])+(1-probypos[zneg])));
	LLH<-LLHpos+LLHneg;
	return(LLH);
}


## sampling functions
acceptDelta<-{};
sample_Delta <- function(u,K,Delta,sdDprop,aD,bD,f,Q,LLHu,cholQ=NULL,Dmat=NULL){

	# Q<-QfromDelta(Delta,dist_mat,AS,f=f);

	Delta_prop=rtnorm(1,mean=Delta,sd=sdDprop,lower=0,upper=threshold);
	out<-QfromDelta(Delta_prop,dist_mat,AS,f=f);
	Qprop<-out[[1]]
	Dmatprop<-out[[2]]

	hasting_term=dtnorm(Delta,mean=Delta_prop,sd=sdDprop,lower=0,upper=threshold,log=T)-dtnorm(Delta_prop,mean=Delta,sd=sdDprop,lower=0,upper=threshold,log=T);
	LLHuprop<-lik.ugivQ(dimension,u,Qprop,K[1])
	LLHproposal <- LLHuprop[[1]]+lik.Delta(Delta_prop,aD,bD); 
	LLH <- 	       LLHu[[1]]+lik.Delta(Delta,aD,bD);
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
sample_f <- function(u,K,Delta,logsdfprop,f,mf,sdlf,Q,LLHu,cholQ=NULL,Dmat=NULL){
	# Q<-QfromDelta(Delta,dist_mat,AS,f=f);
	f_prop<-rlnorm(1,log(f),logsdfprop);

	Qprop<-Qfromf(Dmat,f=f_prop);

	hasting_term=dlnorm(f,log(f_prop),logsdfprop,log=T)-dlnorm(f_prop,log(f),logsdfprop,log=T);
	LLHuprop<-lik.ugivQ(dimension,u,Qprop,K[1],cholQ)
	LLHproposal <- LLHuprop+lik.f(f_prop,mf,sdlf); 
	LLH <- 	       LLHu+lik.f(f,mf,sdlf);
	lnr <- LLHproposal-LLH+hasting_term;
	cat("f:",f,"LLH:",LLH,"f prop:",f_prop,"LLHproposal:",LLHproposal,"h_t",hasting_term,"lnr",lnr);

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
acceptDeltaf<-{};
sampleDeltaf<-function(dimension,u,K,Delta,logsdfprop,sdDprop,f,mf,sdlf,aD,bD,Q,LLHu,cholQ=NULL,Dmat=NULL){
	Delta_prop=rtnorm(1,mean=Delta,sd=sdDprop,lower=0,upper=threshold);
	f_prop<-rlnorm(1,log(f),logsdfprop);

	Qprop<-QfromDelta(Delta_prop,dist_mat,AS,f=f_prop)[[1]];

	hasting_term=dtnorm(Delta,mean=Delta_prop,sd=sdDprop,lower=0,upper=threshold,log=T)-dtnorm(Delta_prop,mean=Delta,sd=sdDprop,lower=0,upper=threshold,log=T);
	hasting_term=hasting_term+dlnorm(f,log(f_prop),logsdfprop,log=T)-dlnorm(f_prop,log(f),logsdfprop,log=T);

	LLHu[[2]]<-u%*%Q%*%u
	LLHu[[1]] <- -1/2*(Ku*LLHu[[2]]-(dimension*log(Ku)+LLHu[[3]])+LLHu[[4]]); # only update LLH with Ku

	LLHuprop<-lik.ugivQ(dimension,u,Qprop,K[1])
	LLHproposal <- LLHuprop[[1]]+lik.Delta(Delta_prop,aD,bD)+lik.f(f_prop,mf,sdlf);; 
	LLH <- 	       LLHu[[1]]+lik.Delta(Delta,aD,bD)+lik.f(f,mf,sdlf);
	lnr <- LLHproposal-LLH+hasting_term;
	cat("Delta:",Delta,"f:",f,"LLH:",LLH,"Delta prop:",Delta_prop,"f prop:",f_prop,"LLHproposal:",LLHproposal,"h_t",hasting_term,"lnr",lnr);

	if(lnr>=log(runif(1))) {
		Delta <- Delta_prop;
		f <- f_prop;
		Q<-Qprop;
		LLHu<-LLHuprop;
		acceptDeltaf<<-c(acceptDeltaf,1);
		cat(" accept 1\n");
	}else{
		acceptDeltaf<<-c(acceptDeltaf,0);
		cat("accept 0\n");
	}
	return(list(Delta,f,Q,LLHu));
}

sampleKu <- function(dim,Q,u,Kushape,Kuscale) {
	pos.shape <- (0.5*(dim) + Kushape);
	pos.scale <- (0.5*as.numeric(u %*% Q %*% u) + Kuscale^(-1))^(-1);
	Ku <- rgamma(n=1, shape=pos.shape, scale=pos.scale);
	return(Ku);
}

Id<-diag.spam(1,dimension);
sample_u <- function(dimension,Q,K,y,cholQ=NULL){
	R <- K[1]*Q+Id;
	# center <- y;
	# u <- rmvnorm.canonical(n=1, b=y, Q=R);
	u <- rmvnorm.canonical(n=1, b=y, Q=R,Rstruct=cholQ);
	
	return(drop(u));
}

sample_y<-sample_y_direct; # defined in sample_y_direct

samplebeta <- function(zpos,zneg,matrix,yprime,a,b) {
	yp.positive <- yprime;
	yp.negative <- yprime; 
	yp.positive[-zpos] <- 0; 
	yp.negative[-zneg] <- 0; 
	a.pos <- (as.vector(t(matrix) %*% yp.positive) + a); 
	b.pos <- (as.vector(t(matrix) %*% yp.negative) + b); 
	beta <- rbeta(n=ncol(matrix), shape1=a.pos, shape2=b.pos); 
	return(beta);
}

########### initialisation
# ## initialisation specific to true data
# zpos <- which(data$status==1);
# zneg <- which(data$status==0);
# zNA <- which(data$status==9);

## initialisation specific to generated data
zpos <- which(z.r==1);
zneg <- which(z.r==0);
zNA <- which(z.r==9);

## general
data.insp <- data$collector;
inspectors <- unique(data.insp[c(zpos,zneg)]);
inspector <- matrix(0,dimension,length(inspectors))
for (i in 1:length(inspectors)) {
  inspector[,i] <- (data.insp==inspectors[i]);
}
inspector <- as.spam(inspector);

Delta=0;
f=0.2;
Ku<-10;
K<-c(Ku,Kv.r)
out<-QfromDelta(Delta,dist_mat,AS,f);
Q<-out[[1]]
Dmat<-out[[2]]

u<-rep(0,dimension);
y<-rep(0,dimension);
LLHu<-lik.ugivQ(dimension,u,Q,Ku)
cholQ<-LLHu[[5]];

## sampler
sampled<-as.matrix(mat.or.vec(nbsimul+1,8));
sampled[1,1]<-Delta;
sampled[1,2]<-0;
sampled[1,3]<-f;
sampled[1,4]<-0;
sampled[1,5]<-Ku;
sampled[1,6]<-lik.zgivy(y,zpos,zneg,bivect);
sampled[1,7]<-lik.zgivw(u,zpos,zneg,bivect);
sampled[1,8]<-lik.ygivw(y,u);

dev.new()
par(mfcol=c(2,5))
# Rprof();
for (i in 1:nbsimul) {
	cat("\n",i,", ");
	## test of Ku,Delta,f fitted on z
	y<-sample_y_direct(u,zpos,zneg,zNA,bivect);
	# plot_reel(data$easting,data$northing,y,main="y")
	# cat("likugivQ after y sampling",lik.ugivQ(dimension,u,Q,K[1]),"\n");

	u<-sample_u(dimension,Q,K,y,cholQ);
	# plot_reel(data$easting,data$northing,u,main="u")
	# cat("likugivQ after u sampling",lik.ugivQ(dimension,u,Q,K[1]),"\n");
	##

	Ku<-sampleKu(dimension,Q,u,Kushape,Kuscale); 
	K[1]<-Ku;
	# LLHKuu<-lik.ugivQ(dimension,u,Q,Ku,cholQ)
	cat("Ku:",Ku);

	# out<-sample_f(u,K,Delta,logsdfprop,f,mf,sdlf,Q,LLHKuu,cholQ=cholQ,Dmat=Dmat); 
	# f<-out[[1]];
	# Q<-out[[2]];
	# LLHfu<-out[[3]];
	# cholQ<-out[[4]];


	# out<-sample_Delta(u,K,Delta,sdDprop,aD,bD,f,Q,LLHfu,Dmat=Dmat); 
	# Delta<-out[[1]];
	# Q<-out[[2]];
	# LLHDu<-out[[3]];
	# cholQ<-out[[4]];
	# Dmat<-out[[5]];

	# co sampling Delta and f
	out<-sampleDeltaf(dimension,u,K,Delta,logsdfprop,sdDprop,f,mf,sdlf,aD,bD,Q,LLHu);
	Delta<-out[[1]];
	f<-out[[2]];
	Q<-out[[3]];
	LLHu<-out[[4]];
	cholQ<-LLHu[[5]];

	# if((i)%%20==0){
	# 	# # adapt sampling of Delta
	# 	# rateaccept<-mean(tail(acceptDelta,20))
	# 	# if(rateaccept<0.333){
	# 	# 	sdDprop<-0.9*sdDprop;
	# 	# 	cat("update of sdDprop to:",sdDprop);
	# 	# }else if(rateaccept>0.666){
	# 	# 	sdDprop<-1.1*sdDprop;
	# 	# 	cat("update of sdDprop to:",sdDprop);
	# 	# }
	# 	# # adapt sampling of f
	# 	# rateaccept<-mean(tail(acceptf,20))
	# 	# if(rateaccept<0.333){
	# 	# 	logsdfprop<-0.9*logsdfprop;
	# 	# 	cat("update of logsdfprop to:",logsdfprop);
	# 	# }else if(rateaccept>0.666){
	# 	# 	logsdfprop<-1.1*logsdfprop;
	# 	# 	cat("update of logsdfprop to:",logsdfprop);
	# 	# }
	# 	# adapt sampling of f
	# 	rateaccept<-mean(tail(acceptDeltaf,20))
	# 	if(rateaccept<0.333){
	# 		logsdfprop<-0.9*logsdfprop;
	# 	 	sdDprop<-0.9*sdDprop;
	# 		cat("\nupdate -");
	# 	}else if(rateaccept>0.666){
	# 		logsdfprop<-1.1*logsdfprop;
	# 	 	sdDprop<-1.1*sdDprop;
	# 		cat("\nupdate +");
	# 	}
	# 	cat("sdDprop:",sdDprop,"logsdfprop",logsdfprop,"\n");
	# }
	LLHyu<-lik.ygivw(y,u);
	LLHy<-lik.zgivy(y,zpos,zneg,bivect);
	LLH<-lik.zgivw(u,zpos,zneg,bivect);
	cat("LLHyu:",LLHyu,"LLHy:",LLHy,"LLH:",LLH);

	sampled[i+1,1]<-Delta;
	sampled[i+1,2]<-0;
	sampled[i+1,3]<-f;
	sampled[i+1,4]<-LLHu[[1]];
	sampled[i+1,5]<-Ku;
	sampled[i+1,6]<-LLHy;
	sampled[i+1,7]<-LLH;
	sampled[i+1,8]<-LLHyu;
}
simulname<-paste(name,"_tr",threshold,"_D",Delta.r,"_f",f.r,"Ku",Ku.r,"_",nbsimul,sep="")
dump("sampled",file=paste(simulname,".txt",sep=""))
# Rprof(NULL);
cat("\n")
dev.new()
par(mfrow=c(2,4))
plot(sampled[-1,1],xlab="Delta",type="l")
abline(h=Delta.r)
plot(sampled[-1,2],xlab="LLHDelta",type="l")
plot(sampled[-1,3],xlab="f",type="l")
abline(h=f.r)
plot(sampled[-1,4],xlab="LLHf",type="l")
plot(sampled[-1,5],xlab="Ku",type="l")
abline(h=Ku.r)
plot(sampled[-1,6],xlab="LLH z given y",type="l")
plot(sampled[-1,7],xlab="LLH z given u",type="l")
plot(sampled[-1,8],xlab="LLH y given u",type="l")
dev.print(device=png,paste(simulname,".png",sep=""),width=800,height=400)

dev.set(3)
hist(sampled[(nbsimul/2):nbsimul,3],main="f posterior")
abline(v=mean(sampled[(nbsimul/2):(nbsimul),3]))
abline(v=f.r,col=4)
hist(sampled[(nbsimul/2):(nbsimul),1],main="Delta posterior")
abline(v=mean(sampled[(nbsimul/2):(nbsimul),1]))
abline(v=Delta.r,col=4)
hist(sampled[(nbsimul/2):(nbsimul),5],main="Ku posterior")
abline(v=mean(sampled[(nbsimul/2):(nbsimul),5]))
abline(v=Ku.r,col=4)
hist(sampled[(nbsimul/2):(nbsimul),7],main="LLH z given u")
abline(v=mean(sampled[(nbsimul/2):(nbsimul),7]))
dev.print(device=pdf,paste(simulname,"prior_post.pdf",sep=""))

dev.set(2)
plot_reel(data$easting,data$northing,u,main="u final")
plot_reel(data$easting,data$northing,y,main="y final")
dev.print(device=pdf,paste(simulname,"map.pdf",sep=""))

dev.set(4)
hist(sampled[(nbsimul/2):(nbsimul),4],main="LLHf (u given Q)")

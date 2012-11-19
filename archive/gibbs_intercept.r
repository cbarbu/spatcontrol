### Entering data into R

### Parameters
name<-"MM_exp_eps0.01_intercept_fT"
n <- 100; # number of simulations
set.seed(333141)
threshold <- 50;
default_kern="exp" # can be "gaussian" or "exp" or "cauchy" or "geometric", see prep_.r for the implementation

# components in the model
use.f<-TRUE # will adapt the spatial auto-correlation shape
use.streets<-TRUE;
use.NormQ <- "mean"; # if Q is normalized ("med", "mean" or FALSE)

# adaptive sampling
highAcceptRate<-0.4 # according to Gelman1996
lowAcceptRate<-0.15 # idem

# initialization values
t <- 3;
ku <- 0.1;
f <- 20
T <- 1

### Hyperparameters

t <- c(0,0.01); # intercept of spatial component
k <- c(1,100);  # hyperparameters of Ku
b <- c(18,2); 	# hyperparameters of the detection quality
fprior=15 # prior characteristic distance in meters
sdlf<-1 # standard deviation in logarithm for characteristic distance prior
logsdfprop<-0.1 # initial standard deviation for the MH sampling of f
mT<-1; # prior for effect of streets
sdlT<-1; # standard deviation in logarithm for effect of streets 
logsdTprop<-0.1 # initial standard deviation for the MH sampling of T
epsilon<-0.01

# setwd("C:/Users/ahon/Desktop/Project/Current");

data <- read.csv("mm.csv",header=T);
# colnames(data) <- c("status","collector","easting","northing");

### Packages

library(spam);
library(msm);
library(lattice);
library(sp);
trellis.par.set(sp.theme());

# WARNING some functions may be redundant in following source, so order 
# of calling matters
source("spam_complement.r")
source("DeltaSampling.r")
source("functions_intercept.r");

### Data cleaning

## Read data

data <- data[order(data$easting,data$northing),];
data <- data[-c(538,809,10714),];
dimension <- nrow(data);

## z

pos <- which(data$status==1);
neg <- which(data$status==0);
na <- which(data$status==9);
z <- list(pos=pos, neg=neg, na=na);

visudata<-data$status
visudata[visudata==9]<-0.5


## distance matrix construction
spam.options(nearestdistnnz=c(9058076,400))
Dmat <- nearest.dist(x=data[,c("easting","northing")], y=NULL, method="euclidian", delta=threshold, upper=NULL);

## Nota:
# diag.spam() is quite slow we then save the position of the diagonal terms in the 
# spam entries of Dmat:
diagDmat<-which(Dmat@entries==0) 
# it allows to do Dmat@entries[diagDmat]<-rep(1,length(diagDmat))
# in spite of diag.spam(Dmat)<-1
# 6 times faster, can be used also for the diagonal of Q

## streets matrices construction
spam.options(nearestdistnnz=c(13764100,400))
SB <- nearest.dist(x=cbind(data$block_num,rep(0,length(data$block_num))), method="euclidian", upper=NULL,delta=0.1)
SB@entries<-rep(1,length(SB@entries))
dmt<-Dmat
dmt@entries<-rep(1,length(dmt@entries))# [dmt@entries!=0]<-1 # 1 only when Dmat not 0
SB@entries<-rep(1,length(SB@entries))
SB<-logic_and_spam(SB,dmt);
AS<-dmt-SB; # get 1 whereever the distances matrix is defined(under threshold) and not same block
AS<-as.spam(AS)

## Covariance matrix construction
# affectation of the kernel
if(default_kern == "exp"){
	Kernel<-expKernel;
}else if(default_kern == "gaussian"){
	Kernel<-gaussianKernel;
}else if(default_kern == "cauchy"){
	Kernel<-cauchyKernel;
}else if(default_kern == "geometric"){
	Kernel<-cauchyKernel;
}else{
	stop("Error, unknown kernel:",default_kern,"Please change default_kern")
}
Q<-QfromfT(Dmat,AS,SB,f,T);
# W<-Dmat
# W@entries[which(W@entries!=0)] <- 1/W@entries[which(W@entries!=0)];
# rsumW <- as.vector(W %*% rep(1,dimension));
# Q <- -1*W;
# diag.spam(Q) <- rsumW+epsilon;

cholQ<-chol.spam(Q)
AutMats <- list(AS=AS,SB=SB,Dmat=Dmat)

## Inspector matrix

inspectors <- unique(data$collector[c(z$pos,z$neg)]);
I <- matrix(0,length(inspectors),dimension)
for (i in 1:length(inspectors)) {
  I[i,] <- (data$collector==inspectors[i]);
}
I <- as.spam(I);

## technical transformations of the parameters
mf<-1+log(fprior);
hyper <- list(t=t, k=k, b=b,f=c(mf,sdlf),T=c(mT,sdlT));

### Gibbs sampler

## Initialize Sampler
##  technical initialization
# # random
# gen <- generator(dim=dimension, hyper=hyper, Q=Q, I=I, z=z, s=1, ku=ku, t=t);
# state <- gen$state;

# fix initialization
u<-rep(0,dimension)
y1<- rep(0,dimension)
sy0<- rep(0,dimension)
beta <- rep(1,nrow(I))
LLHu<-llh.ugivQ(dimension,u,Q,ku,cholQ=cholQ)
state <- list(ku=ku, x=c(u,t), sy0=sy0, y1=y1, beta=beta,f=f,T=T,Q=Q,cholQ=cholQ,LLHu=LLHu);
MHsd<-list(T=logsdTprop,f=logsdfprop)

R <- makeRuint(dim=dimension, hyper=hyper, Q=Q, ku=ku);
cholR <- chol.spam(R, memory=list(nnzcolindices=4e6));

adaptOK<-FALSE

## Storage

ksamples <- c();
tsamples <- c();
lsamples <- c();
msamples <- c();
fsamples <- c();
Tsamples <- c();
usamplesSum <- rep(0,dimension)
sy0samplesSum <- rep(0,dimension)
y1samplesSum <- rep(0,dimension)
bsamples <- matrix(0, nrow=n, ncol=length(inspectors));
mI<-0

## Sampling loop
Rprof()
for (iter in 1:n) {
  smpl <- gibbs(dim=dimension, hyper=hyper, I=I, z=z, s=1, state=state,MHsd=MHsd,AutMats=AutMats);

  if((iter)%%20==0 && !adaptOK){
	  MHsd <- updateMHparam(MHsd);
	  adaptOK<-MHsd$OK
  }
  ksamples[iter] <- smpl$ku;
  fsamples[iter] <- smpl$f;
  Tsamples[iter] <- smpl$T;
  tsamples[iter] <- smpl$x[(dimension+1)];
  usamplesSum <- usamplesSum+smpl$x[1:dimension]
  sy0samplesSum <- sy0samplesSum+smpl$sy0
  y1samplesSum <- y1samplesSum+smpl$y1
  lsamples[iter] <- lik(dim=dimension, hyper=hyper, Q=smpl$Q, I=I, state=smpl);
  # lsamples[iter] <- llh.ugivQ(dimension,u,smpl$Q,smpl$ku);
  # msamples[iter] <- moran(dim=dimension, W=W, y1=(smpl$y1));
  bsamples[iter,] <- smpl$beta;
  state <- smpl;
  cat("n=",iter,"Intercept:", tsamples[iter],"ku=", ksamples[iter],"\n");
  cat("L=", lsamples[iter],"Mean.Beta=", mean(bsamples[iter,]),"\n")
  # mI <- paste("Moran's I=", msamples[iter]);
}
Rprof(NULL)
save.image(file="EndSampleImage.img")

par(mfrow=c(1,4))
plot(ksamples,type="l")
plot(tsamples,type="l")
plot(lsamples,type="l")
plot(apply(bsamples,1,mean),type="l")

dev.new()

par(mfrow=c(2,3))
plot_reel(data$easting,data$northing,visudata,base=0,top=1,main="data")
plot_reel(data$easting,data$northing,smpl$x[1:dimension],base=-1,top=1,main="u")
plot_reel(data$easting,data$northing,smpl$sy0,base=-1,top=1,main="sy0")
plot_reel(data$easting,data$northing,smpl$y1,base=0,top=1,main="y1")
plot_reel(data$easting,data$northing,usamplesSum/n,base=-1,top=1,main="est u")
estimateb<-apply(bsamples,2,mean)
plot_reel(data$easting,data$northing,t(I)%*%estimateb,base=0,top=1,main="est beta")
dev.new()
hist(estimateb)

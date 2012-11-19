### Entering data into R
### parameters
## K
K.hyper <- c(1,1,1,1);
Kushape=K.hyper[1]
Kuscale=K.hyper[2]
Kvshape=K.hyper[3]
Kvscale=K.hyper[4]
p <- 0.5; # the sample parameter

Ku<-1;
Kv<-1;
Delta<-10;
f<-1
inity<- 0

nbsimul <- 1000;
freqsave=10

# #### original data
# # # data <- read.csv("data.csv",header=T);
# # # colnames(data) <- c("locality","status","collector","easting","northing");
#  
# ## good data
# data <- read.csv("DB_mm_blocks_05May2011.csv",header=T);
# # be sure to change in original file :
# # localidad -> locality
# # colector -> collector
# # ; -> ,
# 
# # cleaning of inspectors
# data$collector[data$collector==" Jorge Ampuero"] <- "Jorge Ampuero"
# data$collector[data$collector==" Jorge A"] <- "Jorge Ampuero"
# data$collector[data$collector=="Jose Velasquez "] <- "Jose Velasquez"
# data$collector[data$collector==" Hugo Vilcahuaman"] <- "Hugo Vilcahuaman"
# data$collector[data$collector=="Julio cesar Condori"] <- "Julio Cesar Condori"
# data$collector[data$collector=="Manuel Tamayo "] <- "Manuel Tamayo"
# 
# data$collector<-factor(data$collector)
# 
# ### Packages for this section
# 
# library(spam);
# powerboost("off") # to be put "on" on production calculus
# 
# ### Data cleaning
# 
# data <- data[order(data$locality, data$easting, data$northing),];
# dimension <- nrow(data);
# 
# #### Data generation
# source("pseudo_data_generation.r")

### sampling functions
sampleK <- function(dim,x,Kushape,Kuscale,Kvshape,Kvscale) {
  u <- x[1:dim];
  eta <- x[(dim + 1:dim)];
  v <- eta - u;
  u.pshape <- (0.5*(dim-1) + Kushape);
  u.pscale <- (0.5*as.numeric(u %*% (Q %*% u)) + Kuscale^(-1))^(-1);
  v.pshape <- (0.5*dim + Kvshape);
  v.pscale <- (0.5*as.numeric(v %*% v) + Kvscale^(-1))^(-1);
  Ku <- 1 # rgamma(n=1, shape=u.pshape, scale=u.pscale);
  Kv <- 1 # rgamma(n=1, shape=v.pshape, scale=v.pscale);
  K <- c(Ku,Kv);
  return(K);
}

samplex <- function(dim,Q,K,y) {
  center <- c(rep(0,dim),y);
  R <- makeR(dim,Q,K);
  # cholR<-chol(R)
  x <- rmvnorm.canonical(n=1, b=center, Q=R, Rstruct=cholR); # exploit the stability of the structure of R, at least 10 times faster.
  # cholR has to be computed at the initiation with the same structure
  return(drop(x));
}

makeR <- function(dim,Q,K) {
  Ku <- K[1];
  Kv <- K[2];
  R <- Ku*Q;
  diag.spam(R) <- diag.spam(R) + Kv;
  R <- cbind(R, diag.spam(-1*Kv,dim,dim));
  R <- rbind(R, cbind(diag.spam(-1*Kv,dim,dim), diag.spam((Kv+1), dim, dim)))
  return(R);
}

library(msm) # for the truncated normal distribution
sampley <- function(dim,ut,yprime) {
  center <- ut[1:dim] + ut[(dim+1)];
  lwbd <- rep(0,dim);
  lwbd[(yprime==0)] <- -Inf;
  upbd <- rep(0,dim);
  upbd[(yprime==1)] <- Inf;
  y <- rtnorm(n=dim, mean=center, sd=1, lower=lwbd, upper=upbd);
  return(y);
}

sampleyprime <- function(dim,ut,betaprime,zpos,zneg,zNA) {
  u <- ut[1:dim];
  t <- ut[(dim+1)];
  yprime <- rep(0,dim);
  yprime[zpos] <- 1;
  p <- (1-betaprime[zneg])*(pnorm(q=0, mean=(u[zneg]+t), sd=1));
  q <- (1 - pnorm(q=0, mean=(ut[zneg]+t), sd=1));
  Pzneg <- p/(p + q);
  PzNA <- pnorm(q=0, mean=(ut[zNA]+t), sd=1);
  yprime[zneg] <- rbinom(length(zneg),1,Pzneg);
  yprime[zNA] <- rbinom(length(zNA),1,PzNA);
  return(yprime);
}

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

# likelihood of z | u,b
llhzgivenbetax <- function(dim,zpos,zneg,ut,betaprime){
  x<-ut[1:dim]+ut[1+dim];
  llhpos<-sum(log((1-pnorm(0,x[zpos],1))*betaprime[zpos]));
  llhneg<-sum(log((1-pnorm(0,x[zneg],1))*(1-betaprime[zneg])+ pnorm(0,x[zneg],1)));
  llh<-llhpos+llhneg;
  return(llh);
}

### Initialization
# Covariance matrix construction
spam.options(nearestdistnnz=c(18104083,400))
dist_mat <-nearest.dist(x=data[,c("easting","northing")], y=NULL, method="euclidian", delta=threshold, upper=NULL);          

y<-rep(inity,dimension);
K <- c(Ku,Kv);
Q<-QfromDelta(Delta,dist_mat,AS,f=f);

R <- makeR(dimension,Q,K);
diag.spam(R) <- diag.spam(R) + 0.000001;
cholR <- chol(R, memory=list(nnzcolindices=5e6));
x <- drop(rmvnorm.canonical(n=1, b=c(rep(0,dimension),y), Q=R, Rstruct=cholR));

LLH <- 	       LLHDeltaQx(Delta     ,y.r,x,     Q,    K,muDelta,sdDelta);
# x<-c(rep(0,length(y)),y)
LLH <- 	       LLHDeltaQx(Delta     ,y.r,x,     Q,    K,muDelta,sdDelta);

n=min(freqsave,nbsimul);

K<-c(K,Delta,LLH);
write.table(t(K), "Ksamples.txt", sep="\t",col.names=FALSE,row.names=FALSE)

Ksamples <- matrix(0,n,4);
fin <- c();

j=0;
graphics.off()
name <- paste("u.r for mu.r",mu.r,", Ku.r",Ku.r,", tr",threshold,", f",f.r,", Delta",Delta.r)
plot_reel(data$easting,data$northing,u.r,main=name)
dev.new()
# par(mfrow=c(1,2))
scan()
for (i in 1:nbsimul) {
	j=j+1;

	Delta_prop=rtnorm(1,mean=Delta,sd=sdSampleDelta,lower=0,upper=threshold);
	Qprop<-QfromDelta(Delta_prop,dist_mat,AS);
	Q<-QfromDelta(Delta,dist_mat,AS,f=f);

	x_prop <- samplex(dimension,Qprop,K,y.r); 

	hasting_term=log(dtnorm(Delta,mean=Delta_prop,sd=sdSampleDelta,lower=0,upper=threshold)/dtnorm(Delta_prop,mean=Delta,sd=sdSampleDelta,lower=0,upper=threshold));

	LLHproposal <- LLHDeltaQx(Delta_prop,y.r,x_prop,Qprop,K,muDelta,sdDelta);
	LLH <- 	       LLHDeltaQx(Delta     ,y.r,x,     Q,    K,muDelta,sdDelta);
	lnr <- LLHproposal-LLH+hasting_term;
	cat("Delta:",Delta,"LLH:",LLH,"Delta prop:",Delta_prop,"LLHproposal:",LLHproposal,"h_t",hasting_term,"lnr",lnr);

	if(lnr>=log(runif(1))) {
		Delta <- Delta_prop;
		x<-x_prop
		acceptDelta<<-c(acceptDelta,1);
		cat(" accept 1\n");
	}else{
		acceptDelta<<-c(acceptDelta,0);
		cat("accept 0\n");
	}
	plot_reel(data$easting,data$northing,x[1:dimension],main=name)

	scan()

  print(i);
}
# possibly in another R session :
Ksamples<-read.table("Ksamples.txt")
names(Ksamples)<-c("Ku","Kv","Delta")
dev.new()
par(mfrow=c(1,3))
plot(Ksamples$Ku)
plot(Ksamples$Kv)
plot(Ksamples$Delta)
dev.print(device=pdf,"Ku_Kv_Delta.pdf")
dev.new()
plot(Ksamples$Ku,Ksamples$Kv)
dev.print(device=pdf,"Ku_Kv.pdf")



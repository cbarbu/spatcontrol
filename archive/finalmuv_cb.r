### Entering data into R
set.seed(500)
graphics.off()

nbsimul <- 1000;
save <- 10;
write <- 10;
avg <- 100;
grid <- 100*(1:120);
burn <- 5e4;
name<-"finaluv_noNA"

# setwd("C:/Users/ahon/Desktop/Project/")

# data <- read.csv("data.csv",header=T);
# colnames(data) <- c("locality","status","collector","northing","easting");
data.base <- read.csv("DB_mm_blocks_05May2011.csv",header=TRUE);
data.base<-data.base[data.base$status!=9,];
# restriction to high inspectors' efficiency zone
# data<-data.base[data.base$northing>8185000&data.base$easting<233400,]
# data<-data.base[data.base$northing<8185500&data.base$easting>233100,] # bugs with bad inspectors 
data<-data.base[data.base$northing>8185000&data.base$northing<8185500&data.base$easting>233250&data.base$easting<233750,]

### cleaning of inspectors
data$collector[data$collector==" Jorge Ampuero"] <- "Jorge Ampuero"
data$collector[data$collector==" Jorge A"] <- "Jorge Ampuero"
data$collector[data$collector=="Jose Velasquez "] <- "Jose Velasquez"
data$collector[data$collector==" Hugo Vilcahuaman"] <- "Hugo Vilcahuaman"
data$collector[data$collector=="Julio cesar Condori"] <- "Julio Cesar Condori"
data$collector[data$collector=="Manuel Tamayo "] <- "Manuel Tamayo"

data$collector<-factor(data$collector)

### Packages for this section
library(spam);

spam.options(cholsymmetrycheck=FALSE, safemode=c(FALSE,FALSE,FALSE)) ## to speed up things once everything is ok, check same results
source("spam_complement.r")

# remove houses alone in their block (including houses alone because other parts of the block are out of the map)
 SB <- nearest.dist(x=cbind(data$block_num,rep(0,length(data$block_num))), method="euclidian", upper=NULL,delta=0.1)
SB@entries<-rep(1,length(SB@entries))
nbvsb<-SB%*%rep(1,length(data$block_num))-1
isolated<-which(nbvsb==0);

if(length(isolated)>0){
	data<-data[-isolated,];
	cat("remove:",isolated,"\n");
}

### Data cleaning
data <- data[order(data$locality, data$easting, data$northing),];
dimension <- nrow(data);

zpos <- which(data$status==1);
zneg <- which(data$status==0);
zNA <- which(data$status==9);
select<-c(zneg,zpos)

data.insp <- data$collector;
inspectors <- unique(data.insp[c(zpos,zneg)]);
inspector <- matrix(0,dimension,length(inspectors))
for (i in 1:length(inspectors)) {
  inspector[,i] <- (data.insp==inspectors[i]);
}
inspector <- as.spam(inspector);

### Covariance matrix construction
dist_mat <-nearest.dist(x=data[,c("easting","northing")], y=NULL, method="euclidian", delta=threshold, upper=NULL);          
dist_mat<-as.spam(dist_mat)
source("DeltaSampling.r")
source("sample_y_direct.r")

threshold <- 100;

spam.options(nearestdistnnz=c(9058076,400))
Q <-nearest.dist(x=data[,c("easting","northing")], y=NULL, method="euclidian", delta=threshold, upper=NULL);          
diag.spam(Q) <- 1;
Q <- -1*(1/Q);
diag.spam(Q) <- 0;
rsumQ <- Q %*% rep(1,dimension);
diag.spam(Q) <- -1*rsumQ;

par(mfrow=c(1,2))
# f with log normal prior
mf<-1# 0.05;
sdlf<-1 # 2
lik.f<-function(f,mf,sdlf){
	LLH<-dlnorm(f,log(mf),sdlf,log=TRUE)
	return(LLH);
}
xabs<-seq(-10,2,0.1)
plot(xabs,lik.f(10^xabs,mf,sdlf),main="f prior",xlab="f",ylab="density")
logsdfprop<-0.1

# T
mT<-1;
sdlT<-1;
lik.T<-lik.f
xabs<-seq(-10,2,0.1)
plot(xabs,lik.T(10^xabs,mT,sdlT),main="T prior",xlab="T (log10)",ylab="density")
logsdTprop<-0.3

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

makeR <- function(dim,Q,K) {
  Ku <- K[1];
  Kv <- K[2];
  R <- Ku*Q;
  diag.spam(R) <- diag.spam(R) + Kv;
  R <- cbind(R, diag.spam(-1*Kv,dim,dim));
  R <- rbind(R, cbind(diag.spam(-1*Kv,dim,dim), diag.spam((Kv+1), dim, dim)));
  return(R);
}

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

sampley <- function(dim,x,yprime) {
  center <- x[(dim+(1:dim))];
  lwbd <- rep(0,dim);
  lwbd[(yprime==0)] <- -Inf;
  upbd <- rep(0,dim);
  upbd[(yprime==1)] <- Inf;
  y <- rtnorm(n=dim, mean=center, sd=1, lower=lwbd, upper=upbd);
  return(y);
}

sampleyprime <- function(dim,x,betaprime,zpos,zneg,zNA) {
  w <- x[(dim+(1:dim))];
  yprime <- rep(0,dim);
  yprime[zpos] <- 1;
  p <- (1-betaprime[zneg])*(pnorm(q=0, mean=(w[zneg]), sd=1));
  q <- (1 - pnorm(q=0, mean=(w[zneg]), sd=1));
  Pzneg <- p/(p + q);
  PzNA <- pnorm(q=0, mean=(w[zNA]), sd=1);
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

### Initialization

library(msm)

K.hyper <- rep(0.1,4); 
a <- 100; 
b <- 1;
f<-1
Delta<-0
T<-1

K <- c(rgamma(1,K.hyper[1],K.hyper[2]), rgamma(1,K.hyper[3],K.hyper[4]));
R <- K[1]*Q;
diag.spam(R) <- diag.spam(R) + K[1]*(1e-5);
cholR <- chol(R);
u <- rnorm(dimension,0,1);
u <- backsolve(cholR,u);
w <- rnorm(dimension,u,sqrt((K[2])^(-1)));
x <- c(u, w);
y <- rnorm(dimension,x[(dimension+(1:dimension))],1); # don't we work on real data now ?
yprime <- (y>0);
beta <- rbeta(ncol(inspector),a,b);
betaprime <- as.vector(inspector %*% beta);
rm(u,w);
cholR<-NULL;

Ksamples <- matrix(0,write,2);
usamples <- matrix(0,write,length(grid));
wsamples <- matrix(0,write,length(grid));
ysamples <- matrix(0,write,length(grid));
bsamples <- matrix(0,write,ncol(inspector));

ypsum <- rep(0,dimension);

write.table(t(K), "Ksamples.txt", sep="\t", col.names=F, row.names=F);
write.table(t(x[grid]), "usamples.txt", sep="\t", col.names=F, row.names=F);
write.table(t(x[(dimension+(grid))]), "wsamples.txt", sep="\t", col.names=F, row.names=F);
write.table(t(y[grid]), "ysamples.txt", sep="\t", col.names=F, row.names=F);
write.table(t(beta), "bsamples.txt", sep="\t",col.names=FALSE,row.names=FALSE);

j=0;
scan()

par(mfcol=c(3,5))
for (i in 1:nbsimul) {
	cat(" i:",i);
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


	K <- sampleK(dimension,Q,x,K.hyper);
	x <- samplex(dimension,Q,K,y,cholR=cholR);
	y <- sampley(dimension,x,yprime);
	yprime <- sampleyprime(dimension,x,betaprime,zpos,zneg,zNA);
	# y<-sample_y_direct(x[1:dimension],zpos,zneg,zNA,betaprime);
	# yprime<-as.numeric(y>0);
	beta <- samplebeta(zpos,zneg,inspector,yprime,a,b);
	# beta<-rep(1,inspector@dimension[2])
	betaprime <- as.vector(inspector %*% beta);
	if((i%%save)==0){
		j = j+1;
		Ksamples[j,] <- K;
		usamples[j,] <- x[grid];
		wsamples[j,] <- x[(dimension+grid)];
		ysamples[j,] <- y[grid];
		bsamples[j,] <- beta;

		plot_reel(data$easting,data$northing,x[1:dimension],main="x")
		plot_reel(data$easting,data$northing,x[(1+dimension):(2*dimension)],main="y")
		plot_reel(data$easting[select],data$northing[select],(betaprime[select]*2)-1,top=1,base=0,main="inspectors")
	}
	if(j==write){
		write.table(Ksamples, "Ksamples.txt", sep="\t", append=T, col.names=F, row.names=F);
		write.table(usamples, "usamples.txt", sep="\t", col.names=F, row.names=F);
		write.table(wsamples, "wsamples.txt", sep="\t", append=T, col.names=F, row.names=F);
		write.table(ysamples, "ysamples.txt", sep="\t", append=T, col.names=F, row.names=F);
		write.table(bsamples, "betasamples.txt", sep="\t", append=T, col.names=F, row.names=F);
		j=0;
	}
	if(((i%%avg)==0) && (i>=burn)) {
		ypsum <- ypsum + yprime;
	}
}

dev.new()
par(mfrow=c(1,2))
K<-read.table("Ksamples.txt")
plot(K[,1],type="l")
plot(K[,2],type="l")




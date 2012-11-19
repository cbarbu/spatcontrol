### Entering data into R

# setwd("C:/Users/ahon/Desktop/Project/")

# data <- read.csv("data.csv",header=T);
# colnames(data) <- c("locality","status","collector","northing","easting");
data.base <- read.csv("DB_mm_blocks_05May2011.csv",header=TRUE);
data.base <-data.base[data.base$status!=9,]

# restriction to high inspectors' efficiency zone
data<-data.base[data.base$northing>8185000&data.base$easting<233400,] # bugs with good inspectors
# data<-data.base[data.base$northing<8185500&data.base$easting>233100,] # bugs with bad inspectors 
# data<-data.base[data.base$northing>8185000&data.base$northing<8185500&data.base$easting>233250&data.base$easting<233750,]

# data<-data.base;

### cleaning of inspectors
data$collector[data$collector==" Jorge Ampuero"] <- "Jorge Ampuero"
data$collector[data$collector==" Jorge A"] <- "Jorge Ampuero"
data$collector[data$collector=="Jose Velasquez "] <- "Jose Velasquez"
data$collector[data$collector==" Hugo Vilcahuaman"] <- "Hugo Vilcahuaman"
data$collector[data$collector=="Julio cesar Condori"] <- "Julio Cesar Condori"
data$collector[data$collector=="Manuel Tamayo "] <- "Manuel Tamayo"

data$collector<-factor(data$collector)

### Packages for this section

library(msm);

library(spam);
spam.options(cholsymmetrycheck=FALSE, safemode=c(FALSE,FALSE,FALSE)) ## to speed up things once everything is ok, check same results

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

data.insp <- data$collector;
inspectors <- unique(data.insp[c(zpos,zneg)]);
inspector <- matrix(0,dimension,length(inspectors))
for (i in 1:length(inspectors)) {
  inspector[,i] <- (data.insp==inspectors[i]);
}
inspector <- as.spam(inspector);

### Covariance matrix construction

threshold <- 100;

spam.options(nearestdistnnz=c(9058076,400))
Q <-nearest.dist(x=data[,c("easting","northing")], y=NULL, method="euclidian", delta=threshold, upper=NULL);          
diag.spam(Q) <- 1;
# Q <- -1*(1/Q);
Q <- -exp(-0.05*Q);
diag.spam(Q) <- 0;
rsumQ <- Q %*% rep(1,dimension);
diag.spam(Q) <- -1*rsumQ;

sampleku <- function(dim,Q,u,K.hyper) {
  Ku.a <- K.hyper[1];
  Ku.b <- K.hyper[2];
  pos.shape <- (0.5*(dim-1) + Ku.a);
  pos.scale <- (0.5*as.numeric(u %*% Q %*% u) + Ku.b^(-1))^(-1);
  ku <- rgamma(n=1, shape=pos.shape, scale=pos.scale);
  return(ku);
}

makeR <- function(dim,Q,ku) {
  R <- ku*Q;
  diag.spam(R) <- diag.spam(R) + 1;
  return(R);
}

sampleu <- function(dim,Q,ku,y,cholR=NULL) {
  u <- rnorm(n=dim, mean=0, sd=1);
  center <- y;
  R <- makeR(dim,Q,ku);
  if(is.null(cholR)){
	  cholR <- chol.spam(R, memory=list(nnzcolindices=4e6));
  }else{
	  cholR <- update.spam.chol.NgPeyton(cholR,R);
  }
  center <- backsolve(cholR, forwardsolve(cholR, center));
  u <- backsolve(cholR,u);
  u <- u + center;
  cholR <<- cholR;
  return(u);
}

sampley <- function(dim,u,yprime) {
  center <- u;
  lwbd <- rep(0,dim);
  lwbd[(yprime==0)] <- -Inf;
  upbd <- rep(0,dim);
  upbd[(yprime==1)] <- Inf;
  y <- rtnorm(n=dim, mean=center, sd=1, lower=lwbd, upper=upbd);
  return(y);
}

sampleyprime <- function(dim,u,betaprime,zpos,zneg,zNA) {
  yprime <- rep(0,dim);
  yprime[zpos] <- 1;
  p <- (1-betaprime[zneg])*(pnorm(q=0, mean=(u[zneg]), sd=1));
  q <- (1 - pnorm(q=0, mean=(u[zneg]), sd=1));
  Pzneg <- p/(p + q);
  PzNA <- pnorm(q=0, mean=(u[zNA]), sd=1);
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

K.hyper <- rep(0.1,2);
a <- 10;
b <- 1;

ku <- rgamma(1,K.hyper[1],K.hyper[2]);
R <- ku*Q;
diag.spam(R) <- diag.spam(R) + ku*(1e-5);
cholR <- chol(R);
u <- rnorm(dimension,0,1);
u <- backsolve(cholR,u);
y <- rnorm(dimension,u,1);
yprime <- (y>0);
beta <- rbeta(ncol(inspector),a,b);
# beta <- rep(1,ncol(inspector));

### Sampler

n <- 1000;
save <- 10;
write <- 10;
avg <- 100;
grid <- 100*(1:120);
burn <- 5e4;
cholR=NULL;

kusamples <- matrix(0,write,1);
usamples <- matrix(0,write,length(grid));
ysamples <- matrix(0,write,length(grid));
bsamples <- matrix(0,write,ncol(inspector));

ypsum <- rep(0,dimension);

write.table(t(ku), "kusamples.txt", sep="\t", col.names=F, row.names=F);
write.table(t(u[grid]), "usamples.txt", sep="\t", col.names=F, row.names=F);
write.table(t(y[grid]), "ysamples.txt", sep="\t", col.names=F, row.names=F);
write.table(t(beta), "bsamples.txt", sep="\t",col.names=FALSE,row.names=FALSE);
write.table(t(u), "usamplesfull.txt", sep="\t", col.names=F, row.names=F);

j=0;
k=0;

for (i in 1:n) {
	cat("i:",i);
	ku <- sampleku(dimension,Q,u,K.hyper);
	cat("ku:",ku);
	u <- sampleu(dimension,Q,ku,y,cholR);
	cat("sdu",sd(u));
	y <- sampley(dimension,u,yprime);
	betaprime <- as.vector(inspector %*% beta);
	yprime <- sampleyprime(dimension,u,betaprime,zpos,zneg,zNA);
	beta <- samplebeta(zpos,zneg,inspector,yprime,a,b);
	if((i%%save)==0){
		j = j+1;
		kusamples[j,] <- ku;
		usamples[j,] <- u[grid];
		ysamples[j,] <- y[grid];
		bsamples[j,] <- beta;
		cat("j:",j);
	}
	if(j==write){
		k=k+1;
		write.table(kusamples, "kusamples.txt", sep="\t", append=T, col.names=F, row.names=F);
		write.table(usamples, "usamples.txt", sep="\t", append=T,col.names=F, row.names=F);
		write.table(t(u), "usamplesfull.txt", sep="\t", append=T,col.names=F, row.names=F);
		write.table(ysamples, "ysamples.txt", sep="\t", append=T, col.names=F, row.names=F);
		write.table(bsamples, "bsamples.txt", sep="\t", append=T, col.names=F, row.names=F);
		j=0;
		cat("k:",k);
	}
	if(((i%%avg)==0) && (i>=burn)) {
		ypsum <- ypsum + yprime;
	}
	cat("\n");
}

source("spam_complement.r")
 par(mfrow=c(2,2))
plot_reel(data$easting,data$northing,data$status,main="Data")
plot_reel(data$easting,data$northing,u,main="u")
plot_reel(data$easting,data$northing,y,main="y")
plot_reel(data$easting,data$northing,betaprime,main="Insp")




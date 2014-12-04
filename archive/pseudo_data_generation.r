
source("functions_intercept.r")
if(!exists("GLOBALSETPARAMETERS")){
	source("parameters_sampler.r")
}else{
	if(GLOBALSETPARAMETERS==TRUE){
		source("parameters_sampler.r")
	}
}
source("inspectors_detection.r")
## data
source("spam_complement.r")
source("import_data.r")

## functions
spam.options(nearestdistnnz=c(9058076,400))
dist_mat <-nearest.dist(x=db[,c("easting","northing")], y=NULL, method="euclidian", delta=threshold, upper=NULL);          
# dist_mat<-as.spam(dist_mat)
## Nota:
# diag.spam() is quite slow we then save the position of the diagonal terms in the 
# spam entries of Dmat:
diag(dist_mat)<- rep(0,dim(dist_mat)[1])
# diag(dist_mat)<- rep(-1,dim(dist_mat)[1])
# diagDmat<-which(dist_mat@entries==-1) 
# dist_mat@entries[diagDmat]<-rep(0,length(diagDmat))
# it allows to do Dmat@entries[diagDmat]<-rep(1,length(diagDmat))
# in spite of diag.spam(Dmat)<-1
# 6 times faster, can be used also for the diagonal of Q

# ## remove isolated
# isolated<-which(apply_by_row_not_null.spam(dist_mat,min)>threshold);
# if(length(isolated)>0){
# 	db<-db[-isolated,]
# 	dist_mat <-nearest.dist(x=db[,c("easting","northing")], y=NULL, method="euclidian", delta=threshold, upper=NULL);          
# }
# dist_mat@entries[dist_mat@entries==0]<-1 # add 1m of distances between A/B houses

# dist_mat<-as.spam(dist_mat)
dimension <- nrow(db);

source("DeltaSampling.r")
library(msm) # for the truncated normal distribution

spam.options(nearestdistnnz=c(13764100,400))
SB <- nearest.dist(x=cbind(db$block_num,rep(0,length(db$block_num))), method="euclidian", upper=NULL,delta=0.1)
SB@entries<-rep(1,length(SB@entries))
dmt<-dist_mat
dmt@entries<-rep(1,length(dmt@entries))# [dmt@entries!=0]<-1 # 1 only when dist_mat not 0
SB@entries<-rep(1,length(SB@entries))
SB<-as.spam(SB*dmt);

AS<-dmt-SB; # get 1 whereever the distances matrix is defined(under threshold) and not same block
AS<-as.spam(AS)

K.r<-c(Ku.r,Kv.r)

Dmat<-dist_mat
if(use.streets){
	Q.r<-QfromfT(Dmat,AS,SB,f=f.r,T=T.r);
}else{
	Q.r<-QfromfT(Dmat,AS,SB,f=f.r,T=1);
}

## generation of cofactors
gen_c.map<-function(nbfact,nbpoints,rates=rep(0.5,nbfact)){
	c.map={}
	for(fact in 1:nbfact){
		prob<-c(1-rates[fact],rates[fact])
		c.col<-sample(0:1,nbpoints,replace=TRUE,prob=prob)
		c.map<-cbind(c.map,c.col);
	}
	return(c.map);
}
# gen_c.map(3,10,rates=c(0.1,0.9,0.5))
if(use.cofactors){
	if(make.map.cofactors){
		c.map<-gen_c.map(nbfact.gen,nrow(Q.r),rep(prob.fact,nbfact.gen));

		# c.map<-cbind(rep(1,nrow(c.map)),c.map); # add a fake parameter always present~v
		c.map.plus<-pseudo_inv(c.map)
		c.map<-as.spam(c.map);
		# Sc<-c.map%*%t(c.map)
		# NB: the rank is given by qr(Qc)$rank and if the rank is < nrow(Qc) the det is null
		# as it is the case here
		# Qc<-as.spam(solve.spam(Sc+diag.spam(epsilon,nrow(Sc)))) # to allow the generation of data
		# rm(Sc)

		## compare generation by c.val and generation by c.map%*%c.val
		# generation by c.val
	}else{
		c.map<-as.matrix(db[,cofs])
	}
	c.val1<-rnorm(nbfact.gen,0,sd=sqrt(1/Kc.r)); 
c.comp1<-c.map%*%c.val1
# generation by c.map%*%c.val
# c.comp2<-drop(rmvnorm.spam(n=1,mu=rep(0,dimension),Sigma=(Sc.r/(Kc.r))));
# c.comp3<-drop(rmvnorm.prec(n=1,mu=rep(0,dimension),Q=(Qc*Kc.r)));
# c.val2<-correct*c.map.plus%*%c.comp2
# check similarity c.comp1/c.comp2
graphics.off()
#  par(mfrow=c(2,3))
#  plot(c.comp3)
#  plot(c.map.plus%*%c.comp3)
#  plot(c.map%*%c.map.plus%*%c.comp3)
#  plot(c.comp2)
#  plot(c.val2)
#  plot(c.map%*%c.val2)
#  dev.new()
#  cat("sd val (",sd(c.val1),",",sd(c.val2),") sd comp:",sd(c.comp1),",",sd(c.comp2),",",sd(c.map%*%c.val2),"\n");
# plot(c.comp2,c.map%*%c.val2)
c.val.r<-c.val1
c.r<-c.comp1;
# c.val.r<-c.map.plus%*%c.r
cat("c values:",c.val.r,"\n");
#  ## plot c.map
#  par(mfrow=c(2,2))
#  for(comp in 1:nbfact.gen){
#  	select<-rep(0,nbfact.gen);
#  	select[comp]<-c.val.r[comp]
#  	nameplot<-paste("c[",comp,"]",sep="");
#  	plot_reel(db$X,db$Y,c.map%*%select,main=nameplot)
#  }
#  plot_reel(db$X,db$Y,c.r,main="c.r")
}
dev.new()

# fixing x and y
## this systematically anulate the spatial effect of u.r by v.r to get y at the mean
# R.r <- makeR(dim=dimension, Q=Q.r, K=c(Ku.r,Kv.r));
# cholR.r <- chol.spam(R.r, memory=list(nnzcolindices=5e6));
# x.r <- rmvnorm.canonical(n=1, b=(c(rep(0,dimension),rep(m,dimension))), Q=R.r, Rstruct=cholR.r);
# u.r<-x.r[1:dimension]
# v.r<-x.r[(dimension+1):(2*dimension)]-x.r[1:dimension]
# y <- rnorm(n=dimension, mean=x.r[c(dimension+(1:dimension))], sd=1);

## if want to change given_rnorm
# given_rnorm<-2*db$status-1 # rnorm(nrow(Q.r))
given_rnorm<- rnorm(nrow(Q.r))
# dump("given_rnorm",file="given_rnorm.r")
# source("given_rnorm.r")
# u.r <- rmvnorm.canonical.pseudo(n=1, b=rep(mu.r,dimension), Q=Ku.r*Q,given_rnorm=given_rnorm);
# Qu.samp<-Q.r+diag.spam(epsilon,dimension);
cholQ.r<-chol(Q.r)
## generate spatially random
u.r <- drop(rmvnorm.prec.pseudo(n=1, mu=rep(mu.r,dimension), Q=Ku.r*Q.r,given_rnorm=given_rnorm));
u.r<-u.r-mean(u.r)+mu.r
## generate arround real distribution
# u.r <- drop(rmvnorm.canonical(n=1,1*(db$status-0.5), Q=Ku.r*Q.r+diag.spam(Kv.r,dimension),given_rnorm=given_rnorm));
# u.r <- drop(rmvnorm.prec.pseudo(n=1, mu=rep(mu.r,dimension), Q=Ku.r*Q.r,given_rnorm=given_rnorm));

Qv<-diag.spam(Kv.r,nrow(Q.r))
w.r<-rep(0,dimension);
if(use.cofactors){
	w.r<-w.r+c.r;
}
if(use.v.gen){
	v.r <- drop(rmvnorm.prec.pseudo(n=1, mu=rep(mv.r,dimension), Q=Qv));
	w.r<-w.r+u.r+v.r
}else{
	w.r<-w.r+u.r
}
y.r <- rnorm(n=dimension, mean=w.r, sd=1);
x.r<-c(u.r,w.r)

db.insp <- db$collector;

inspectors <- unique(db.insp[which(db$status<9)]);
inspector <- matrix(0,dimension,length(inspectors))
for (i in 1:length(inspectors)) {
  inspector[,i] <- (db.insp==inspectors[i]);
}
inspector <- as.spam(inspector);

# beta.r<-merge(data.frame(inspectors),inspectors_detection,by.x="inspectors",by.y="inspectors")[[2]];
# beta.r<-rbeta(length(inspectors),1,1)
b.r<-rep(beta.r,length(inspectors))
bivect <- as.vector(inspector %*% b.r);
# par(mfrow=c(1,1))

par(mfrow=c(2,3))
nameplot <- paste("u.r for mu.r",mu.r,", Ku.r",Ku.r,", tr",threshold,", f",f.r,", T",T.r)
plot_reel(db$X,db$Y,u.r,main=nameplot)
print(summary(u.r))
if(use.cofactors){
name_plot <- paste("c.r for Kc.r",Kc.r)
plot_reel(db$X,db$Y,c.r,main=name_plot)
}
if(use.v.gen){
name_plot <- paste("v.r for Ku.r",Ku.r,"Kv.r",Kv.r,"mv.r",mv.r)
plot_reel(db$X,db$Y,v.r,main=name_plot)
}
nameplot <- paste("y.r for Ku.r",Ku.r,"Kv.r",Kv.r,"mu.r",mu.r,"mv.r",mv.r)
plot_reel(db$X,db$Y,y.r,main=nameplot)

plot_reel(db$X,db$Y,2*(bivect-0.5),main="bivect")
z.r<-generate_z(y.r,bivect,which(db$status==9)); 

nameplot <- paste("z.r for Ku.r",Ku.r,"Kv.r",Kv.r,"mu.r",mu.r,"mv.r",mv.r)
plot_reel(db$X,db$Y,z.r*2-1,main=nameplot)
dump("z.r",file="generated_z.r");

# print("LLH reference:")
# print(LLHDeltaQx(Delta.r,y.r,x.r,Q.r,K.r,muDelta,sdDelta))[[1]]
# LLHDeltaQx(Delta.r,y.r,rep(0,2*dimension),Q.r,K.r,muDelta,sdDelta)
# R.r<-makeR(dimension,Q.r,K);
# R <- Q.r;
# diag.spam(R) <- diag.spam(R) + 1;
# R <- cbind.spam(R, diag.spam(-1,dimension));
# R <- rbind.spam(R, cbind.spam(diag.spam(-1,dimension), diag.spam((2), dimension)))
# det(R,memory=list(nnzcolindices=10e7))
# det(R)
# R <- Q.r;
# R <- bdiag.spam(R,diag.spam((2), 1))
# det(R)

# print(summary(v.r))
# print(summary(y))

# # to visualize problems with AS
# varQ<- -Q
# diag(varQ)<-0
# u_init<-rep(0,length(db$status))
# u_init[sample(length(db$status),30)]<-1
# 
# u<- varQ%*%u_init
# par(mfrow=c(1,2))
# plot_reel(db$X,db$Y,u_init,main="u_init")
# plot_reel(db$X,db$Y,u,main="u")
# 
# dev.new()
# plot(db$X,db$Y,col=db$block_num,asp=1)
# text(db$X[which(u_init==1)],db$Y[which(u_init==1)],label=(1:length(u_init))[which(u_init==1)])
# 
# dev.new()
# selection<-(1070:1080)
# par(mfrow=c(1,3))
# plot_reel(db$X[selection],db$Y[selection],u_init[selection],main=nameplot)
# plot_reel(db$X[selection],db$Y[selection],u[selection],main=nameplot)
# plot(db$X[selection],db$Y[selection],col=db[selection,]$block_num,cex=3,pch=16,asp=1)
# text(db$X[selection],db$Y[selection],label=(1:length(u_init))[selection])
## corresponding yprime
# yprime <- y>0
# 
# plot(db$X,db$Y,col=1,pch=15,asp=1,main="yprime")
# lines(db$X[yprime],db$Y[yprime],col=7,type="p",pch=15,cex=0.6)

# ## fixing beta
# # need where are the inspectors
# zNA <- which(db$status==9);
# znoNA <-which(db$status==1 | db$status==0)
# db.insp <- db$collector;
# inspectors <- unique(db.insp[znoNA]); 
# inspector <- matrix(0,dimension,length(inspectors))# matrix houses on rows and 1 in the column corresponding to the inspector, only 0 if no inspector
# for (i in 1:length(inspectors)) {
#   inspector[,i] <- (db.insp==inspectors[i]);
# }
# inspector <- as.spam(inspector);
# 
# # sampling of beta itself
# a <- 35;
# b <- 4;
# beta <- rbeta(ncol(inspector),a,b); # detection rate per inspector
# betaprime <- as.vector(inspector %*% beta); # detection per house
# nameplot <- paste("beta for a",a,"b",b)
# # plot_reel(db$X,db$Y,betaprime,main=nameplot,base=0.1,top=0.9)
# 
# 
# # sampling the db: z 
# probzinspect <- yprime*betaprime; # given the structure of betaprime the NA have probability 0
# z<-rbinom(length(db[,1]),1,probzinspect);   
# zpos <- which(z==1)
# zneg <- which(z==0)
# # dev.new()
# plot(db$X,db$Y,col=1,pch=15,asp=1,main="z")
# lines(db$X[zNA],db$Y[zNA],col=4,type="p",pch=15,cex=0.2)
# lines(db$X[zpos],db$Y[zpos],col=7,type="p",pch=15,cex=0.6)
# 
# dev.new()
# # scan()
# par(mfrow=c(2,3))
# # for(Ku.r in c(0.005,0.5,5,500)){
# # 	for(Kv.r in c(0.005,0.5,5,500)){
# for(m in seq(-2,2,1)){
# 		R.r <- makeR(dim=dimension, Q=Q, K=c(Ku.r,Kv.r));
# 		x.r <- rmvnorm.canonical(n=1, b=(c(rep(m,dimension),rep(m,dimension))), Q=R.r, Rstruct=cholR.r);
# 		y <- rnorm(n=dimension, mean=x.r[c(dimension+(1:dimension))], sd=1);
# 
# 		nameplot <- paste("Ku.r",Ku.r,"Kv.r",Kv.r,"m",m)
# 		print(nameplot)
# 		plot_reel(db$X,db$Y,y,main=nameplot)
#   	}
# #   }

# # autocorrelation of the pseudo-data
# houses_XY2<- db[c("easting","northing")]; # not unicode and number of bugs
# pres_abs <- y # get presence absence data 
# library("spdep");
# temp=system.time(dnb <- dnearneigh(as.matrix(houses_XY2), 0, 50));
# temp=system.time(lw <- nb2listw(dnb, zero.policy=TRUE));
# temp=system.time(mt <- moran.test(pres_abs, lw, zero.policy=TRUE,adjust.n=FALSE));
# print(mt)

# # sensibility of the precision P of u+v to changes of Ku when accomodated by Kv
# Kur<-1
# Kvr<-1
# V<-1/(1/Kvr + 1/(Kur*diag(Q)))
# plot(V,type="l")
# 
# M<-mean(diag(Q))
# 
# 	color=2
# for( Ku in c(0.5,1,10,100)){
# 	Kv<-max(0.001,1/(1/Kvr+1/(Kur*M)-1/(Ku*M)))
# 	V<-1/(1/Kv + 1/(Ku*diag(Q)))
# 	cat("Ku",Ku,"Kv",Kv,"mean(V)",mean(V))
# 	lines(V,col=color)
# 	color<-color+1
# 	scan()
# }
# # Ku is only usefull for the variance arround the mean of the variance of the us. A variance of the variance of the u.
# # Honestly, who cares? Can be great to accomodate for any kind of strange patterns but useless otherwise

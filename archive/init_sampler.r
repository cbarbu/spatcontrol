## general
if(use.insp){
	data.insp <- db$collector;
	inspectors <- unique(data.insp[which(db$status!=9)]);
	inspector <- matrix(0,dimension,length(inspectors))
	for (i in 1:length(inspectors)) {
		inspector[,i] <- (data.insp==inspectors[i]);
	}
	inspector <- as.spam(inspector);
	beta<-rep(1,inspector@dimension[2])
	bivect <- as.vector(inspector %*% beta);
}else{
	bivect<-rep(priorinspquality,dimension)
}


# starting values
K<-c(Ku,Kv,Kc);
if(use.streets){
	Q<-QfromfT(dist_mat,AS,SB,f,T);
}else{
	Q<-QfromfT(dist_mat,AS,SB,f,T=1);
}

cat("at Q build\n T:",T,"f:",f,"Ku:",Ku,"Kv:",Kv,"\n")
cholQ<-chol(Q);

u<-rep(0,dimension);
y<-rep(0,dimension);
yprime <- (y>0);
if(use.v){
w <- rnorm(dimension,u,sqrt((K[2])^(-1)));
v<-w-u;
x <- c(u, w);
}else if(use.intercept){
	x<-c(u,int)
	w<-u+int
}
if(use.cofactors){
	c.val<-rep(0,length(cofs));
	c.comp<-c.map%*%c.val
}else{
	LLHc<-0
}
# starting values
if(start.on.true && use.generated){
	T<-T.r
	f<-f.r
	Delta.r=0;
	u<-u.r;
	v<-v.r;
	y<-y.r;
	yprime <- (y>0);
	w <- w.r;
	x <- c(u, w);
	beta<-beta.r
}

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

######
if(use.intercept){
	R <- makeRuint(dim=dimension, hyper=hyper, Q=Q, ku=Ku);
	# cholR <- chol.spam(R, memory=list(nnzcolindices=4e6),);
	cholR <<- get.cholMat(R,cholR);
}else{
	cholR<-NULL;
}

write.table(t(u), "usamples.txt", sep="\t",col.names=FALSE,row.names=FALSE)
write.table(t(y), "ypsamples.txt", sep="\t",col.names=FALSE,row.names=FALSE)

## sampler
ItTestNum<-gibbsit(NULL,NminOnly=TRUE);
cat("init ItTestNum:",ItTestNum,"\n",file="convergence_tests.txt")
if(final.run){
	# beginEstimate comes after adjustment of the sampling variances
	nbsimul<-ItTestNum;
	beginEstimate<-nbsimul
	use.autostop<-TRUE
}else{
	use.autostop<-FALSE;
}
nbtraced=17;
LLHu<-llh.ugivQ(dimension,u,Q,K[1])
LLH<-llh.zgivw(w,zpos,zneg,bivect)
if(use.insp){
	LLHb<-sum(dbeta(beta,abeta,bbeta,log=TRUE));
}else{
	LLHb<-0;
}
if(use.cofactors){
	LLHc<-sum(dnorm(c.val,0,sqrt(1/Kc),log=TRUE));
	c.valsamp<-as.matrix(mat.or.vec(nbsimul+1,nbfact.gen));
}else{
	LLHc<-sum(dnorm(rep(0,nbfact.gen),mean=0,sd=sqrt(1/Kc),log=TRUE));
}
LLHv<-sum(dnorm(v,0,Kv,log=TRUE));
if(fit.spatstruct){
	sampled<-as.matrix(mat.or.vec(nbsimul+1,nbtraced));
	sampled[1,1]<-T;
	sampled[1,2]<-LLHu
	sampled[1,3]<-f;
	sampled[1,4]<-LLHu
	sampled[1,5]<-Ku;
	sampled[1,6]<-llh.zgivy(y,zpos,zneg,bivect);
	sampled[1,7]<-LLH;
	sampled[1,8]<-llh.ygivw(y,w);
	sampled[1,9]<-0
	sampled[1,10]<-Kv;
	sampled[1,11]<-mean(u);
	sampled[1,12]<-Kc;
	sampled[1,13]<-LLHv;
	sampled[1,14]<-LLHc;
	sampled[1,15]<-LLHb;
	sampled[1,16]<-0;
	sampled[1,17]<-intercept;
}else{
	grid.stab<-seq(1,length(w),ceiling(length(w)/5))# values of the field tested for stability, keep 5 values
	nbtraced<-2*(2+length(grid.stab))+4
	spacer<-(2+length(grid.stab))


	sampled<-as.matrix(mat.or.vec(nbsimul+1,nbtraced));
	sampled[1,1]<-mean(u)
	sampled[1,2]<-sd(u)
	sampled[1,3:spacer]<-u[grid.stab]
	sampled[1,spacer+1]<-mean(w)
	sampled[1,spacer+2]<-sd(u)
	sampled[1,(spacer+3):(2*spacer)]<-w[grid.stab]
	LLHu<-llh.ugivQ(dimension,u,Q,K[1])
	sampled[1,(2*spacer)+1]<-llh.ugivQ(dimension,u,Q,K[1])
	sampled[1,(2*spacer)+2]<-llh.ygivw(y,w);
	sampled[1,(2*spacer)+3]<-llh.zgivy(y,zpos,zneg,bivect);
	sampled[1,(2*spacer)+4]<-llh.zgivw(w,zpos,zneg,bivect);
}



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

adaptOK<-FALSE
if(use.cosamplingfT){
	firstAdaptOK<-FALSE
}else{
	firstAdaptOK<-TRUE
}



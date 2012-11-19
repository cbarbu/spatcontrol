# get a sample
source("pseudo_data_generation.r")

Q.r<-QfromDelta(Delta.r,dist_mat,AS,f.r,addeps=0);

# precision matrix of x, x~N(0,S)
makeS <- function(dim,Q,K) {
 Ku <- K[1];
 Kv <- K[2];
 S <- Ku*Q;
 diag.spam(S) <- diag.spam(S) + Kv;
 S <- cbind(S, diag.spam(-1*Kv,dim,dim));
 S <- rbind(S, cbind(diag.spam(-1*Kv,dim,dim), diag.spam(Kv, dim, dim)))
 return(S);
}

# likelihood of original Delta
S<-makeS(dimension,Q.r,K.r)
exp_part<- x.r%*%S%*%x.r
epsilon<- 0.00000000001
cholS<-chol(S+diag.spam(epsilon,2*dimension))
detS<-2*(determinant(cholS)$modulus)
cat("exp_part",exp_part,"detS",detS,"LLH",-(exp_part+detS)/2);

# likelihood of x.r function of parameters "up" in the model
Deltas<-seq(0,threshold,10)
exp_part<-rep(0,length(Deltas))
detS<-rep(0,length(Deltas))
for(i in 1:length(Deltas)){
	Q<-QfromDelta(Deltas[i],dist_mat,AS,f.r,addeps=0);
	S<-makeS(dimension,Q,K.r)
	cholS<-chol(S+diag.spam(epsilon,2*dimension))
	exp_part[i]<- x.r%*%S%*%x.r
	detS[i]<-2*(determinant(cholS)$modulus)
	cat("exp_part",exp_part[i],"detS",detS[i],"LLH",-(exp_part[i]+detS[i])/2,"\n");
	cat("det: mod",determinant(cholS)$modulus,"sign",determinant(cholS)$sign);
}

LLH<- -(detS + exp_part)/2;

plot(Deltas,detS)
plot(Deltas,exp_part)
plot(c(min(Deltas),max(Deltas)),c(max(LLH),max(LLH-100)),type="n")
lines(Deltas,LLH,type="b")
# plot(Deltas,LLH,type="b")

# the problem is that without an additional term on the diagonal for Q at least, no way to calculate the LLH as Q is not inversible and its determinant is then null and the log of this is infinitly negative
# so we would need a way to get a proportional estimate of the likelihood of u~(0,Q) for det(Q) null and changing Q and the question is if adding a term to the diagonal is legitimate


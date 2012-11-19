# get a sample
source("pseudo_data_generation.r")

samplex <- function(dim,Q,K,y,cholR=NULL){
  center <- c(rep(0,dim),y);
  R <- makeR(dim,Q,K);
  # cholR<-chol(R)
  x <- rmvnorm.canonical(n=1, b=center, Q=R, Rstruct=cholR); # exploit the stability of the structure of R, at least 10 times faster.
  # cholR has to be computed at the initiation with the same structure
  return(drop(x));
}
# test llh of samples with same parameters
dimension<-nrow(Q.r)
nbsamp=100
LLHsame<-rep(0,nbsamp);
R.r <- makeR(dimension,Q.r,K.r);
cholR<-chol(R.r)
center<-c(rep(0,dimension),y.r);
for(it in 1:nbsamp){
	x <- rmvnorm.canonical(n=1, b=center, Q=R.r,Rstruct=cholR); # exploit the stability of the structure of R, at least 10 times faster.
	mu_r <- drop(solve.spam(cholR, center)) 
	X<-drop(x-mu_r);
	LLHsame[it]<-drop(t(X)%*%R.r%*%X)*(-1/2);
}
hist(LLHsame)
medSame<-median(LLHsame)
abline(v=medSame)
cat("medSame:",medSame,"\n");

LLHdiff<-rep(0,nbsamp);
Delta<-0
f<- 0.5
Q<-QfromDelta(Delta,dist_mat,AS,f);
scan()
R <- makeR(dimension,Q,K.r);
cholR<-chol(R.r)
for(it in 1:nbsamp){
	# x <- samplex(dimension,Q,K.r,y.r,cholR);
	x <- rmvnorm.canonical(n=1, b=center, Q=R,Rstruct=cholR); 
	mu_r <- drop(solve.spam(cholR, center)) 
	X<-drop(x-mu_r);
	LLHdiff[it]<-drop(t(X)%*%R%*%X)*(-1/2);
}
hist(LLHdiff)
medDiff<-median(LLHdiff)
abline(v=medDiff)
cat("medDiff:",medDiff,"\n");

##  results
# for Ku.r=Kv.r=1
# threshold = 100
# Delta.r=10 Delta = 50
# in all case the difference on the simplified LLH is minor but:
# f = 0.5, the diff is even inverse of the expected
# f = 1  the diff is null
# f = 5 the diff is weak but as expected
# Delta.r=50 Delta = 10
# f = 0.5, the diff is null
# f = 1  the diff is null
# f = 5  the diff is null

# threshold = 50
# Delta.r=10 Delta = 50
# f = 0.5  the diff is null
# f = 5  the diff is null
# Delta.r=50 Delta = 10
# f = 0.5  the diff is null
# f = 5  the diff is null

# for Ku.r=Kv.r=100
# threshold = 100
# Delta.r=10 Delta = 50
# f = 0.5, null or weak
# f = 1  the diff is null
# f = 5 the diff weak or null

# for Ku.r=Kv.r=0.1
# threshold = 100
# Delta.r=10 Delta = 50
# f = 1 null
# Delta.r=10 Delta = 50
# f = 1 null
# Delta.r=100 Delta = 10
# f = 1 null
# f = 5 weak or null

# study of influence of identif of f
# threshold = 100
# for Ku.r=Kv.r=1
# Delta = 0
# f.r = 1 f = 5 -> null
# f.r = 5 f = 1 -> null
# f.r = 0.5 f = 5 -> null
# for Ku.r=Kv.r=100
# f.r = 0.5 f = 5 -> null

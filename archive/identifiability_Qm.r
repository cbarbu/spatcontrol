# get a sample
source("pseudo_data_generation.r")

# let us use Q.r as the mean matrix
# the following are the residual of the fitting of y.r by the spatial mean
resV<-(Q.r-diag(1,dimension))%*%y.r
# and this the sum of the square of the residuals we gone use as likelihood
resS<- t(resV)%*%resV

# the other possibility is to use Qm with a normalisation by row
Qw<- -Q.r
diag(Qw)<-0;
sbyrow<-rsum(Qw)
matNorm<-0*Q.r
diag(matNorm)<-1/sbyrow
Qm<-matNorm%*%Qw
# plot(rsum(Qm))
# obviously the symmetry is largely lost but for our purpose we don't care

# an other possibility is to ponderate the residual count by the isolation of the house, then
# resV<-MatSum%*%(I-MatNorm%*%Qm) # where MatSum is the matricial inverse of MatNorm
matSum<-0*Q.r
diag(matSum)<-sbyrow
Id<-matSum%*%matNorm
resVmat<-matSum%*%Id-Qw 
# and this is Q.r !!!! at least I understood the hypothesis behind Q.r !!!!!
# Qm is the matrix to get the residuals of the spatial fit, ponderated by 
# the isolation of the house
# - extrapolation by simple normalized mean ponderated by weights
# inv(matSum)%*%Qw
# - we take the residuals of the fit of y by this mean
# y-inv(matSum)%*%Qw%*%y = (Id-inv(matSum)%*%Qw)%*%y
# - we ponderate the weight of the residuals in the sum of residuals by
# the isolation of the house:
# matSum%*%(Id-inv(matSum)%*%Qw)%*%y
# = (matSum%*%Id-Qw)%*%y
# = Q.r%*%y
# we then are completely legitimate in using Q.r%*%y as the likelihood 
# of u, it is its construction! and the sum of the square of the residuals is:
res<- Q.r%*%y.r
sumSqRes <- t(res)%*%res
# that has to be minimized (NB: -1/2*sumSqRes probably give a loglikelihood according to a normal law
cat("sum of squares Delta.r",sumSqRes,"LLH?",-sumSqRes/2,"\n");

# LLH? fn of Delta
Deltas<-seq(0,100,10)

detQs<-rep(0,length(Deltas))
sumSqRes<-rep(0,length(Deltas))
for(i in 1:length(Deltas)){
	Q<-QfromDelta(Deltas[i],dist_mat,AS,f.r);
	res<- Q%*%y.r
	sumSqRes[i] <- t(res)%*%res
	detQs[i]<-abs(determinant(Q)$modulus)
}
plot(Deltas,sumSqRes)
plot(Deltas,detQs)
plot(Deltas,detQs+sumSqRes)

# this is not working as is, increasing Delta make directly a smaller sum of squares
# should be corrected by the determinant of Q



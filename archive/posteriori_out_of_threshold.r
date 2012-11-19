# assess the importance what is cut using the threshold at 100m 
# a posteriori using meanf and meanT
# for each house maximize the share of autocorrelation from above the threshold
# by counting how many houses are not included, affecting them the weight at 100m
# and see how the sum of these weight compare with the sum of the weight within the threshold
load("EndVisu.img")
source("spam_complement.r")
source("parameters_sampler.r")
source("functions_intercept.r")
source("DeltaSampling.r")
T<-0.02
f<-1
Kernel<-cauchyKernel

#### sum of weigths under threshold
ones<-rep(1,dim(Q)[1])
Q<- SpecificMultiply.spam(1/T,Kernel(T,Dmat,f),SB);#the precision matrix of u 
Q<- SpecificMultiply.spam(1/T,Kernel(T,dist_mat_large,f),SBL);#the precision matrix of u 
diag(Q)<-0
weightIn<-Q%*%ones

#### maximization of weights above threshold
# weights for a larger neighborhood
bigThreshold<-3*threshold
dist_mat_large <-nearest.dist(x=data[,c("easting","northing")], y=NULL, method="euclidian", delta=bigThreshold, upper=NULL);          
diag(dist_mat_large)<-0
dist_mat_large<-as.spam(dist_mat_large)

dist_mat_large_only<-dist_mat_large
dist_mat_large_only@entries[dist_mat_large@entries<threshold]<-rep(0,length(which(dist_mat_large@entries<threshold)))
dist_mat_large_only<-as.spam(dist_mat_large_only)

SBL <- nearest.dist(x=cbind(data$block_num,rep(0,length(data$block_num))), method="euclidian", upper=NULL,delta=0.1)
SBL@entries<-rep(1,length(SBL@entries))
diag(SBL)<-0
dmtL<-dist_mat_large_only
dmtL@entries<-rep(1,length(dmtL@entries))
SBL<-SBL*dmtL
SBL<-as.spam(SBL)

QL<- SpecificMultiply.spam(1/T,Kernel(T,dist_mat_large_only,f),SBL);#the precision matrix of u 
diag(QL)<-0
weightInL<-QL%*%ones

# nb of neighbors out large
Q1<-dist_mat_large
Q1@entries<-rep(1,length(Q1@entries))
nbInL<-Q1%*%ones
nbOutL<-dim(Q)[1]-nbIn

# max influence of these neighbors
MaxWeightOut<-nbOutL*Kernel(T,bigThreshold,f)+weightInL

#### max of neglected influence
neglectedInf<-MaxWeightOut/(MaxWeightOut+weightIn)
cat("mean of neglected influence was a posteriori:",mean(neglectedInf),"\n")
hist(neglectedInf)


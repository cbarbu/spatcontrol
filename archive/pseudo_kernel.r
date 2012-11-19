u<-rep(0,dimension)
zknown<-which(block_data$status!=9)
u[zposb]<-valPos
u[znegb]<-valNeg
nbzknown<-length(zknown)
nbzNAb<-length(zNAb)
### calculate a u only for unknown given 0 or 1
## subset Q
QgetNA<-as.matrix(Qdiff[zNAb,zknown])
## normalize the lignes to 1 the weights
sumLigne<-QgetNA%*%rep(1,nbzknown)
existNeigh<-which(sumLigne!=0) # to avoid dividing by zero	
# NormalizingMat<-mat.or.vec(nbzNAb,nbzNAb)
# diag(NormalizingMat)[existNeigh]<-1/sumLigne[existNeigh]
# 
# QgetNA<-NormalizingMat%*%QgetNA
# MeanKnown<-mean(u[zknown]) # rate of positives in 0/1

u.zNAb<-QgetNA%*%u[zknown]
if(length((u[zNAb])[-existNeigh])>0){
	u.zNAb[-existNeigh]<-uOfNoNeigh
}
u[zNAb]<-u.zNAb
est.u.b<-u


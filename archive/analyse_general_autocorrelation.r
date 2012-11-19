load("EndVisu_good.img")
source("spam_complement.r")
source("functions_intercept.r")
source("DeltaSampling.r")

distances = seq(5,125,20)
nbRep<-10
# change parameters estimates
meanf<-1.05*9.34
meanT<-0.7*0.27
est.mu<- -1.12
est.Ku<- 0.38
use.NormQ <- FALSE
epsilon <- 0.003
# est.Kv <- 212.97; # rgamma(n=1, shape=K.hyper[3], scale=K.hyper[4]);

est.yp<-(est.y>0)
# distances = seq(15,140,15)

# data<-data.base
## autocorrelation in true status for observed locations
sel<-(data$status!=9)
mats_neigh<-gen.mats.neigh(distances,data$easting[sel],data$northing[sel],data$block_num[sel])
mIrefest<-structured.moransI(distances,mats_neigh,est.yp[sel],nb_rep_sign=0);

plot.structured.moransI(distances,mIrefest);

## autocorrelation in observed status for observed locations

dev.new()
par(mfcol=c(2,5))

est.Q<-QfromfT(Dmat,AS,SB,f=meanf,T=meanT)
est.detection<-inspector%*%est.beta
noNA<-which(data$status!=9)

MIs<-generated.morans.struct(distances,mats_neigh,nbRep,est.detection=est.detection,est.Q=est.Q,est.Ku=est.Ku,est.mu=est.mu,est.c.comp=est.c.comp,est.Kv=est.Kv,true.val=data$status)
stop()

est.Q<-QfromfT(Dmat,AS,SB,f=meanf,T=1)
MIsGap<-generated.morans.struct(distances,mats_neigh,nbRep,est.detection=est.detection,est.Q=est.Q,est.Ku=est.Ku,est.mu=est.mu,est.c.comp=est.c.comp,est.Kv=est.Kv,est.v=est.v,true.val=data$status,trueStatus=FALSE)

# plot share gap effect
medGapEffect<-apply(MIsGap$MI2-MIsGap$MI3,2,median)
medBase<-apply(MIs$MI2-MIs$MI3,2,median)
# medNoCof<-apply(MIsNoCof$MI2-MIsNoCof$MI3,2,median)
# medNoInsp<-apply(MIsNoInsp$MI2-MIsNoInsp$MI3,2,median)
# medNoCofInsp<-apply(MIsNoCofNoInsp$MI2-MIsNoCofNoInsp$MI3,2,median)
dev.new()
plot(medBase,ylim=c(0,max(medBase)),type="l",xaxt="n",ylab="IS-ID",xlab="Distance (m)")
get.med.position.axis(distances)
lines(medGapEffect,lty=2)
# lines(medNoCof,lty=3)
# lines(medNoInsp,lty=4)
# lines(medNoCofInsp,lty=2)

# plot correlation structure in est.Q
# careful, this cannot apply as is to the whole map, too big
dev.new()
par(mfcol=c(2,2))
est.Q<-QfromfT(Dmat,AS,SB,f=meanf,T=meanT)
mats_neigh<-gen.mats.neigh(distances,data$easting,data$northing,data$block_num)
MIsGap<-generated.morans.struct(distances,mats_neigh,nbRep,est.detection=est.detection,est.Q=est.Q,est.Ku=est.Ku,est.mu=est.mu,est.c.comp=est.c.comp,est.Kv=est.Kv,est.v=est.v,true.val=data$status)
StructCorrel(distances,Q=est.Ku*est.Q,mats_neigh=mats_neigh)
CovMat<-solve(est.Ku*est.Q)
SC<-StructCorrel(distances,CovMat=CovMat,mats_neigh=mats_neigh)
plot(med_position,sqrt(SC[-1,1]),ylim=c(0,max(sqrt(SC[-1,]))))
lines(med_position,sqrt(SC[-1,2]))
lines(med_position,sqrt(SC[-1,3]))
plot(med_position,SC[-1,1],ylim=c(0,max(SC[-1,])))
lines(med_position,SC[-1,2])
lines(med_position,SC[-1,3])
NormVar<-1/diag(CovMat)
NormMat<-as.spam(0*CovMat)
diag(NormMat)<-NormVar
CorMat<-NormMat%*%CovMat
SCor<-StructCorrel(distances,CovMat=CorMat,mats_neigh=mats_neigh)
plot(med_position,sqrt(SCor[-1,1]),ylim=c(0,1))
lines(med_position,sqrt(SCor[-1,2]))
lines(med_position,sqrt(SCor[-1,3]))
plot(med_position,SCor[-1,1],ylim=c(0,1))
lines(med_position,SCor[-1,2])
lines(med_position,SCor[-1,3])

# get autocorrelation

# autocorrelation on whole map
dev.new()
par(mfcol=c(2,4))
period<- "" # keep the all map, for the test
source("import_data.r")
est.detectionFull<-rep(mean(est.detection),dim(data)[1])

dist_matFull <-nearest.dist(x=data[,c("easting","northing")], y=NULL, method="euclidian", delta=threshold, upper=NULL);          
diag(dist_matFull)<- rep(-1,dim(dist_matFull)[1])
diagDmat<-which(dist_matFull@entries==-1) 
dist_matFull@entries[diagDmat]<-rep(0,length(diagDmat))

mats_neighQ<-gen.mats.neigh(c(0,threshold),data$easting,data$northing,data$block_num)
DmatFull<-dist_matFull
SBFull<-mats_neighQ[[2]][["SBr"]]
ASFull<-mats_neighQ[[2]][["ASr"]]


# data.base$oanimal<-as.integer((data.base$CO==1 | data.base$AV==1 | data.base$GA==1 | data.base$OV==1 | data.base$otros.animales!="-1"))
c.mapFull<-as.matrix(data[,cofs])
est.c.compFull<-c.mapFull%*%est.c.val

sel<-which(data$status!=9)
mats_neighM<-gen.mats.neigh(distances,data$easting[sel],data$northing[sel],data$block_num[sel])

est.QFull<-QfromfT(DmatFull,ASFull,SBFull,f=meanf,T=meanT)
MIs<-generated.morans.struct(distances,mats_neighM,nbRep,est.detection=est.detectionFull,est.Q=est.QFull,est.Ku=est.Ku,est.mu=est.mu,est.c.comp=est.c.compFull,est.Kv=est.Kv,true.val=data$status)

# MIsNoCof<-generated.morans.struct(distances,mats_neighM,nbRep,est.QFull,est.Ku,est.mu,0*est.c.compFull,est.Kv,0*est.c.compFull,est.detectionFull,true.val=data$status)

est.QFull<-QfromfT(DmatFull,ASFull,SBFull,f=meanf,T=1)
MIsGap<-generated.morans.struct(distances,mats_neighM,nbRep,est.detection=est.detectionFull,est.Q=est.QFull,est.Ku=est.Ku,est.mu=est.mu,est.c.comp=est.c.compFull,est.Kv=est.Kv,true.val=data$status)

# plot share gap effect
medGapEffect<-apply(MIsGap$MI2-MIsGap$MI3,2,median)
medBase<-apply(MIs$MI2-MIs$MI3,2,median)
# medNoCof<-apply(MIsNoCof$MI2-MIsNoCof$MI3,2,median)
dev.new()
plot(medBase,ylim=c(0,max(medBase)),type="l",xaxt="n",ylab="IS-ID",xlab="Distance (m)")
get.med.position.axis(distances)
lines(medGapEffect,lty=2)
# lines(medNoCof,lty=2)




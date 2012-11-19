# to dump, use in the R session that create the results something like
# ExpSSampled<-smallsampled
# GroupWeightExp<-GroupWeight
# CRExp<-CR
# dump(c("ExpSSampled","GroupWeightExp","CRExp"),file="ExpDump.R")

nameSimul<-"mergegraphs_Pau2007_cau_epsKc0.01StreetsTr100_Full_Kpseudoflat_flatInsp_noNorm_broadfT_on_Pau2007_rep300_0-600-15"
lisFit<-1
distances = seq(0,600,15)
city<-"PAUCARPATA"
subsetcity<-0
period<- "fall.2007" 
nbRep<-300

threshold<-100
use.NormQ<-FALSE
epsilon <- 0.01

use.map.gen<-FALSE
use.NA<-TRUE
value.NA<-9
Models<-c("Exp","Gau","Cau","Geo")

library(locfit)
plotFits<-function(x,listFits,topLabel=TRUE,...){
	pred<-data.frame(x)
	med<-list()
	maxDensity<-0
	for(model in names(listFits)){
		pred[[model]]<-predict(listFits[[model]],x)
		maxDensity<-max(maxDensity,max(pred[[model]]))
		med[[model]]<-x[which(pred[[model]]==max(pred[[model]]))]
	}

	plot(1,1,type="n",xlim=c(min(x),max(x)),ylim=c(0,maxDensity),...)
	num<-1
	legendNames<-c()
	legendLty<-c()
	for(model in names(listFits)){
		lines(x,pred[[model]],lty=num)
		if(topLabel){
			mtext(model,side=3,at=med[[model]])
		}
		legendNames<-c(legendNames,model)
		legendLty<-c(legendLty,num)
		num<-num+1
	}
	legend("topleft",legendNames,lty=legendLty)
	return(invisible(pred))
}
source("ExpDump.R")
source("CauDump.R")
source("GauDump.R")
source("GeoDump.R")

par(mfrow=c(2,2))
# basic fit
Tfit<-list()
Tfit$Exp<-locfit(~ExpSSampled[,1],xlim=c(0,15),alpha=lisFit)
Tfit$Cau<-locfit(~CauSSampled[,1],xlim=c(0,15),alpha=lisFit)
Tfit$Gau<-locfit(~GauSSampled[,1],xlim=c(0,15),alpha=lisFit)
Tfit$Geo<-locfit(~GeoSSampled[,1],xlim=c(0,15),alpha=lisFit)
plotFits(seq(0,1.2,0.01),Tfit,ylab="density",xlab="Direct effect of streets (sigma)")

# fit on logs
Tfitlog<-list()
Tfitlog$Exp<-locfit(~log(ExpSSampled[,1]),alpha=lisFit)
Tfitlog$Cau<-locfit(~log(CauSSampled[,1]),alpha=lisFit)
Tfitlog$Gau<-locfit(~log(GauSSampled[,1]),alpha=lisFit)
Tfitlog$Geo<-locfit(~log(GeoSSampled[,1]),alpha=lisFit)

plotFits(seq(-5,1,0.01),Tfitlog,ylab="density",xlab="Direct effect of streets log(sigma)")

ffit<-list()
ffit$Exp<-locfit(~ExpSSampled[,3],alpha=lisFit)
ffit$Cau<-locfit(~CauSSampled[,3],alpha=lisFit)
ffit$Gau<-locfit(~GauSSampled[,3],alpha=lisFit)
ffit$Geo<-locfit(~GeoSSampled[,3],alpha=lisFit)
plotFits(seq(0,50,0.01),ffit,ylab="density",xlab="Characteristic distance (d)")

ffitLog<-list()
ffitLog$Exp<-locfit(~log(ExpSSampled[,3]),alpha=lisFit)
ffitLog$Cau<-locfit(~log(CauSSampled[,3]),alpha=lisFit)
ffitLog$Gau<-locfit(~log(GauSSampled[,3]),alpha=lisFit)
ffitLog$Geo<-locfit(~log(GeoSSampled[,3]),alpha=lisFit)
plotFits(seq(0,4,0.01),ffitLog,ylab="density",xlab="Characteristic distance (d)")

# then should show the best fit of the autocorrelation for the three kernels
# should use Ku*exp(f) of even something integrative, typically the shape of 
# of the auto-correlation of the true infestation, given est.u and others

est.Ku<-list()
est.Ku$Exp<-mean(ExpSSampled$Ku)
est.Ku$Cau<-mean(CauSSampled$Ku)
est.Ku$Gau<-mean(GauSSampled$Ku)
est.Ku$Geo<-mean(GeoSSampled$Ku)

est.f<-list()
est.f$Exp<-mean(ExpSSampled$f)
est.f$Cau<-mean(CauSSampled$f)
est.f$Gau<-mean(GauSSampled$f)
est.f$Geo<-mean(GeoSSampled$f)

est.T<-list()
est.T$Exp<-mean(ExpSSampled$T)
est.T$Cau<-mean(CauSSampled$T)
est.T$Gau<-mean(GauSSampled$T)
est.T$Geo<-mean(GeoSSampled$T)

est.mu<-list()
est.mu$Exp<-mean(ExpSSampled$mu)
est.mu$Cau<-mean(CauSSampled$mu)
est.mu$Gau<-mean(GauSSampled$mu)
est.mu$Geo<-mean(GeoSSampled$mu)

est.Kv<-list()
est.Kv$Exp<-mean(ExpSSampled$Kv)
est.Kv$Cau<-mean(CauSSampled$Kv)
est.Kv$Gau<-mean(GauSSampled$Kv)
est.Kv$Geo<-mean(GeoSSampled$Kv)

source("functions_intercept.r")
# look at the Moran's I for "TrueStatus"
source("spam_complement.r")
source("DeltaSampling.r")

source("import_data.r")
cat("import ok\n")

dist_mat <-nearest.dist(x=data[,c("easting","northing")], y=NULL, method="euclidian", delta=threshold, upper=NULL);          
diag(dist_mat)<-0
Dmat<-dist_mat
dmtr<-dist_mat
dmtr@entries<-rep(1,length(dmtr@entries))
mats_neighQ<-gen.mats.neigh(c(0,threshold),data$easting,data$northing,data$block_num)
SB<-mats_neighQ[[2]][["SBr"]]
AS<-mats_neighQ[[2]][["ASr"]]

### look at the structure of the weights with distance
# mats_neigh<-gen.mats.neigh(c(0,distances),data$easting,data$northing,data$block_num)
# stNeigh<-StructNeigh(c(0,distances),mats_neigh)
# stNeigh$SBratio<-stNeigh$nbNeighSB/stNeigh$nbNeighTotal
# plot(stNeigh$med_position,stNeigh$SBratio,type="l")
# dev.new()
# par(mfrow=c(1,3))
# x<-stNeigh$med_position
# val<-data.frame(x)
# val$Exp<-expKernel(1,x,est.f$Exp)
# val$Cau<-cauchyKernel(1,x,est.f$Cau)
# val$Gau<-gaussianKernel(1,x,est.f$Gau)
# val$Geo<-geometricKernel(1,x,est.f$Geo)
# plot(x,x,type="n",ylim=c(0,1),ylab="Brut weigth within same block")
# lines(x,val$Exp,lty=1)
# lines(x,val$Gau,lty=2)
# lines(x,val$Cau,lty=3)
# lines(x,val$Geo,lty=4)
# plot(x,x,type="n",ylim=c(0,1),ylab="Brut weight across street(s)")
# lines(x,val$Exp*est.T$Exp,lty=1)
# lines(x,val$Gau*est.T$Gau,lty=2)
# lines(x,val$Cau*est.T$Cau,lty=3)
# lines(x,val$Geo*est.T$Geo,lty=4)
# plot(x,x,type="n",ylim=c(0,1),ylab="Mean brut weight according to SB/AS distribution")
# lines(x,val$Exp*(est.T$Exp*(1-stNeigh$SBratio)+stNeigh$SBratio),lty=1)
# lines(x,val$Gau*(est.T$Gau*(1-stNeigh$SBratio)+stNeigh$SBratio),lty=2)
# lines(x,val$Cau*(est.T$Cau*(1-stNeigh$SBratio)+stNeigh$SBratio),lty=3)
# lines(x,val$Geo*(est.T$Geo*(1-stNeigh$SBratio)+stNeigh$SBratio),lty=4)
# ###

est.Q<-list()
est.Q$Exp<-QfromfT(Dmat,AS,SB,f=est.f$Exp,T=est.T$Exp,kern=expKernel)
est.Q$Cau<-QfromfT(Dmat,AS,SB,f=est.f$Cau,T=est.T$Cau,kern=cauchyKernel)
est.Q$Gau<-QfromfT(Dmat,AS,SB,f=est.f$Gau,T=est.T$Gau,kern=gaussianKernel)
est.Q$Geo<-QfromfT(Dmat,AS,SB,f=est.f$Geo,T=est.T$Geo,kern=geometricKernel)
cat("est Q ok\n")

# ### evolution of the sum of the cumulative sum of the spatial weights with distance
# dev.new()
# par(mfrow=c(2,2))
# PartVectRat<-list()
# for(model in Models){
# 	dimension<-dim(est.Q[[model]])[1]
# 	PartVect<-mat.or.vec(length(distances),dimension)
# 	tot<-as.spam(mat.or.vec(dimension,dimension))
# 	for (i in 1:(length(distances))){
# 		dmtr<-mats_neigh[[i+1]][[1]]
# 		diag(dmtr)<-0
# 		tot<-tot+dmtr
# 		Qdist<-as.matrix(dmtr)*as.matrix(est.Q[[model]])
# 		PartVect[i,]<-apply(Qdist,1,sum)
# 	}
# 	PartVectRat[[model]]<- -PartVect %*% (1/(diag(est.Q[[model]])-epsilon))/dimension
# 	plot(stNeigh$med_position,PartVectRat[[model]],type="l")
# 	mtext(model,side=3)
# }
# dev.print(device=pdf,"PartVect.pdf")

# par(mfrow=c(2,2))
# for(model in Models){
# 	plot(stNeigh$med_position,cumsum(PartVectRat[[model]]),type="l",ylim=c(0,1),xlab="Distance (m)",ylab="Mean cumulative participation")
# 	abline(h=0.95,lty=2)
# 	mtext(model,side=3)
# }
# dev.print(device=pdf,"PartVectCumul.pdf")

#

# est.WSB<-list()
# for(model in Models){
# finalWeight<- -est.Q[[model]]
# diag(finalWeight)<-0
# sumW<-apply_by_row_not_null.spam(finalWeight,sum)
# sumWSB<-apply_by_row_not_null.spam(SB*finalWeight,sum)
# est.WSB[[model]]<-mean(sumWSB/sumW)
# }

### look at the structure of the partial correlation with distance
### can only be done for relatively small subset (fall2007)
# according to Schaefer2005 the partial correlation is given by
library("corpcor")
source("bugfixMultipleMatSpam.R")
# mSB<-as.matrix(SB)
est.prSB<-list()
dev.new()
par(mfcol=c(3,4))
for(model in Models){
	b<-decompose.invcov(as.matrix(est.Ku[[model]]*est.Q[[model]])) # list(matrix of partial correlations,vector variances)
	pr.spam<-as.spam(b$pr*dmtr) # keep only entries matching Dmat
	plot(Dmat@entries,pr.spam@entries,pch=".",xlab="Distance (m)",ylab="Estimated partial correlation",col=(dmtr+2*SB)@entries)

	# check the relationship first neighbor across street dist
	# and partial correlation
	SBnoself<-SB
	diag(SBnoself)<-0
	MinSBdist<-apply_by_row_not_null.spam(as.spam(Dmat*SBnoself),min)
	MaxSBpr<-apply_by_row_not_null.spam(b$pr*SBnoself,max)
	MinASdist<-apply_by_row_not_null.spam(Dmat*AS,min)
	MaxASpr<-apply_by_row_not_null.spam(pr.spam*AS,max)
	plot(MinASdist,MaxASpr,ylim=c(0,1),xlim=c(0,100),pch=".")
	lines(MinSBdist,MaxSBpr,type="p",col=4,pch=".")

	# get in general link distance to partial correlation
	DiagQ<-diag(est.Q[[model]])[Dmat@colindices]
	stDiagQ<-round(10*(DiagQ-min(DiagQ))/(max(DiagQ)-min(DiagQ)))
	colDiag<-heat.colors(11)[stDiagQ+1]
	plot(Dmat@entries,pr.spam@entries,pch=".",xlab="Distance (m)",ylab="Estimated partial correlation",col=colDiag)

	diag(b$pr)<-0
	# plot(apply(b$pr,2,sum),pch=3)

	sumprSB<- (b$pr*SB) %*% rep(1,dim(SB)[1])
	sumpr<-b$pr %*% rep(1,dim(b$pr)[1])
	est.prSB[[model]]<-mean(sumprSB/sumpr)
}
stop()
# dev.print(device=png,"PartialCorrelationStructure.png",width=1600)
# ###

# cov<-solve(est.Ku$Exp*est.Q$Exp)
# a<-decompose.cov(cov)
# r2<-a$r*a$r
# r2d0<-r2
# diag(r2d0)<-0
# sumCor<-apply(r2d0,2,sum)
# sumCorSB<-apply(SB*r2d0,2,sum)
# mean(sumCorSB/sumCor)
# 
# small<-which(cov<0.01)
# diag.cov<-diag(cov)
# diag(cov)<-0

dev.new()
par(mfcol=c(2,4))
sel<-which(data$status!=9)
mats_neigh<-gen.mats.neigh(distances,data$easting[sel],data$northing[sel],data$block_num[sel])
MIs<-list()
MIs$Exp<-generated.morans.struct(distances,mats_neigh,nbRep,est.Q=est.Q$Exp,est.Ku=est.Ku$Exp,est.mu=est.mu$Exp,est.Kv=est.Kv$Exp,trueStatus=FALSE,true.val=data$status)
mtext("Exp",side=3)
MIs$Cau<-generated.morans.struct(distances,mats_neigh,nbRep,est.Q=est.Q$Cau,est.Ku=est.Ku$Cau,est.mu=est.mu$Cau,est.Kv=est.Kv$Cau,trueStatus=FALSE,true.val=data$status)
mtext("Cau",side=3)
MIs$Gau<-generated.morans.struct(distances,mats_neigh,nbRep,est.Q=est.Q$Gau,est.Ku=est.Ku$Gau,est.mu=est.mu$Gau,est.Kv=est.Kv$Gau,trueStatus=FALSE,true.val=data$status)
mtext("Gau",side=3)
MIs$Geo<-generated.morans.struct(distances,mats_neigh,nbRep,est.Q=est.Q$Geo,est.Ku=est.Ku$Geo,est.mu=est.mu$Geo,est.Kv=est.Kv$Geo,trueStatus=FALSE,true.val=data$status)
mtext("Geo",side=3)

save.image(file="MIs.img")

ymin1<- Inf
ymin2<- Inf
ymax1<- -Inf
ymax2<- -Inf
for(model in Models){
	for(column in c(1,length(distances)-1)){
		min1<-quantile(MIs[[model]][["MI1"]][,column],probs=c(0.025))
		ymin1<-min(ymin1,min1)
		min2<-quantile(MIs[[model]][["MI2"]][,column]-MIs[[model]][["MI3"]][,column],probs=c(0.025))
		ymin2<-min(ymin2,min2)
		max1<-quantile(MIs[[model]][["MI1"]][,column],probs=c(0.975))
		ymax1<-max(ymax1,max1)
		max2<-quantile(MIs[[model]][["MI2"]][,column]-MIs[[model]][["MI3"]][,column],probs=c(0.975))
		ymax2<-max(ymax2,max2)
	}
	cat(model,min1,max1,min2,max2,"\n")
}
par(mfcol=c(2,4))
for(model in Models){
replot.gen.struct(distances,MIs[[model]],ylim1=c(0,ymax1),ylim2=c(ymin2,ymax2))
mtext(model,side=3)
}
# dev.print(device=pdf,"FitSimulMoran.pdf")

# dev.new()
# par(mfcol=c(2,4))
# mats_neigh<-gen.mats.neigh(distances,data$easting,data$northing,data$block_num)
# MIs$Exp<-generated.morans.struct(distances,mats_neigh,nbRep,est.Q=est.Q$Exp,est.Ku=est.Ku$Exp,est.mu=est.mu$Exp,est.Kv=est.Kv$Exp,trueStatus=TRUE)
# mtext("Exp",side=3)
# MIs$Cau<-generated.morans.struct(distances,mats_neigh,nbRep,est.Q=est.Q$Cau,est.Ku=est.Ku$Cau,est.mu=est.mu$Cau,est.Kv=est.Kv$Cau,trueStatus=TRUE)
# mtext("Cau",side=3)
# MIs$Gau<-generated.morans.struct(distances,mats_neigh,nbRep,est.Q=est.Q$Gau,est.Ku=est.Ku$Gau,est.mu=est.mu$Gau,est.Kv=est.Kv$Gau,trueStatus=TRUE)
# mtext("Gau",side=3)
# MIs$Geo<-generated.morans.struct(distances,mats_neigh,nbRep,est.Q=est.Q$Geo,est.Ku=est.Ku$Geo,est.mu=est.mu$Geo,est.Kv=est.Kv$Geo,trueStatus=TRUE)
# mtext("Geo",side=3)

# # from same block
# dev.new()
# par(mfrow=c(1,2))
# GW<-list()
# GW$Exp<-GroupWeightExp
# GW$Cau<-GroupWeightCau
# GW$Gau<-GroupWeightGau
# GW$Geo<-GroupWeightGeo
# 
# rSBfitSpat<-list()
# rSBfitMean<-list()
# for(model in names(GW)){
# 	# spatial variety for estimated f/T
# 	normWeights<- -est.Q[[model]]
# 	diag(normWeights)<-0
# 	sumT<-normWeights%*%rep(1,length(data$status))
# 	sumSB<-(normWeights*SB)%*%rep(1,length(data$status))
# 	ratioSB<-sumSB/sumT
# 	rSBfitSpat[[model]]<-locfit(~ratioSB,xlim=c(0,1),alpha=lisFit)
# 
# 	# posterior of the mean
# 	rSBfitMean[[model]]<-locfit(~GW[[model]][,1],xlim=c(0,1),alpha=lisFit)
# }
# plotFits(seq(0,1,0.01),rSBfitMean)
# plotFits(seq(0,1,0.01),rSBfitSpat)



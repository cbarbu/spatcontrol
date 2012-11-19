# in case reloading
library("LaplacesDemon")
source("morans_functions.r")
source("spam_complement.r")
source("prep_sampling.r")
source("DeltaSampling.r")

#### parameters
PercentHousesRemoved<-5 # % of houses removed
nbRepeat<-5 # nb of test of PercentHousesRemoved houses
randomRemoved<-FALSE # resample PercentHousesRemoved or avoid previously removed (then limit nbRepeat to 100/PercentHousesRemoved)

#### init the kernel_fit_space.r
## general before any simulation
Ku<-mean(smallsampled$Ku)
Kv<-mean(smallsampled$Kv)
K<-c(Ku,Kv,Kc);
c.val<-est.c.val
f<-meanf
T<-meanT
Q.est<-QfromfT(dist_mat,AS,SB,meanf,meanT);

#### initialisation before each simulation
## choice of the data to be analyzed
# # choice of only one block
# block<-2008 
# sel<-which(data$block_num==block)
# choice of everything
sel<-(1:length(data$status))

# actual subsetting
block_data<-data[sel,]


## technical initialisation
QualityResult<-list()
known<-which(block_data$status<9)# choice only in known status
block_data$TrueStatus<-block_data$status # ref when choice n% of the houses to be set as unknown
if(!randomRemoved){
	nbHousesPerSlice<-round(PercentHousesRemoved*length(known)/100)
	nbRepeat<-min(ceiling(length(known)/nbHousesPerSlice),nbRepeat)
	TotalToGuess<-sample(known,min(round(nbRepeat*PercentHousesRemoved*length(known)/100),length(known)))
	cat("Will remove",nbRepeat,"times",nbHousesPerSlice,"for a total of",length(TotalToGuess))
}

## calculus
for(numRepeat in 1:nbRepeat){
	if(randomRemoved){
		ToGuess<-sample(known,round(PercentHousesRemoved*length(known)/100))
	}else{
		lowBound<-1+(numRepeat-1)*nbHousesPerSlice
		upBound<-min(lowBound+nbHousesPerSlice,length(TotalToGuess))
		ToGuess<-TotalToGuess[lowBound:upBound]
		cat("\ntesting",lowBound,"to",upBound,"of",length(TotalToGuess),"\n")
	}
	block_data$status[ToGuess]<-rep(9,length(ToGuess))

	## ploting the chosen data
	par(mfrow=c(2,2))
	block_visudata<-block_data$pos
	block_visudata[block_data$insp==0]<-0.5
	plot_reel(block_data$easting,block_data$northing,block_visudata,base=0,top=1)

	block_u_pred <- u_pred[sel]
	plot_reel(block_data$easting,block_data$northing,block_u_pred,base=0,top=max(block_u_pred))
	text(block_data$easting,block_data$northing,labels=round(block_u_pred,2),pos=4)

	# technical inititialisation according to the chosen dataset
	starter<-1
	dimension<-length(sel)
	w<-est.w[sel]
	u<-est.u[sel]
	bivect<-est.detection[sel]
	zposb<-which(block_data$status[sel]==1)
	znegb<-which(block_data$status[sel]==0)
	zNAb<-which( block_data$status[sel]==9)
	y<-as.integer(rnorm(length(w),mean=w,sd=1)>0)
	y[zposb]<-1
	Q<-Q.est[sel,sel]
	R <- makeR(length(u),Q,K);
	cholR <- chol.spam(R);
	cholQ<-chol.spam(Q)
	c.comp<-drop(c.map[sel,]%*%c.val)
	grid.stab<-seq(1,length(w),ceiling(length(w)/5))# values of the field tested for stability, keep 5 values

	ItTestNum<- gibbsit(NULL,NminOnly=TRUE);
	beginEstimate<-1
	AdaptOK<-TRUE

	nbsimul<-ItTestNum+beginEstimate
	nbtraced<-2*(2+length(grid.stab))+4
	spacer<-(2+length(grid.stab))

	sampledb<-as.matrix(mat.or.vec(nbsimul+1,nbtraced));
	sampledb[1,1]<-mean(u)
	sampledb[1,2]<-sd(u)
	sampledb[1,3:spacer]<-u[grid.stab]
	sampledb[1,spacer+1]<-mean(w)
	sampledb[1,spacer+2]<-sd(u)
	sampledb[1,(spacer+3):(2*spacer)]<-w[grid.stab]
	LLHu<-llh.ugivQ(dimension,u,Q,K[1])
	sampledb[1,(2*spacer)+1]<-llh.ugivQ(dimension,u,Q,K[1])
	sampledb[1,(2*spacer)+2]<-llh.ygivw(y,w);
	sampledb[1,(2*spacer)+3]<-llh.zgivy(y,zposb,znegb,bivect);
	sampledb[1,(2*spacer)+4]<-llh.zgivw(w,zposb,znegb,bivect);


	sum.u.b<-rep(0,length(u))
	sum.v.b<-rep(0,length(u))
	sum.w.b<-rep(0,length(u))
	sum.y.b<-rep(0,length(u))
	sum.yp.b<-rep(0,length(u))

	Rprof()
	source("kernel_fit_space.r")
	Rprof(NULL)
	## checking on one block
	# # par(mfrow=c(2,1))
	# # plot_reel(block_data$easting,block_data$northing,block_visudata,base=0,top=1)
	# u_pred.b<-pnorm(est.u.b,0,1)
	# plot_reel(data$easting[sel],data$northing[sel],pnorm(c.comp,0,1))
	# text(block_data$easting,block_data$northing,labels=round(pnorm(c.comp,0,1),2),pos=4)
	# plot_reel(data$easting[sel],data$northing[sel],u_pred.b)
	# text(block_data$easting,block_data$northing,labels=round(u_pred.b,2),pos=4)

	# ## general analysis
	# par(mfrow=c(2,2))
	# hist(est.w)
	# hist(est.w[data$status==9],add=T,col=4)
	# 
	# hist(pnorm(est.w))
	# hist(pnorm(est.w[data$status==9]),add=T,col=4)
	# 
	# hist(est.w.b)
	# hist(est.w.b[data$status==9],add=T,col=4)
	# 
	# hist(pnorm(est.w))
	# hist(pnorm(est.w[data$status==9]),add=T,col=4)

	# est.detection.b<-inspector[sel,]%*%est.beta
	# plot.prob.to.observed((pnorm(est.w.b,0,1)*est.detection.b)[ToGuess],block_data$TrueStatus[ToGuess])
	# plot(est.w.b[ToGuess],(est.w[sel])[ToGuess])
	# abline(a=0,b=1)
	# hist(est.yp.b[ToGuess])
	# plot.prob.to.observed(u_pred*est.detection,visudata)
	# plot.prob.to.observed(u_pred*est.detection,visudata)
	# ## quality of prediction for to be guessed
	QualityResult[[numRepeat]]<-plot.prob.to.observed(est.yp.b[ToGuess],block_data$TrueStatus[ToGuess],xlim=c(0,1),ylim=c(0,1))
}

#### analysis of the results
## transform the initial big list in digestible chunks
OverAllQuality<-QualityResult[[1]]
asMatQualityResult<-mat.or.vec(nbRepeat,length(QualityResult[[1]][[1]]))
asMatQualityResult[1,]<-QualityResult[[1]][[2]]
for(numRepeat in 2:nbRepeat){
	OverAllQuality[[2]]<-QualityResult[[numRepeat]][[2]]+OverAllQuality[[2]]
	asMatQualityResult[numRepeat,]<-QualityResult[[numRepeat]][[2]]
}

OverAllQuality[[2]]<-OverAllQuality[[2]]/nbRepeat # the mean probability of getting 1 for each class of prediction

## plot
# base plot
plot(OverAllQuality[[1]],OverAllQuality[[2]],xlim=c(0,1),ylim=c(0,1))
abline(a=0,b=1)
# confidence interval
if(randomRemoved){
	# simply use the standard deviation
	sdQualityResult<-rep(0,length(QualityResult[[1]][[1]]))
	for(numPoint in 1:length(QualityResult[[1]][[1]])){
		sdQualityResult[numPoint]<-sd(asMatQualityResult[,numPoint])
	}
	errbar(OverAllQuality[[1]],OverAllQuality[[2]],OverAllQuality[[2]]+sdQualityResult,OverAllQuality[[2]]-sdQualityResult,add=TRUE)
}else{
	# use a binomial confidence interval as we look at the confidence interval of
	# the probability underlying CountPos 1 among CountTot draws
	BinomAnalysis<-binom.confint(x=CountPos,n=CountTot,conf.level=0.95,methods="exact")
	errbar(OverAllQuality[[1]],BinomAnalysis$mean,BinomAnalysis$upper,BinomAnalysis$lower)
}

# making the same for the omitted data

# making the same for the omitted data, taking into account only the same block


## use estimates from visualize_sampler.r and other analyses before in full_sampler.r
## to evaluate the quality of the prediction

# in case reloading
library("LaplacesDemon")
library("binom")
# library("truncnorm")
# library("realtimeadapt")

source("morans_functions.r")
source("spam_complement.r")
source("DeltaSampling.r")
source("functions_intercept.r")

#### parameters
max.nb.simul<-100000
PercentHousesRemoved<-5 # % of houses removed
nbRepeat<-20 # nb of test of PercentHousesRemoved houses
randomRemoved<-FALSE # resample PercentHousesRemoved or avoid previously removed (then limit nbRepeat to 100/PercentHousesRemoved)
same.map<-FALSE
subsetFromNew<-FALSE
perBlockSample<-FALSE # Nota: imply randomRemoved if nbRepeat>1
if(perBlockSample){
	randomRemoved<-TRUE
}

#### init the kernel_fit_space.r
## general before any simulation
if(!same.map){
	## parameters
	source("parameters_sampler.r")
	GLOBALSETPARAMETERS<-FALSE
	use.cofactors<-FALSE
	# changes in the map to use:
	period<-"nofall.2007"
	# get general parameters from file estimated.txt
	source("estimated.txt")
	Ku<-est.Ku
	Kv<-est.Kv
	Kc<-1
	K<-c(Ku,Kv,Kc);
	f<-meanf
	T<-meanT
	# technique initialisation
	source("pseudo_data_generation.r")
	Q.est<-QfromfT(dist_mat,AS,SB,f,T);
	est.mean.beta<-mean(est.beta)
	est.detection<-rep(est.mean.beta,dimension)
	use.autostop<-TRUE
	adaptOK<-TRUE
}else{
	source("estimated.txt")
	Ku<-est.Ku
	Kv<-est.Kv
	Kc<-1
	K<-c(Ku,Kv,Kc);
	c.val<-est.c.val
	f<-meanf
	T<-meanT 
	Q.est<-QfromfT(dist_mat,AS,SB,f,T);
	est.detection<-est.beta
}
source("prep_sampling.r")

#### initialisation before each simulation
## choice of the data to be analyzed
# # choice of only one block
# block<-2008 
# sel<-which(db$block_num==block)
# choice of everything

# actual subsetting
if(subsetFromNew){
	## prior map

	## new map
	newdata<-read.csv("knownSubset.csv")
	newknowns<-newdata[which(newdata$opened == 1),]
	db$X<-db$easting
	db$Y<-db$northing
	block_data<-subsetAround(db,newknowns$unicode,threshold)
}else{
	sel<-(1:length(db$status))
	block_data<-db[sel,]
}

## technical initialisation
QualityResult<-list()
known<-which(block_data$status<9)# choice only in known status
block_data$TrueStatus<-block_data$status # ref when choice n% of the houses to be set as unknown
nbHousesPerSlice<-ceiling(PercentHousesRemoved*length(known)/100)
nbRepeat<-min(ceiling(length(known)/nbHousesPerSlice),nbRepeat)
if(!randomRemoved){
	TotalToGuess<-sample(known,min(round(nbRepeat*PercentHousesRemoved*length(known)/100),length(known)))
	cat("Will remove",nbRepeat,"times",nbHousesPerSlice,"for a total guessed of",length(TotalToGuess),"among",length(known),"known\n")
}

## ploting the chosen data
	par(mfrow=c(2,2))
	block_visudata<-block_data$pos
	block_visudata[block_data$insp==0]<-0.5
	plot_reel(block_data$easting,block_data$northing,block_visudata,base=0,top=1)

# general initialisation according to the chosen dataset
# Nota: the variables don't really need to be reinitialized
dimension<-length(sel)
if(same.map){
	w<-est.w[sel]
	u<-est.u[sel]
}else{
	w<-rep(0,dimension)
	u<-rep(0,dimension)
}
bivect<-est.detection[sel]
y<-as.integer(rnorm(length(w),mean=w,sd=1)>0)
Q<-Q.est[sel,sel]
R <- makeRuv(dimension,Q,K);
cholR <- chol.spam(R,memory=list(nnzcolindices=300000));
cholQ<-chol.spam(Q)
if(use.cofactors){
c.comp<-drop(c.map[sel,]%*%c.val)
}
grid.stab<-seq(1,length(w),ceiling(length(w)/5))# values of the field tested for stability, keep 5 values
## calculus
est.yp.b.total<-rep(0,length(block_data$status))
est.w.b.total<-rep(0,length(block_data$status))
est.u.b.total<-rep(0,length(block_data$status))
est.y.b.total<-rep(0,length(block_data$status))
est.sd.w.total<-rep(0,length(block_data$status))

for(numRepeat in 1:nbRepeat){
	if(randomRemoved){
		if(perBlockSample){
			keptAsKnown<-c()
			for(idBlock in levels(as.factor(block_data$block_num))){
				KnownInBlock<-intersect(which(block_data$block_num==idBlock),known)
				nbKnownInBlock<-length(KnownInBlock)
				if(nbKnownInBlock>0){
					nbMaxToKeep<-max(round((100-PercentHousesRemoved)/100*nbKnownInBlock),1)
					linesHouses<-sample(KnownInBlock,nbMaxToKeep)
					keptAsKnown<-c(keptAsKnown,linesHouses)
					cat("keep:",nbMaxToKeep,"out of",nbKnownInBlock,"in block",idBlock,"\n")
				}
			}
			ToGuess<-not.in(keptAsKnown,dimension)
			TotalToGuess<-ToGuess

			cat("will keep as known",length(keptAsKnown),"out of",length(known),"known in",dimension,"initially, leaving",length(ToGuess),"to predict\n")
		}else{
			ToGuess<-sample(known,round(PercentHousesRemoved*length(known)/100))
		}
	}else{
		lowBound<-1+(numRepeat-1)*nbHousesPerSlice
		upBound<-min(lowBound+nbHousesPerSlice-1,length(TotalToGuess))
		ToGuess<-TotalToGuess[lowBound:upBound]
		cat("\ntesting",lowBound,"to",upBound,"of",length(TotalToGuess),"\n")
	}
	block_data$status<-block_data$TrueStatus
	block_data$status[ToGuess]<-rep(9,length(ToGuess))
plot(block_data$easting,block_data$northing,col=(block_data$block_num%%5+1))
sel<-which(block_data$status!=9)
lines(block_data$easting[sel],block_data$northing[sel],col="blue",pch=13,type="p")

	# technical inititialisation 
	starter<-1
	zposb<-which(block_data$status==1)
	znegb<-which(block_data$status==0)
	zNAb<-which( block_data$status==9)

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

	# Rprof()
	source("kernel_fit_space.r")
	# Rprof(NULL)
	## checking on one block
	# # par(mfrow=c(2,1))
	# # plot_reel(block_data$easting,block_data$northing,block_visudata,base=0,top=1)
	# u_pred.b<-pnorm(est.u.b,0,1)
	# plot_reel(db$easting[sel],db$northing[sel],pnorm(c.comp,0,1))
	# text(block_data$easting,block_data$northing,labels=round(pnorm(c.comp,0,1),2),pos=4)
	# plot_reel(db$easting[sel],db$northing[sel],u_pred.b)
	# text(block_data$easting,block_data$northing,labels=round(u_pred.b,2),pos=4)

	# ## general analysis
	# par(mfrow=c(2,2))
	# hist(est.w)
	# hist(est.w[db$status==9],add=T,col=4)
	# 
	# hist(pnorm(est.w))
	# hist(pnorm(est.w[db$status==9]),add=T,col=4)
	# 
	# hist(est.w.b)
	# hist(est.w.b[db$status==9],add=T,col=4)
	# 
	# hist(pnorm(est.w))
	# hist(pnorm(est.w[db$status==9]),add=T,col=4)

	# est.detection.b<-inspector[sel,]%*%est.beta
	# plot.prob.to.observed((pnorm(est.w.b,0,1)*est.detection.b)[ToGuess],block_data$TrueStatus[ToGuess])
	# plot(est.w.b[ToGuess],(est.w[sel])[ToGuess])
	# abline(a=0,b=1)
	# hist(est.yp.b[ToGuess])
	# plot.prob.to.observed(u_pred*est.detection,visudata)
	# plot.prob.to.observed(u_pred*est.detection,visudata)
	# ## quality of prediction for to be guessed
	# QualityResult[[numRepeat]]<-plot.prob.to.observed(est.yp.b[ToGuess],block_data$TrueStatus[ToGuess],xlim=c(0,1),ylim=c(0,1))
	est.yp.b.total[ToGuess]<-est.yp.b[ToGuess]
	est.yp.b.total[ToGuess]<-est.yp.b[ToGuess]
	est.u.b.total[ToGuess]<-est.u.b[ToGuess]
	est.w.b.total[ToGuess]<-est.w.b[ToGuess]
	est.sd.w.total[ToGuess]<-est.sd.w[ToGuess]
	est.y.b.total[ToGuess]<-est.y.b[ToGuess]
}
save.image("EndPredictLoop.img")
# dump(c("QualityResult","est.yp.b.total","est.u.b.total","est.w.b.total","est.y.b.total"),file="QualityResults.r")
dump(c("est.yp.b.total","est.u.b.total","est.w.b.total","est.y.b.total"),file="QualityResults.r")

#### analysis of the results
# ## transform the initial big list in digestible chunks
# nbRepeat<-length(QualityResult)
# QualityClasses<-QualityResult[[1]][[1]]
# asMatQualityResult<-mat.or.vec(nbRepeat,length(QualityResult[[1]][[1]]))
# asMatCountPos<-mat.or.vec(nbRepeat,length(QualityResult[[1]][[1]]))
# asMatCountTot<-mat.or.vec(nbRepeat,length(QualityResult[[1]][[1]]))
# for(numRepeat in 1:nbRepeat){
# 	asMatQualityResult[numRepeat,]<-QualityResult[[numRepeat]][[2]]
# 	asMatCountPos[numRepeat,]<-QualityResult[[numRepeat]][[3]]
# 	asMatCountTot[numRepeat,]<-QualityResult[[numRepeat]][[4]]
# }
# CountPos<-apply(asMatCountPos,2,sum)
# CountTot<-apply(asMatCountTot,2,sum)
# ObservedRatePos<-CountPos/CountTot
# 
# ## plot
# # base plot
# par(mfrow=c(1,2))
# plot(QualityClasses,ObservedRatePos,xlim=c(0,1),ylim=c(0,1))
# abline(a=0,b=1)
# # confidence interval
# if(randomRemoved){
# 	# simply use the standard deviation
# 	sdQualityResult<-rep(0,length(QualityResult[[1]][[1]]))
# 	for(numPoint in 1:length(QualityResult[[1]][[1]])){
# 		sdQualityResult[numPoint]<-sd(asMatQualityResult[,numPoint],na.rm=TRUE)
# 	}
# 	errbar(QualityClasses,ObservedRatePos,ObservedRatePos+sdQualityResult,ObservedRatePos-sdQualityResult,add=TRUE)
# }else{
# 	# use a binomial confidence interval as we look at the confidence interval of
# 	# the probability underlying CountPos 1 among CountTot draws
# 	BinomAnalysis<-binom.confint(x=CountPos,n=CountTot,conf.level=0.95,methods="exact")
# 	errbar(QualityClasses,BinomAnalysis$mean,BinomAnalysis$upper,BinomAnalysis$lower,add=T)
# }

# making the same for the omitted data

# making the same for the omitted data, taking into account only the same block

#### new better method, using directly est.prob.b.pos (even if equivalent)
# estimated prob for each house 
est.prob.b.pos<-est.yp.b.total*est.detection
nbClasses<-10
meanQualityByClass<-c()
meanProbInf<-c()
CountPos2<-c()
CountTot2<-c()

source("functions_intercept.r")
for(numClass in 1:nbClasses){
	selClass<-intersect(which(est.prob.b.pos<numClass/nbClasses & est.prob.b.pos>= (numClass-1)/nbClasses),TotalToGuess)

	meanQualityByClass[numClass]<-mean(est.prob.b.pos[selClass])
	CountPos2[numClass]<-length(which(block_data$TrueStatus[selClass]==1))
	TruelyObserved<-intersect(which(block_data$TrueStatus!=9),selClass)
	CountTot2[numClass]<-length(TruelyObserved)
	meanProbInf[numClass]<-mean(block_data$TrueStatus[TruelyObserved])
}

dev.new(width=3.5,height=3.5)
op <- par(mar = c(4,4,0.5,0.5))
plot(meanQualityByClass,meanProbInf,xlim=c(0,1),ylim=c(0,1),xlab="Predicted",ylab="Observed",asp=1,xaxs="i",yaxs="i")
abline(a=0,b=1)
BinomAnalysis<-binom.confint(x=CountPos2,n=CountTot2,conf.level=0.95,methods="exact")
errbar(meanQualityByClass,BinomAnalysis$mean,BinomAnalysis$upper,BinomAnalysis$lower,add=T)

dev.new(width=4.5,height=4.5)
op <- par(mar = c(4,4,1,1))
plot_reel(block_data$easting,block_data$northing,block_visudata,base=0,top=1)

# dev.print(device=pdf,"predict_quality.pdf")

# nice plotting of the predicted probability surface
zNoNA<-c(zneg,zpos)
plot_reel(block_data$easting,block_data$northing,est.prob.b.pos,base=0,top=0.7)
block_data_visu_used<-block_data$TrueStatus
block_data_visu_used[TotalToGuess]<- 8
sel<-block_data_visu_used==1 | block_data_visu_used==0
lines(block_data$easting[sel],block_data$northing[sel],col="blue",pch=13,type="p")

out<-grid.from.kernel(block_data$easting,block_data$northing,est.prob.b.pos,Kernel,T=T,f,steps=150,tr=threshold)
dist.weight<-matrix(as.vector(as.matrix(out$z)),nrow=length(out$xs))

par(mfrow=c(1,2))
image(x=out$xs,y=out$ys,z=dist.weight,asp=1,col=heat.colors(100),xlab="x (m)",ylab="y (m)")
pmat<-persp(out$xs,out$ys,dist.weight,asp=1,zlim=c(0,1),
	col=make.col.persp(dist.weight,color.function=heat.colors),
	phi=20,theta=-30,
	border=NULL, # set cell borders: NULL->border ; NA-> no limits
	xlab="x",ylab="y",zlab="influence")

#### plotting data/prediction/quality of prediction
block_visudata[block_data$TrueStatus==9]<-0
dev.new(width=12,height=4)
op <- par(mar = c(5,5,3,3),mfrow=c(1,3),cex.lab=1.5)
plot_reel(block_data$easting,block_data$northing,block_visudata,base=0,top=1,main="Infested households")
plot_reel(block_data$easting,block_data$northing,est.prob.b.pos,base=0,top=0.7,main="Prediction of infestation")
plot(meanQualityByClass,meanProbInf,xlim=c(0,1),ylim=c(0,1),xlab="Predicted observation prevalence",ylab="Observed infestation prevalence",asp=1,xaxs="i",yaxs="i",main="Quality of prediction")
abline(a=0,b=1)
errbar(meanQualityByClass,BinomAnalysis$mean,BinomAnalysis$upper,BinomAnalysis$lower,add=T)
printdev(device=pdf,"QualityPrediction.pdf")

## plot map of the predicted quality of the prediction
dev.new(width=12,height=4)
par(mfrow=c(1,3))
plot_reel(block_data$easting,block_data$northing,block_visudata,base=0,top=1,main="Infested households")
plot_reel(block_data$easting[TotalToGuess],block_data$northing[TotalToGuess],est.prob.b.pos[TotalToGuess],base=0,top=est.mean.beta,main="Prediction of infestation")
plot_reel(block_data$easting[TotalToGuess],block_data$northing[TotalToGuess],est.sd.w.total[TotalToGuess],base=-5,top=5,main="estimated quality of the prediction")

# Credible interval for the estimate probability of finding something:
est.prob.inf.fromw<-pnorm(0,est.w.b.total,1,lower.tail=FALSE)
est.prob.inf.fromw.min<-pnorm(0,est.w.b.total-est.sd.w.total,1,lower.tail=FALSE)
est.prob.inf.fromw.max<-pnorm(0,est.w.b.total+est.sd.w.total,1,lower.tail=FALSE)

plot(est.prob.inf.fromw[TotalToGuess]~est.yp.b.total[TotalToGuess])
lines(est.prob.inf.fromw.max[TotalToGuess]~est.yp.b.total[TotalToGuess],type="p")
lines(est.prob.inf.fromw.min[TotalToGuess]~est.yp.b.total[TotalToGuess],type="p")
abline(a=0,b=1)

# geometric mean of probability of observation
predicted.rates<-est.prob.b.pos[TotalToGuess]
observed<-block_visudata[TotalToGuess]
prob.obs<-rep(-1,length(observed))
prob.obs[observed==1]<-predicted.rates[observed==1]
prob.obs[observed==0]<-(1-predicted.rates[observed==0])

geom.mean.prob.obs<-exp(mean(log(prob.obs)))

# geometric mean only for guessed as good prediction
sd.prediction<-est.sd.w.total[TotalToGuess]
good.pred<-which(sd.prediction<quantile(sd.prediction,prob=0.95))
good.predicted.rates<-predicted.rates[good.pred]
good.observed<-observed[good.pred]
good.prob.obs<-rep(-1,length(good.observed))
good.prob.obs[good.observed==1]<-good.predicted.rates[good.observed==1]
good.prob.obs[good.observed==0]<-(1-good.predicted.rates[good.observed==0])

good.geom.mean.prob.obs<-exp(mean(log(good.prob.obs)))
cat("Prob of good prediction:",geom.mean.prob.obs,"only for small variances:",good.geom.mean.prob.obs,"\n")

# % of pos catched for small fractions of the population
cut.rate<-10/100
risky<-which(predicted.rates>cut.rate)
nb.pos.caught<-sum(observed[risky])
nb.pos.total<-sum(observed)
cat("by checking",length(risky)/length(predicted.rates)*100,"% of the houses, we would catch",nb.pos.caught,"out of",nb.pos.total,"positives (",100*nb.pos.caught/nb.pos.total,"%)\n")

# R2
rsquared<-1-sum((observed-predicted.rates)^2)/sum((observed-mean(observed))^2)
# MacFadden
prob.obs.null.model<-rep(-1,length(observed))
prob.obs.null.model[observed==1]<-mean(observed)
prob.obs.null.model[observed==0]<-(1-mean(observed))

MacFadden<-sum(log(prob.obs))/sum(log(prob.obs.null.model))

# Adjusted count
num.correct<-length(which(predicted.rates>0.5 & observed==1))+length(which(predicted.rates<=0.5 & observed==0))
num.most.frequent.outcome<-max(length(which(observed==1)),length(which(observed==0)))
adjusted.count<-(num.correct-num.most.frequent.outcome)/(length(observed)-num.most.frequent.outcome)

cat("R2:",rsquared,"; MacFadden:",MacFadden,"; adjusted.count:",adjusted.count,"\n")

# R2 for the predicted vs observed by class:
not.nan<-which(!is.nan(meanProbInf))
rsquared<-1-sum((meanQualityByClass[not.nan]-meanProbInf[not.nan])^2)/sum((meanProbInf[not.nan]-mean(meanProbInf[not.nan]))^2)


# #### quality of prediction by city block
# # get the mean of infestation by block in predicted households
# for(num.block in levels(as.factor(block_data$block_num))){
# 	predicted.in.this.block<-intersect(which(block_data$block.num),TotalToGuess)
# 	mean.this.block<-mean(block_data$TrueStatus[predicted.in.this.block])
# 
# 	}

#### prediction of infestation post spraying
dev.new()
par(mfrow=c(2,3))
plot_reel(block_data$easting[zNA],block_data$northing[zNA],est.yp.b[zNA],base=0,top=1)

out<-grid.from.kernel(block_data$easting[zNA],block_data$northing[zNA],est.yp.b[zNA],Kernel,T=T,f,steps=150,tr=threshold)
dist.weight<-matrix(as.vector(as.matrix(out$z)),nrow=length(out$xs))
library(fields)
image.plot(x=out$xs,y=out$ys,z=dist.weight,asp=1,col=heat.colors(100),xlab="x (m)",ylab="y (m)")

## example of prediction at a house level
draw.infested.post.spray<-rbinom(length(est.yp.b[zNA]),1,prob=est.yp.b[zNA])
plot_reel(block_data$easting[zNA],block_data$northing[zNA],draw.infested.post.spray,base=0,top=1)

# =>1119 infested houses post spraying on Paucarpata: need to do that on cyclo2 and then integrate the 
# probability to open the door we may go down quite a bit but still 1119 households + after spraying is huge
# nevertheless if the opendoor model is write, the probability to open would be almost 100% when infested so may go much much down
# Nota: in sprayed house, we can count that even if the guys didn't find something they sprayed so the population has been shot down

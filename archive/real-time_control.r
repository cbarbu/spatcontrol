## simulate real-time control on a fix map of infestation (spraying)
# the pseudo-kernel is quite efficient and allows to quickly see how the true kernel may success
# the interesting parameters are valNeg and uOfNoNeigh
# a positive valNeg encourage search arround spots even if negative (but search more arround positive houses)
# a uOfNoNeigh>valNeg allows jumps to unexplored areas when to much negatives arround
# 
# in case reloading
library("LaplacesDemon")
library("binom")
source("morans_functions.r")
source("spam_complement.r")
source("prep_sampling.r")
source("DeltaSampling.r")
graphics.off()
set.seed(5)

#### parameters
reportRate<-5/100
housesExploredPerWave<-2 # multiplied by the number of denuncia init
nbRepeat<-150 # nb of test of PercentHousesRemoved houses
randomInit<-TRUE # resample PercentHousesRemoved or avoid previously removed (then limit nbRepeat to 100/PercentHousesRemoved)
max.nb.simul<-100
epsilon<-1
pseudo.kernel<-TRUE
# for pseudo.kernel
valNeg<- -0.5
valPos<- 1
uOfNoNeigh<-valNeg
f<-meanf
T<-meanT

#### init the kernel_fit_space.r
## general before any simulation
Ku<-mean(smallsampled$Ku)
Kv<-mean(smallsampled$Kv)
K<-c(Ku,Kv,Kc);
c.val<-est.c.val
Q.est<-QfromfT(dist_mat,AS,SB,meanf,meanT);

#### initialisation before each simulation
## choice of the data to be analyzed
# # choice of only one block
# block<-2008 
# sel<-which(data$block_num==block)
# choice of everything
sel<-(1:length(data$status))
dimension<-length(sel)

# actual subsetting
block_data<-data[sel,]
block_data$TrueStatus<-data[sel,]$status # ref when choice n% of the houses to be set as unknown
percentagePos<-length(which(block_data$TrueStatus==1))/length(which(block_data$TrueStatus!=9))

visublockdata<-block_data$TrueStatus
visublockdata[visublockdata==9]<-0.5
plot_reel(block_data$easting,block_data$northing,visublockdata,base=0)

## technical initialisation
QualityResult<-list()
positives<-which(block_data$TrueStatus==1)

known<-which(block_data$TrueStatus<9)# choice only in known status
## denuncia based
# notDenunciaPos<-sample(positives,round(length(positives)*(1-reportRate)))
nbDenuncias<-round(length(positives)*reportRate)
denuncia<-sample(positives,nbDenuncias)
## encuesta based
# denuncia<-sample((1:dimension),nbDenuncias)
nbNextWave<-nbDenuncias*housesExploredPerWave
if(!randomInit){
	nbHousesPerSlice<-round(PercentHousesRemoved*length(known)/100)
	nbRepeat<-min(ceiling(length(known)/nbHousesPerSlice),nbRepeat)
	TotalToGuess<-sample(known,min(round(nbRepeat*PercentHousesRemoved*length(known)/100),length(known)))
	cat("Will remove",nbRepeat,"times",nbHousesPerSlice,"for a total of",length(TotalToGuess))
}

numSimul<-0
for(uOfNoNeigh in seq(valNeg,1,0.5)){
	numSimul<-numSimul+1
	# uOfNoNeigh<-valNeg
block_data$status<-9
block_data$status[denuncia]<-block_data$TrueStatus[denuncia]
ToGuess<-which(block_data$status==9)
currentlyKnown<-rep(0,dimension)
currentlyKnown[denuncia]<-1

dev.new()
par(mfrow=c(1,3))
visublockdata<-block_data$status
visublockdata[visublockdata==9]<-0.5
plot_reel(block_data$easting,block_data$northing,visublockdata,base=0)

## calculus
est.yp.b.total<-rep(0,length(block_data$status))
est.w.b.total<-rep(0,length(block_data$status))
est.u.b.total<-rep(0,length(block_data$status))
est.y.b.total<-rep(0,length(block_data$status))

# initialisation for first guess
est.w.b<-rep(0,dimension)
est.u.b<-rep(0,dimension)
est.yp.b<-rep(0,dimension)
bivect<-est.detection[sel]

film.known<-mat.or.vec(nbRepeat+1,dimension)
film.guessed<-mat.or.vec(nbRepeat+1,dimension)
film.knownValues<-mat.or.vec(nbRepeat+1,dimension)
film.known[1,]<-currentlyKnown
film.knownValues[1,]<-block_data$status
film.guessed[1,]<-est.yp.b
nbPos<-rep(0,nbRepeat+1)
nbTry<-rep(0,nbRepeat+1)
nbKnown<-rep(0,nbRepeat+1)

for(numRepeat in 1:nbRepeat){
	# technical inititialisation according to the chosen dataset
	starter<-1
	w<-est.w.b
	u<-est.u.b
	zposb<-which(block_data$status==1)
	znegb<-which(block_data$status==0)
	zNAb<-which(block_data$status==9)
	y<-as.integer(rnorm(length(w),mean=w,sd=1)>0)
	y[zposb]<-1
	Q<-Q.est[sel,sel]
	R <- makeR(length(u),Q,K);
	if(!pseudo.kernel){
		cholR <- chol.spam(R,nnzcolindices=7e5)
		cholQ<-chol.spam(Q,nnzcolindices=7e5)
	}else{
		Qdiff<- SpecificMultiply.spam(1/T,Kernel(T,Dmat,f),SB);
		Qdiff<-Qdiff[sel,sel]
	}
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

	# Rprof()
	if(!pseudo.kernel){
		source("kernel_fit_space.r")
	}else{
		source("pseudo_kernel.r")
	}
	plot_reel(block_data$easting,block_data$northing,visublockdata,base=0)
	lines(block_data$easting[HighestRisk],block_data$northing[HighestRisk],type="p",col="blue")
	# Rprof(NULL)
	plot_reel(block_data$easting,block_data$northing,est.u.b,base=min(uOfNoNeigh,valNeg))

	## choose next to be unveiled
	HighestRisk<-(ToGuess[order(est.u.b[ToGuess],decreasing=TRUE)])[1:min(nbNextWave,length(ToGuess))]
	block_data$status[HighestRisk]<-block_data$TrueStatus[HighestRisk]
	currentlyKnown[HighestRisk]<-1
	visublockdata<-block_data$status
	visublockdata[visublockdata==9]<-0.5
	plot_reel(block_data$easting,block_data$northing,visublockdata,base=0)
	lines(block_data$easting[HighestRisk],block_data$northing[HighestRisk],type="p",col="blue")
	# readline()
	ratePosLast<-length(which(block_data$status[HighestRisk]==1))/length(which(block_data$status[HighestRisk]!=9))
	
	film.known[numRepeat+1,]<-currentlyKnown
	film.knownValues[numRepeat+1,]<-visublockdata
	film.guessed[numRepeat+1,]<-est.yp.b
	nbPos[numRepeat+1]<-length(which(block_data$status==1))
	nbKnown[numRepeat+1]<-length(which(block_data$status!=9))
	nbTry[numRepeat+1]<-sum(currentlyKnown)
	# cat("nbPos",nbPos,"\n")
	# cat("nbKnown",nbKnown,"\n")
	cat("%known pos",(100*nbPos/length(which(block_data$TrueStatus==1)))[numRepeat+1],"\n")
	cat("%pos in known",(100*nbPos/nbKnown)[numRepeat+1],"brut:",percentagePos,"\n")
	cat("%Tried",(100*nbTry/dimension)[numRepeat+1],"\n")

	ToGuess<-which(currentlyKnown==0)
	if(length(ToGuess)==0)
		break()
}

par(mfrow=c(1,3))
plot(nbPos/nbKnown)
ratePosFound<-nbPos/length(which(block_data$TrueStatus==1))
plot(ratePosFound)
rateTry<-nbTry/length(block_data$TrueStatus)
lines(rateTry)
plot(ratePosFound-rateTry)
maxDiffNum<-which((ratePosFound-rateTry)==max(ratePosFound-rateTry))
cat("for max maxDiffNum we found",ratePosFound[maxDiffNum]*100,"% of pos by sampling", 100*rateTry[maxDiffNum],"% of the houses\n")
saved[[numSimul]]<-list(valPos,valNeg,uOfNoNeigh,rateTry,ratePosFound)
}
dev.new()
plot(c(0,1),c(0,1),asp=1,type="n",xlim=c(0,1),ylim=c(0,1))
listValNeg<-rep(0,length(saved))
for(i in 1:length(saved)){
	lines(saved[[i]][[4]],saved[[i]][[5]],lty=i)
	listValNeg[i]<-as.character(saved[[i]][[2]])
}
legend("bottomright",listValNeg,lty=1:length(saved))
dump("saved",file="save_real-time_control.r")



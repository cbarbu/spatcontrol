### general settings
nameSimul<-"Pau_exp_epsKc0.01StreetsInspTr50_Full_Kpseudoflat_noNorm_broadfT_fprior40_pbeta2.18"
use.generated <- FALSE	# TRUE use generated data  FALSE use real data
use.v<-TRUE # use local error
use.insp<-TRUE # use inspectors (fit or use fixed coef)
fit.insp<-TRUE # fit the quality of inspectors
priorinspquality<-1 # quality of inspectors [0-1] if insp not fit but used
use.NA<-TRUE # guess what is in unknown (NA)
value.NA<- 9 # allow to change unknown to other for sensitivity analysis (normal:9 if considered non-infested 0; if considered infested 1)
use.cofactors<-TRUE;
mh.cof<-TRUE; 	# metropolis hastings vs. gibbs sampling for cofactors
use.streets<-TRUE; # use the fragmentation parameter
fit.spatstruct<-TRUE # is the spatial structure (f and T) refitted?
use.f<-TRUE # if FALSE, f should be set carefully
use.NormQ <- "FALSE"; # use of a normalized precision matrix Q ("med", "mean" or FALSE)
use.cosamplingfT<-FALSE;
use.intercept<-FALSE;  # use intercept on the spatial component (x)
intercept<-0;
use.MHx<-FALSE;
begin.MHx<-20;
use.MHK<-FALSE;
logsdKv<-0.05
logsdKu<-0.05
multipleFieldSamp<-1; # if set to 1 normal sampling, if integer bigger than one repet the sampling of u,v,c,K
num.simul.fixKu<- -1 # set to -1 if don't want to fix Ku
# for the parameters of the generated data refer to pseudo_data_generation.r
freqsave=10; # frequence of recording for thined parameters
final.run<-TRUE # if TRUE stopAdaptSampling and beginEstimate and nbsimul
		# hereafter are not used 
nbsimul <-100;
beginEstimate <- nbsimul/10;
stopAdaptSampling<-beginEstimate/2 
visu.progression<-FALSE
# use.autostop<-TRUE; # use final.run

# adaptive sampling
highAcceptRate<-0.4 # according to Gelman1996
lowAcceptRate<-0.15 # idem

s<-1

# seed<-sample(1:2^6,1)
# set.seed(seed)
set.seed(11583)
library("spam")
spam.options(cholsymmetrycheck=FALSE, safemode=c(FALSE,FALSE,FALSE)) ## to speed up things once everything is ok, check same results
powerboost() ## not sure it is usefull after spam.options

## fix sampling parameters
city<-"PAUCARPATA"; # PAUCARPATA MM or generated_map
subsetcity <- 0
period<-"fall.2007" # "", "fall.2007" or "fall.2008" for Paucarpata cyclo I dataset
default_kern="exp" # can be "gaussian" or "exp" or "cauchy" or "geometric", see prep_.r for the implementation
threshold <- 100; # max number of meters to consider neighbourgs
epsilon=0.01 # the value added to the diagonal of Q in order to allow for LLH calculations and data generation
## for priors on Ku/Kv see KuKvPriors.r
Kushape <- 0.001; Kuscale <- 1000; # weakly informative prior
Kvshape <- 0.001; Kvscale <- 1000; # weakly informative prior
Kc<-0.01; # the precision of the cofactors

fprior=40 # MUST be adapted to the expected scale of the autocorrelation
# same unit than the x,y of the points
sdfprior=NA
sdlf<-2 # log(sdfprior^2/fprior^2 +1) # sdfprior standard deviation in the normal scale
# but easier to give the standard deviation on the log scale for priors
# can be downsized if divergence problems
# if both defined sdlf wins

Tprior<-1; # expected is no barrier effect
sdlT<-2; #  Wide spectrum of possibilities for the barrier effect 

meant <- 0 # deprecated
Kt <- 0.01 # deprecated

abeta <- 2; ## (18,2) allow to have the mean at 0.9
bbeta <- 18; ## (1,1) gives flat prior

sdCoupleFact<-0.5
sdx.val<-0.07

## generation parameters
Ku.r <- 1# rgamma(n=1, shape=K.hyper[1], scale=K.hyper[2]);
use.v.gen<-TRUE;
Kv.r <- 10; # rgamma(n=1, shape=K.hyper[3], scale=K.hyper[4]);
Kc.r <- 10;
Delta.r <-0; # <- rtnorm(1,mean=muDelta,sd=sdDelta,lower=0,upper=Inf); # the supplementary distance due to streets in meters
f.r<- 22.28 # characteristic distance
# or simply permit the identification of Delta
T.r=0.25 # taux d'association accross streets over association within blocks
# need to be small enough to allow for the identication of Delta
mu.r<- -1.89
mv.r<-0
beta.r<-1
if(city=="PAUCARPATA"){
	make.map.cofactors=FALSE;
	cofs<-c("CU","PE","oanimal","I.NO","P.NO")
	nbfact.gen<-length(cofs)
}else{
	nbfact.gen<-5;
	prob.fact<-0.1;
	make.map.cofactors=TRUE;
}

# map generation
if(city=="generated_map"){
	use.map.gen<-TRUE;
	nb_houses_per_sidex=7;
	nb_houses_per_sidey=7;
	ratio_street_dist=1;
	nb_blocks_per_side=7;
	inter_houses_space=10;
}else{
	use.map.gen<-FALSE;
}
## initialisation parameters
start.on.true <- FALSE # if set on TRUE the following is not used
T <- 1;
f <- 20;
Ku<-1;
Kv<-1;
int<-0;

# Kernel, NB: specific treatment of spam is 2 to 10 times more efficient
expKernel<-function(T,matdist,f){
	if(class(matdist)=="spam"){
		K<-matdist;
		K@entries<- T*exp(-matdist@entries/f);
	}else{
	K<- T*exp(-matdist/f);
	}
	# K<- T*exp(-log(2)*matdist/f); ## normalized to have 0.5 when dist=f
	return(K);
}
gaussianKernel<-function(T,matdist,f){
	if(class(matdist)=="spam"){
		K<-matdist;
		K@entries<- T*exp(-(matdist@entries*matdist@entries)/f^2);
	}else{
		K<- T*exp(-(matdist*matdist)/f^2);
	}
	# K<- T*exp(-log(2)*(matdist*matdist)/f^2);## normalized to have 0.5 when dist=f
	return(K);
}
# system.time(for(i in 1:10){A<-gaussianKernel(1,Dmat,20)})
cauchyKernel<-function(T,matdist,f){
	if(class(matdist)=="spam"){
		K<-matdist;
		K@entries<- T/(1+(matdist@entries*matdist@entries)/f^2);
	}else{
		K<- T/(1+(matdist*matdist)/f^2);
	}
	return(K);
}
geometricKernel<-function(T,matdist,f){
	if(class(matdist)=="spam"){
		K<-matdist;
		K@entries<- T/(1+matdist@entries/f);
	}else{
		K<- T/(1+matdist/f);
	}
	return(K);
}
# # presentation of the kernels:
# xabs<-seq(0,threshold)
# plot(xabs,expKernel(1,xabs,fprior),type="l",lty=1,xlab="distance (m)",ylab="strength of the association")
# lines(xabs,gaussianKernel(1,xabs,fprior),type="l",lty=2)
# lines(xabs,cauchyKernel(1,xabs,fprior),type="l",lty=3)
# lines(xabs,geometricKernel(1,xabs,fprior),type="l",lty=4)
# limf<-qlnorm(c(0.025,0.975),1+log(fprior),sdlf)
# lines(xabs,expKernel(1,xabs,limf[1]),type="l",lty=1,col="grey")
# lines(xabs,gaussianKernel(1,xabs,limf[1]),type="l",lty=2,col="grey")
# lines(xabs,cauchyKernel(1,xabs,limf[1]),type="l",lty=3,col="grey")
# lines(xabs,geometricKernel(1,xabs,limf[1]),type="l",lty=4,col="grey")
# lines(xabs,expKernel(1,xabs,limf[2]),type="l",lty=1,col="grey")
# lines(xabs,gaussianKernel(1,xabs,limf[2]),type="l",lty=2,col="grey")
# lines(xabs,cauchyKernel(1,xabs,limf[2]),type="l",lty=3,col="grey")
# lines(xabs,geometricKernel(1,xabs,limf[2]),type="l",lty=4,col="grey")
# 
# dev.print(device=pdf,paste("kernels_comparison_f",fprior,"_prior.pdf",sep=""))
# 
# legend("topright",c("Exponential","Gaussian","Cauchy","Geometric"),lty=c(1,2,3,4))

# affectation of the kernel
if(default_kern == "exp"){
	Kernel<-expKernel;
}else if(default_kern == "gaussian"){
	Kernel<-gaussianKernel;
}else if(default_kern == "cauchy"){
	Kernel<-cauchyKernel;
}else if(default_kern == "geometric"){
	Kernel<-geometricKernel;
}else{
	stop("Error, unknown kernel:",default_kern,"Please change default_kern")
}

if(!use.streets){
	use.cosamplingfT<-FALSE
}

## not to be changed after this line ####
GLOBALSETPARAMETERS<-TRUE
use.yprime<-FALSE # 



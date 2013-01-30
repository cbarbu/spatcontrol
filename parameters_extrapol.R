################################
### Parameters
################################
#=======================
# Starting values (kept all along if not fitted)
# important for extrapol.spatautocorel()
# these values be estimated using a 
# representative sample 
#=======================
T <- 0.87#0.3  # efecto de calles (bajo de 1 es barrera)
f <- 8.55 #9 # distancia caracteristica de la relacion espacial
Ku<- 2.15 #0.75 # error espacial ( mas grande es menos el error)
Kv<- 170 # error no espacial ( mas grande es menos el error)
mu <- -0.71 # fuerza inicial del preditor espacial (negativo es probabilidad de positividad menos de 0.5)
# startCof<-c(0.68,0.47,0.21,-1.14,-0.28)
# names(startCof)<-c("CU","PE","oanimal","I.NO","P.NO")
# meanBeta<-0.696 # deteccion de los inspectores

#=======================
# Priors and quantitative parameters of the model
#=======================

## spatial component
fprior=10 # MUST be adapted to the expected scale of the autocorrelation, in the same unit than the x,y of the points
sdlf<-1 # prior deviation for f on the log scale, usually safer to set it to 1
# can be downsized if divergence problems
logsdfprop<-0.1 # initial log sd for proposal

Tprior<-1; # prior on barriers expected is no barrier effect (1)
sdlT<-2; # prior deviation on log scale (2=>wide spectrum of possibilities for the barrier effect)
logsdTprop<-0.1 # initial log sd for proposal

muPrior<- NULL #NULL # 0: prior at 0.5 proba of +, usually not realistic
	   # but ok for the streets article
	   # can be set to anything in ]-Inf;+Inf[
	   # if set to NULL will default to the mu
           # corresponding to the rate of + in the observed

Kushape <- 0.001; Kuscale <- 1000; # parameter of Ku prior

epsilon <- 0.01 # the value added to the diagonal of Q, interpreted as the precision of the spatial component for isolated houses (allow for LLH calculations and data generation). Should not be changed.

## non spatial component
# Kc<-0.01; # prior precision of the cofactors, arround 0
Kc<-1/sqrt(2.5) # gelman's prior scale for coefficients in arm/bayesglm

Kvshape <- 0.001; Kvscale <- 1000; # same for Kv

## inspectors
abeta <- 7; ## (18,2) allow to have the mean at 0.9
bbeta <- 3; ## (1,1) gives flat prior

alpha.poi <- 7
beta.poi <- 3

alpha.poni <- 7
beta.poni <- 3

priorinspquality<- 0.9 # if insp not fitted quality of inspectors (default to 0.9, not 1)

## opening priors
prior.io.mean <-2 ##prior.io.mean used as default/initial value for io
prior.io.var<-2

prior.fo.mean<-1 ##prior.fo.mean used as default/initial value for fo 
prior.fo.var<-2

#=======================
# How to sample? (technical parameters of the sampling)
#=======================
freqsave=20; # frequence of saving to text files, also thinning of "massive" data, like the maps

visu.progression<-FALSE # turn on/off the mapping at each iterations
			# of the results

# speed-up parameters for spam package
spam.options(cholsymmetrycheck=FALSE, safemode=c(FALSE,FALSE,FALSE))
powerboost() ## not sure it is usefull after spam.options

## not to be changed after this line ####
GLOBALSETPARAMETERS<-TRUE


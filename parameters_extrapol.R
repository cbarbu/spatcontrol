################################
### Parameters
################################
#=======================
# Starting values (kept all along if not fitted)
# important for extrapol.spatautocorel()
# these values be estimated using a 
# representative sample 
#=======================
T <- 1 # 0.87 # efecto de calles (bajo de 1 es barrera)
f <- 15 #8.55 # distancia caracteristica de la relacion espacial
Ku<- 0.48 #2.15 # error espacial ( mas grande es menos el error)
Kv<- 170 # error no espacial ( mas grande es menos el error)
mu <- NULL # fuerza inicial del preditor espacial (negativo es probabilidad de positividad menos de 0.5)
	   # set to NULL for starting on the prior
# startCof<-c(0.68,0.47,0.21,-1.14,-0.28)
# names(startCof)<-c("CU","PE","oanimal","I.NO","P.NO")
# meanBeta<-0.696 # deteccion de los inspectores

#=======================
# Priors and quantitative parameters of the model
#=======================

#----------------------------
## spatial component
#----------------------------
### autocorrelation parameters
fprior=9 # MUST be adapted to the expected scale of the autocorrelation, in the same unit than the x,y of the points
sdlf<-1 # prior deviation for f on the log scale, usually safer to set it to 1
# can be downsized if divergence problems
logsdfprop<-0.1 # initial log sd for proposal

Tprior<-0.3; # prior on barriers expected is no barrier effect (1)
sdlT<-1; # prior deviation on log scale (2=>wide spectrum of possibilities for the barrier effect)
logsdTprop<-0.1 # initial log sd for proposal

### positiveness prior (prior for the value of each house
muPriorObs<- "QuseAverage" 
           # Can take different types of values
           # ]-Inf,+Inf[ set the mean for mu prior 
	   #     => 0: prior at 0.5 proba of positive 
           #        usually not realistic can be problematic 
	   #        but ok for the streets article
	   # "QuseAverage" will set the mu
           #        corresponding to the rate of + in the observed, 
	   #        setting the the standard deviation according 
           #        Ku,Kv, and the initial Q matrix
	   #        For the unobserved it is modified by 
	   #        ORMuPriorNonObs hereafter defined
	   # "Qkrigging" will make a simular transform than QuseAverage
           #        but doing first a krigging of the positiveness
	   # NULL will default to "QuseAverage"

ORMuPriorNonObs <- 1 # prior odds ratio of non observed points compared to observed (1 same;]0,1[ less infested, ]1,+Inf[: more infested) 
epsilonProba <- 0.001 # in the initial krigging, minimum proba, 
# also maximum proba = 1-epsilonProba (avoid qnorm -> infinite values)

Kushape <- 0.001; Kuscale <- 1000; # parameter of Ku prior

epsilon <- 1/100 # the value added to the diagonal of Q, interpreted as the precision of the spatial component for isolated houses (allow for LLH calculations and data generation). Should not be changed.
# epsilon <- 1/sqrt(2.5) # the value added to the diagonal of Q, interpreted as the precision of the spatial component for isolated houses (allow for LLH calculations and data generation). Should not be changed.

## non spatial component
# Kc<-0.01; # prior precision of the cofactors, arround 0
Kc<-1/sqrt(2.5) # gelman's prior scale for coefficients in arm/bayesglm

Kvshape <- 0.001; Kvscale <- 1000; # same for Kv

## inspectors
abeta <- 7; ## (18,2) allow to have the mean at 0.9
bbeta <- 3; ## (1,1) gives flat prior
# la qualidad media es dada por abeta/(abeta+bbeta)

alpha.poi <- 7
beta.poi <- 3

alpha.poni <- 7
beta.poni <- 3

priorinspquality<- 0.9 # if insp not fitted quality of inspectors (default to 0.9, not 1)

## opening priors
prior.io.mean <-2 ##prior.io.mean used as default/initial value for io
prior.io.var <- 0.1

prior.fo.mean <- 2 ##prior.fo.mean used as default/initial value for fo 
prior.fo.var <- 0.1

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


################################
### Parameters
################################

#=======================
# Priors and quantitative parameters of the model
#=======================
threshold <- 50; # max number of meters to consider neighbourgs

## spatial component
fprior=40 # MUST be adapted to the expected scale of the autocorrelation, in the same unit than the x,y of the points
sdlf<-2 # prior deviation for f on the log scale
# can be downsized if divergence problems
logsdfprop<-0.1 # initial log sd for proposal

Tprior<-1; # prior on barriers expected is no barrier effect (1)
sdlT<-2; # prior deviation on log scale (2=>wide spectrum of possibilities for the barrier effect)
logsdTprop<-0.1 # initial log sd for proposal

Kushape <- 0.001; Kuscale <- 1000; # parameter of Ku prior

epsilon <- 0.01 # the value added to the diagonal of Q, interpreted as the precision of the spatial component for isolated houses (allow for LLH calculations and data generation). Should not be changed.

## non spatial component
Kc<-0.01; # prior precision of the cofactors, arround 0

Kvshape <- 0.001; Kvscale <- 1000; # same for Kv

## inspectors
abeta <- 1; ## (18,2) allow to have the mean at 0.9
bbeta <- 1; ## (1,1) gives flat prior

priorinspquality<- 0.696 # if insp not fitted quality of inspectors

#=======================
# Starting values (kept all along if not fitted)
#=======================
T <- 1.0 # 0.3119698
f <- 20
Ku<- 1
Kv<- 1

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


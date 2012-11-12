##############
# Example of use of the MCMC 
##############
# As this code is continually developped 
# request for an updated version are welcome at corentin.barbu [gmail]
# get necessary functions
source("extrapol_field.R",local=TRUE)

#==================
# Data preparation
#==================
# simulation name (for use with sec_launch.sh)
nameSimul<-"exp_real_trial"

# effect of streets assessment on Paucarpata (Barbu 2012)
db<-read.csv("OriginalDataPaucarpata.csv")
# Nota: the format of these data is not "perfect" on purpose
#       to remind important "checks" in the cleaning up of data
#       in addition, duplicates should be avoided as the behavior 
#       is undefined if there are duplicated locations

# avoid a number of miscodifications
db<-set_to(db,init=c("NULL"),final=0)

# avoid geographic unknowns
db<-db[which(!is.na(db$easting)),]

# setting the names as required by extrapol.spatautocorel()
db<-changeNameCol(db,"easting","X") # from "easting" to "X"
db<-changeNameCol(db,"northing","Y")
db<-changeNameCol(db,"infested","positive")
db<-changeNameCol(db,"open","observed")
db<-changeNameCol(db,"cityBlockNum","GroupNum")
db<-changeNameCol(db,"inspector","IdObserver")

# overview
plot(db$X,db$Y,col=db$fitSet+1,asp=1)
with(db[which(db$positive==1),],lines(X,Y,col="yellow",type="p"))
legend("bottomleft",c("fitting dataset","validation dataset","infested"),col=c("red","black","yellow"),pch=1)

#==================
## Computations
#==================
# the calculations performed by fit.spatautocorel 
# are to a large extent determined by what you feed it in "db"
# if IdObserver it will account for differences between observers
# if GroupNum it will account for groups such as city-blocks

# # extrapol infestation probabilities in non-observed using parameters 
# # in "parameters_extrapol.R"
# updated<-extrapol.spatautocorel(db=db)

# fit the spatial autocorrelation and extrapol to non-observed 
# in the testing set
set.seed(777) # to be able to reproduce the results
dbFit<-fit.spatautocorel(db=db[which(db$fitSet==1),],cofactors=c("CU","PE","oanimal","I.NO","P.NO"),nbiterations=-1,threshold=50,nocheck=TRUE)

##### visualizing maps
par(mfrow=c(2,3))
plot_reel(dbFit$X,dbFit$Y,dbFit$positive,base=0,top=1,main="Data")
plot_reel(dbFit$X,dbFit$Y,dbFit$p.i,base=0,top=1,main="Probability of being positive")
plot_reel(dbFit$X,dbFit$Y,dbFit$est.u,base=0,top=1,main="Estimated spatial component")
plot_reel(dbFit$X,dbFit$Y,dbFit$est.v+dbFit$est.c,base=0,top=1,main="Estimated non-spatial component")
plot_reel(dbFit$X,dbFit$Y,dbFit$est.obs,base=0,top=1,main="Estimated observation quality")

# #=======================
# # Generation parameters
# #=======================
# 
# ## generation parameters
# Ku.r <- 1# rgamma(n=1, shape=K.hyper[1], scale=K.hyper[2]);
# use.v.gen<-TRUE;
# Kv.r <- 10; # rgamma(n=1, shape=K.hyper[3], scale=K.hyper[4]);
# Kc.r <- 10;
# Delta.r <-0; # <- rtnorm(1,mean=muDelta,sd=sdDelta,lower=0,upper=Inf); # the supplementary distance due to streets in meters
# f.r<- 22.28 # characteristic distance
# # or simply permit the identification of Delta
# T.r=0.25 # taux d'association accross streets over association within blocks
# # need to be small enough to allow for the identication of Delta
# mu.r<- -1.89
# mv.r<-0
# beta.r<-1



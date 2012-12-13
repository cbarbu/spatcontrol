##############
# Example of use of the MCMC 
##############
# As this code is continually developped 
# request for an updated version are welcome at corentin.barbu [gmail]
# get necessary functions
source("spatcontrol/spatcontrol.R",local=TRUE,chdir=TRUE)

#==================
# Data preparation
#==================
# simulation name (for use with sec_launch.sh)
nameSimul<-"Extrapol2009-2012-FA-I-east"

# effect of streets assessment on Paucarpata (Barbu 2012)
# db<-read.csv("OriginalDataPaucarpata.csv")
db<-read.csv("merge_gps_generalRociado.csv")

# don't want CICLO II
db<-db[-which(db$CICLO==2),]

# avoid geographic unknowns
db<-db[which(!is.na(db$X)),]

#### Careful: setting the names is required by spatcontrol
### infestation
db$Iinfestada <- as.numeric(db$I_TRIAT != "0" & db$I_TRIAT != "")
db$Pinfestada <- as.numeric(db$P_TRIAT != "0" & db$P_TRIAT != "")
db$positive <- as.numeric (db$Iinfestada==1 | db$Pinfestada==1)
db$positive[which(is.na(db$positive))]<-0

### sub-set to intersting
## stay in Arequipa neighborghood
plot.id(db$X,db$Y,db$D)
with(db[which(db$positive==1),],lines(X,Y,col="yellow",type="p"))
abline(v=216000)
db<-db[which(db$X>216000),]

## limit to north/south]
db<-db[which(db$D==13 | db$D==10),]

## remove houses too far from infested houses (over 1000m)
# distance from infested houses
sel<-which(db$positive == 1)
matdist<-nearest.dist(x=cbind(db$X,db$Y),y=cbind(db$X[sel],db$Y[sel]),method="euclidian",delta=600,upper=NULL)

# within the threshold or not
minDistInf <- apply_by_row_not_null.spam(matdist,min)
rm(matdist)
plot(db$X,db$Y,asp=1,pch="+",cex=0.2)
with(db[which(!is.na(minDistInf)),],lines(X,Y,col="green",type="p",cex=0.2))
with(db[which(db$positive==1),],lines(X,Y,col="blue",type="p",cex=0.2))

db <- db[which(!is.na(minDistInf)),]

### City blocks
db<-changeNameCol(db,"polygon","GroupNum")

### Observed houses
# db<-changeNameCol(db,"open","observed")
db$observed<-0
db$observed[which(db$No_CARGAS>0)]<-1

### Inspectors ID
# db<-changeNameCol(db,"inspector","IdObserver")
# remove accents and put in upper case
db$BRIGADISTA <- toupper(iconv(db$BRIGADISTA, to="ASCII//TRANSLIT"))
# avoid trailing spaces
db$BRIGADISTA <- gsub(" *$","",db$BRIGADISTA)

# avoid first letter second familly name
db$BRIGADISTA <- gsub(" [A-Z].$","",db$BRIGADISTA)

# individual corrections
db$BRIGADISTA[db$BRIGADISTA == "CARLOS VELAVELA"] <- "CARLOS VELA VELA"
db$BRIGADISTA[db$BRIGADISTA == "CRISTIAN VELASQUEZ"] <- "CRISTHIAN VELASQUEZ"
db$BRIGADISTA[db$BRIGADISTA == "JAIME ROJAS"] <- "JAIME ROJAS VILCA"
db$BRIGADISTA[db$BRIGADISTA == "OSCAR BARRIOS ZEA"] <- "OSCAR BARRIOS"
db$BRIGADISTA[db$BRIGADISTA == "VIDAL MAYTA"] <- "VIDAL MAYTA HUANCA"
db$IdObserver <- as.numeric(factor(db$BRIGADISTA))

# overview
plot(db$X,db$Y,col=2*(db$observed)+1,asp=1,pch=".")
with(db[which(db$positive==1),],lines(X,Y,col="red",type="p",pch="."))
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

# # fit only the cofactors, with a random error term
# # No option yet to fit the intercept so this is probably not what you want to do
# dbFit<-fit.spatautocorel(db=db[which(db$fitSet==1 & db$observed==1),-which(names(db) %in% c("GroupNum","IdObserver","X","Y"))],cofactors=c("CU","PE","oanimal","I.NO","P.NO"),nbiterations=-1,threshold=50,nocheck=FALSE,kern="exp",use.v=TRUE)

# get parameters for sampling 
source("parameters_extrapol.R")

# Full

dbFit<-extrapol.spatautocorel(db=db,nbiterations=-1,cofactors=c("CU","PE","I_NO","P_NO"))
# Nota: for a reasonable fit you should better use nbiterations=500000 or -1 
# for autostopping

samples<-trace.mcmc()
estimates<-posteriors.mcmc(samples=samples,dbFit=dbFit)
summary.spatcontrol(estimates=estimates)


##### visualizing maps
par(mfrow=c(2,3))
plot_reel(dbFit$X,dbFit$Y,dbFit$positive,base=0,top=1,main="Data")
plot_reel(dbFit$X,dbFit$Y,dbFit$p.i,base=0,top=1,main="Probability of being positive")
plot_reel(dbFit$X,dbFit$Y,dbFit$est.u,base=0,top=1,main="Estimated spatial component")
plot_reel(dbFit$X,dbFit$Y,dbFit$est.v+dbFit$est.c,base=0,top=1,main="Estimated non-spatial component")
plot_reel(dbFit$X,dbFit$Y,dbFit$est.obs,base=0,top=1,main="Estimated observation quality")


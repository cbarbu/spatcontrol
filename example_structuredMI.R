source("spatcontrol.R",local=TRUE)

db<-read.csv("JitteredDataPaucarpata.csv",header=TRUE)
# avoid a number of miscodifications
db<-set_to(db,init=c("NULL"),final=0)
# avoid geographic unknowns
db<-db[which(!is.na(db$X)),]

# smaller area for a quick example (still 4000 points)
db<-db[db$fitSet==1,]

# structured spatial autocorrelation analysis
breaks<-seq(15,135,15) # breaks for the distance classes

mats.neigh<-gen.mats.neigh(breaks,db$X,db$Y,db$GroupNum) # distance matrices

sMI<-structured.moransI(mats_neigh=mats.neigh,db$positive,nb_rep_sign=9) # computation of the structured Moran's I

# NB: nb_rep_sign is the number of permutations performed to assess the significance of the IS-ID statistic, using 600 permutations would be the correct order of magnitude

plot.structured.moransI(sMI) # plotting of the correlograms and significance


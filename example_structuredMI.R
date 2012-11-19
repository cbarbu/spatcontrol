source("extrapol_field.R",local=TRUE)

db<-read.csv("OriginalDataPaucarpata.csv")
# avoid a number of miscodifications
db<-set_to(db,init=c("NULL"),final=0)
# avoid geographic unknowns
db<-db[which(!is.na(db$easting)),]

# smaller area for a quick example (still 4000 points)
db<-db[,db$fitSet==1]
breaks<-seq(15,135,15)
Rprof()
mats.neigh<-gen.mats.neigh(breaks,db$easting,db$northing,db$cityBlockNum)

sMI<-structured.moransI(mats_neigh=mats.neigh,db$infested,nb_rep_sign=30)

plot.structured.moransI(sMI)
Rprof(NULL)



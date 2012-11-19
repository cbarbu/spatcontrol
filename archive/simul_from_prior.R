source("RanalysisFunctions.R") # at least for set_to
# maps.tot<-read.csv("roc_p.i_fromIadjustedwithII.csv")

# focusing on hunter for now
maps<-maps.tot[which(maps$D==7),]

par(mfrow=c(1,3))
plot_reel(maps$X,maps$Y,maps$infested,base=0)
plot_reel(maps$X,maps$Y,maps$est_p.i_db,base=0)
plot_reel(maps$X,maps$Y,maps$est_p.i_da,base=0,top=max(maps$est_p.i_da))

# zm() here helps

## import reports
den<-read.csv("DENUNCIAS_2012Jun12.csv",sep=";")
den<-set_to(den)
den<-den[which(den$L!=0),]

den$infestedDen<- (den$ADULTOS>0 | den$NINFAS>0)
den$colonizedDen<- (den$NINFAS>0)
den$unicode_gps<-make_unicode_gps(den$UNICODE)

# get denuncias in maps
maps$den<-as.numeric(maps$unicode %in% den$unicode_gps)

# get positive denuncias in maps
maps$denPos<-as.numeric(maps$unicode %in% den$unicode_gps[which(den$infestedDen)])

# plot denuncias on top of previous
par(mfrow=c(1,3))
plot_reel(maps$X,maps$Y,maps$infested,base=0)
with(maps[which(maps$denPos==1),],lines(X,Y,type="p",pty=16,col="blue"))
plot_reel(maps$X,maps$Y,maps$est_p.i_db,base=0)
with(maps[which(maps$denPos==1),],lines(X,Y,type="p",pty=16,col="blue"))
plot_reel(maps$X,maps$Y,maps$est_p.i_da,base=0,top=max(maps$est_p.i_da))
with(maps[which(maps$denPos==1),],lines(X,Y,type="p",pty=16,col="blue"))



source("spatcontrol/spatcontrol.R",local=TRUE,chdir=TRUE)
db<-read.csv("OriginalDataPaucarpata.csv")
db<-set_to(db,init=c("NULL"),final=0)
db<-db[which(!is.na(db$easting)),]
db<-changeNameCol(db,"easting","X") # from "easting" to "X"
db<-changeNameCol(db,"northing","Y")
db<-changeNameCol(db,"infested","positive")
db<-changeNameCol(db,"open","observed")
db<-changeNameCol(db,"cityBlockNum","GroupNum")
db<-changeNameCol(db,"inspector","IdObserver")

## corresponding random map
# cofactors
Tc.val<-c(0.68,0.47,0.21,-1.14,-0.28)
names(Tc.val)<-c("CU","PE","oanimal","I.NO","P.NO")

# inspectors

# seedNum<-sample(1:10000,1)
# cat("seedNum:",seedNum)
seedNum<-4786
set.seed(seedNum)

nbObs<-length(levels(as.factor(db$IdObserver)))
obs.qual<-rbeta(nbObs,4.5,2)

names(obs.qual)<-levels(as.factor(db$IdObserver))
dbFit<-gen.map(db,mu=-1.3,Ku=0.48,Kv=169,f=9.0,T=0.3,c.val=Tc.val,obs.qual=obs.qual)
dbFit$fitSet<-db$fitSet

par(mfrow=c(1,2))
plot(db$X,db$Y,col=db$fitSet+1,asp=1)
with(db[which(db$positive==1),],lines(X,Y,col="yellow",type="p"))
legend("bottomleft",c("fitting dataset","validation dataset","infested"),col=c("red","black","yellow"),pch=1)

plot(dbFit$X,dbFit$Y,col=dbFit$fitSet+1,asp=1,main=paste("Seed:",seedNum))
with(dbFit[which(dbFit$positive==1),],lines(X,Y,col="yellow",type="p"))
legend("bottomleft",c("fitting dataset","validation dataset","infested"),col=c("red","black","yellow"),pch=1)

write.csv(dbFit,"JitteredDataPaucarpata.csv",row.names=FALSE)


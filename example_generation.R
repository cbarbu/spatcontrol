source("spatcontrol.R",local=TRUE)
db<-read.csv("JitteredDataPaucarpata.csv")

## corresponding random map
# cofactors
Tc.val<-c(0.68,0.47,0.21,-1.14,-0.28)
names(Tc.val)<-c("CU","PE","oanimal","I.NO","P.NO")

# inspectors
nbObs<-length(levels(as.factor(db$IdObserver)))
obs.qual<-rbeta(nbObs,4.5,2)
names(obs.qual)<-levels(as.factor(db$IdObserver))

# spatial map generation
dbFit<-gen.map(db,mu=-1.3,Ku=0.48,Kv=169,f=9.0,T=0.3,c.val=Tc.val,obs.qual=obs.qual)
dbFit$fitSet<-db$fitSet

# visualization
par(mfrow=c(1,2))
plot(db$X,db$Y,col=db$fitSet+1,asp=1,main="Original")
with(db[which(db$positive==1),],lines(X,Y,col="yellow",type="p"))
legend("bottomleft",c("fitting dataset","validation dataset","infested"),col=c("red","black","yellow"),pch=1)

plot(dbFit$X,dbFit$Y,col=dbFit$fitSet+1,asp=1,main="Generated")
with(dbFit[which(dbFit$positive==1),],lines(X,Y,col="yellow",type="p"))
legend("bottomleft",c("fitting dataset","validation dataset","infested"),col=c("red","black","yellow"),pch=1)



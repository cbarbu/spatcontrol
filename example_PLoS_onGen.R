source("extrapol_field.R")
nameSimul<-"FullonGen"
db<-read.csv("JitteredDataPaucarpata.csv")
# avoid a number of miscodifications
db<-set_to(db,init=c("NULL"),final=0)

# avoid geographic unknowns
db<-db[which(!is.na(db$X)),]
db$positive[which(is.na(db$positive))]<-0

plot(db$X,db$Y,col=db$fitSet+1,asp=1)
with(db[which(db$positive==1),],lines(X,Y,col="yellow",type="p"))

legend("bottomleft",c("fitting dataset","validation dataset","infested"),col=c("red","black","yellow"),pch=1)

set.seed(777) # to be able to reproduce the results

# Full
dbFit<-fit.spatautocorel(db=db[which(db$fitSet==1),],nbiterations=-1,threshold=50,nocheck=FALSE,cofactors=c("CU","PE","oanimal","I.NO","P.NO"),kern="exp",use.v=TRUE)

# post analysis
samples<-trace.mcmc()
estimates<-posteriors.mcmc(samples=samples,dbFit=dbFit)
summary.spatcontrol(estimates=estimates)




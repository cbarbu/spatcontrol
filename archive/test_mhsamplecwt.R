
source("extrapol_field.R")
graphics.off()
# set.seed(1)
niterations<-10000
# generate parameters
N<-500
Tc.val<-c(1.2,0.21,-0.09,-6,0.05)
TKv<-5
names(Tc.val)<-c("CU","PE","oanimal","I.NO","P.NO")
c.map<-mat.or.vec(N,length(Tc.val))
names(c.map)<-names(Tc.val)
for(i in 1:length(Tc.val)){
  c.map[sample(1:500,100),i]<-rep(1,100)
}
c.compt<-as.matrix(c.map)%*% Tc.val
sdv<-sqrt(1/TKv)
vt<-rnorm(N,mean=0,sd=sdv)
wt<- c.compt + vt
yt<-rnorm(N,mean=w,sd=1)
values<-as.numeric(y>0)

db<-as.data.frame(cbind(values,c.map))
names(db)<-c("positive",names(Tc.val))

# for now clearly wrong:
# dbFit<-fit.spatautocorel(db=db,cofactors=c("CU","PE","oanimal","I.NO","P.NO"),nbiterations=-1,nocheck=FALSE)

## fit using wt, ok
dev.new()
c.comp<-rep(0,N)
c.val<-rep(0,length(Tc.val))

c.vals<-mat.or.vec(niterations,length(Tc.val))
Kvs<-rep(0,niterations)

sdc.val<-0.03
nit<-1
zNA<-c()
Kc<-0.01 # precision of vals prior
Kvshape <- 0.001; Kvscale <- 1000; # same for Kv
Rprof()
while(nit<=niterations){
	c.all<-mhsamplec(c.val,c.comp,c.map,sdc.val,Kc,0*values,wt,zNA);
	c.val<-c.all[[1]]
	c.comp<-c.all[[2]]
	c.vals[nit,]<-c.val

	Kv<-sampleKv(wt,c.compt,Kvshape,Kvscale)
	Kvs[nit]<-Kv

	nit<-nit+1
}
# traces
par(mfrow=c(2,3))
for(i in 1:length(Tc.val)){
  plot(c.vals[,i])
  abline(h=Tc.val[i],col="green")
}
plot(Kvs)
abline(h=TKv,col="green")
abline(h=1/sd(vt)^2,col="blue")

# posteriors
dev.new()
par(mfrow=c(2,3))
for(i in 1:length(Tc.val)){
  get.estimate(c.vals[-(1:niterations/2),i],name=names(Tc.val)[i],leg=TRUE,true.val=Tc.val[i])
}

get.estimate(Kvs,name="Kv",leg=TRUE,true.val=TKv)

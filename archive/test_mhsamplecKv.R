
graphics.off()
source("extrapol_field.R")
# set.seed(777)
niterations<-20000
# generate parameters
N<-500
Tc.val<-c(1.2,0.21,-0.09,-6,0.05)
TKv<-100
names(Tc.val)<-c("CU","PE","oanimal","I.NO","P.NO")
c.map<-mat.or.vec(N,length(Tc.val))
names(c.map)<-names(Tc.val)
for(i in 1:length(Tc.val)){
  c.map[sample(1:500,100),i]<-rep(1,100)
}
sdc<-sqrt(1/TKv)
c.compt<-as.matrix(c.map)%*% Tc.val
vt<-rnorm(N,mean=0,sd=sdc)
wt<-c.compt + vt
yt<-rnorm(N,mean=wt,sd=1)
values<-as.numeric(yt>0)

db<-as.data.frame(cbind(values,c.map))
names(db)<-c("positive",names(Tc.val))

# for now clearly wrong:
# dbFit<-fit.spatautocorel(db=db,cofactors=c("CU","PE","oanimal","I.NO","P.NO"),nbiterations=-1,nocheck=FALSE)

sample_w_nospat<-function(y,c,Kv){
## classical formula for posterior given
## y ~ N(w,1)
## w ~ N(c,1/Kv)
  var<-1/Kv
  meanPost<-y*var/(1+var)+c/(1+var)
  sdPost<-sqrt(1/(Kv+1))
  w<-rnorm(length(w),mean=meanPost,sd=sdPost)
  return(w)
}

## fit using yt
dev.new()
c.comp<-rep(0,N)
c.val<-rep(0,length(Tc.val))

c.vals<-mat.or.vec(niterations,length(Tc.val))
ws<-mat.or.vec(niterations,length(wt))
Kvs<-rep(0,niterations)

sdc.val<-0.03
w<-rnorm(N,mean=0,sd=1)
v<-w
Kv<-1

nit<-1
zNA<-c()
Kc<-0.01 # precision of vals prior
Kvshape <- 0.001; Kvscale <- 1000; # same for Kv
Rprof()
while(nit<=niterations){
	c.all<-mhsamplec(c.val,c.comp,c.map,sdc.val,Kc,v,yt,zNA);
	c.val<-c.all[[1]]
	c.comp<-c.all[[2]]
	c.vals[nit,]<-c.val

	w <- sample_w_nospat(yt,c.comp,Kv)
	ws[nit,]<-w
	v<-w-c.comp

	Kv<-sampleKv(w,c.comp,Kvshape,Kvscale)
	Kvs[nit]<-Kv

	nit<-nit+1
}
Rprof(NULL)
# distrib fields
par(mfcol=c(2,4))
hist(wt)
hist(w)
hist(vt)
hist(v)
hist(c.compt)
hist(c.comp)
hist(yt)

# traces
dev.new()
par(mfrow=c(2,3))
for(i in 1:length(Tc.val)){
  plot(c.vals[,i],pch=".")
  abline(h=Tc.val[i],col="green")
}
plot(Kvs,pch=".")
abline(h=TKv,col="green")
abline(h=1/sd(vt)^2,col="blue")

dev.new()
par(mfrow=c(4,6))
for(i in 1:max(4*6,length(Tc.val))){
  plot(ws[,i],pch=".")
  abline(h=wt[i],col="green")
}

# posteriors
dev.new()
par(mfrow=c(2,3))
for(i in 1:length(Tc.val)){
  get.estimate(c.vals[-(1:niterations/2),i],name=names(Tc.val)[i],leg=TRUE,true.val=Tc.val[i])
}

get.estimate(Kvs,name="Kv",leg=TRUE,true.val=TKv)


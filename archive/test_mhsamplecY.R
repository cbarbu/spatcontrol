
graphics.off()
source("extrapol_field.R")
# set.seed(777)
niterations<-10000
# generate parameters
N<-500
Tc.val<-c(1.2,0.21,-0.09,-3,0.05)
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
wt<-c.compt # + vt
yt<-rnorm(N,mean=wt,sd=1)
values<-as.numeric(yt>0)

db<-as.data.frame(cbind(values,c.map))
names(db)<-c("positive",names(Tc.val))

# for now clearly wrong:
dbFit<-fit.spatautocorel(db=db,cofactors=c("CU","PE","oanimal","I.NO","P.NO"),nbiterations=-1,nocheck=FALSE)
stop()

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
zpos<-which(values==1)
zneg<-which(values==0)
bivect<-rep(1,length(w))
# Kc<-0.1 # precision of vals prior
Kc<-1/sqrt(2.5) # gelman's prior scale for coefficients in arm/bayesglm
		# a bit too restrictive given the values we have here
Kvshape <- 0.001; Kvscale <- 1000; # same for Kv
while(nit<=niterations){
	c.all<-mhsamplec(c.val,c.comp,c.map,sdc.val,Kc,0*y,y,zNA);
	c.val<-c.all[[1]]
	c.comp<-c.all[[2]]
	c.vals[nit,]<-c.val

	y<-sample_y_direct(c.comp,zpos,zneg,zNA,bivect);

	nit<-nit+1
}
### analysis using bayesglm
library(arm)
db<-as.data.frame(c.map)
names(db)<-names(Tc.val)
db$values<-values
summary(bglm<-bayesglm(values~ CU+PE+oanimal+I.NO+P.NO,data=db,family=binomial(probit),n.iter=100))

# distrib fields
par(mfcol=c(2,4))
hist(wt)
hist(w)
hist(vt)
hist(v)
hist(c.compt)
hist(c.comp)
hist(yt)
hist(y)

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
	namevar<-names(Tc.val)[i]
	get.estimate(c.vals[-(1:niterations/2),i],name=namevar,leg=TRUE,true.val=Tc.val[i])
	abline(v=bglm$coefficients[[namevar]],col="pink3")
}

get.estimate(Kvs,name="Kv",leg=TRUE,true.val=TKv)




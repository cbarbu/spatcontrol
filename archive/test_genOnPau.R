
nameSimul<-"test_fullGen"
graphics.off()
source("extrapol_field.R")
set.seed(777)
niterations<-1000
### generation parameters
N<-500
Tc.val<-c(1.2,0.21,-0.09,-3,0.05)
names(Tc.val)<-c("CU","PE","oanimal","I.NO","P.NO")
TKv<-100
Tmu<- -1
TKu=0.3
Tf=10
TT=1
kern=expKernel
obs.qual=1


### generation
# cofactors
c.map<-mat.or.vec(N,length(Tc.val))
for(i in 1:length(Tc.val)){
  c.map[sample(1:500,100),i]<-rep(1,100)
}

db<-as.data.frame(cbind(values,c.map))
names(db)<-c("positive",names(Tc.val))
# X,Y
db$X<-runif(N)*100 # spread over 100m
db$Y<-runif(N)*100 # spread over 100m
names(c.map)<-names(Tc.val)

# spatial autocorrelation and final values
simul<-gen.map(db,mu=Tmu,Ku=TKu,Kv=TKv,f=Tf,T=TT,c.val=Tc.val,obs.qual=obs.qual)
db$positive<-values<-simul$z
vt<-simul$v
yt<-simul$y
wt<-simul$w
ut<-simul$u
c.compt<-simul$c

# sdc<-sqrt(1/TKv)
# c.compt<-as.matrix(c.map)%*% Tc.val
# vt<-rnorm(N,mean=0,sd=sdc)
# wt<-c.compt + vt
# yt<-rnorm(N,mean=wt,sd=1)
# values<-as.numeric(yt>0)

### analysis using bayesglm
library(arm)
summary(bglm<-bayesglm(positive~ CU+PE+oanimal+I.NO+P.NO,data=db,family=binomial(probit),n.iter=100))

### fitting with full machinery
set.seed(777)
dbFit<-fit.spatautocorel(db=db,cofactors=c("CU","PE","oanimal","I.NO","P.NO"),nbiterations=-1,nocheck=FALSE,use.v=TRUE)
estimates<-posteriors.mcmc(dbFit=dbFit)
summary.spatcontrol(dbFit=dbFit,estimates=estimates)

### fit using basic loop
c.comp<-rep(0,N)
c.val<-rep(0,length(Tc.val))

c.vals<-mat.or.vec(niterations,length(Tc.val))
ws<-mat.or.vec(niterations,length(wt))
Kvs<-rep(0,niterations)

sdc.val<-0.02
w<-rep(0,N);
y<-w
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
set.seed(777)
# cat("rand",rnorm(1));
while(nit<=niterations){
	# cat("c.val",c.val,"mw",mean(w),"Kv",Kv,"my",mean(y),"\n")
	y<-sample_y_direct(w,zpos,zneg,zNA,bivect);
    # cat("meany:",mean(y),"sums zpos",sum(zpos),"neg",sum(zneg),"NA",sum(zNA),"b",mean(bivect))

	v<-sample_v(y-c.comp,Kv)
	# cat("mv",mean(v),"my",mean(y),"mwnotr",mean(c.comp),"Kv",Kv,"\n");

	Kv<-sampleKv(v,Kvshape,Kvscale)
	Kvs[nit]<-Kv

	c.all<-mhsamplec(c.val,c.comp,c.map,sdc.val,Kc,v,y,zNA);
	c.val<-c.all[[1]]
	c.comp<-c.all[[2]]
	c.vals[nit,]<-c.val
	# cat("mc",mean(c.all),"sdc.val",sdc.val,"mv",mean(y),"mwnoc",mean(v),"Kc",Kc,"\n");

	w<-v+c.comp
	ws[nit,]<-w

	nit<-nit+1
}
c.vals<-as.data.frame(c.vals)
names(c.vals)<-names(Tc.val)
estCof<-group.posteriors(c.vals,main="Cofactors' posteriors simple",true.vals=Tc.val)
summary.spatcontrol(estCof)
# # distrib fields
# par(mfcol=c(2,4))
# hist(wt)
# hist(w)
# hist(vt)
# hist(v)
# hist(c.compt)
# hist(c.comp)
# hist(yt)
# hist(y)
# 
# # traces
# dev.new()
# par(mfrow=c(2,3))
# for(i in 1:length(Tc.val)){
#   plot(c.vals[,i],pch=".")
#   abline(h=Tc.val[i],col="green")
# }
# plot(Kvs,pch=".")
# abline(h=TKv,col="green")
# abline(h=1/sd(vt)^2,col="blue")
# 
# dev.new()
# par(mfrow=c(4,6))
# for(i in 1:max(4*6,length(Tc.val))){
#   plot(ws[,i],pch=".")
#   abline(h=wt[i],col="green")
# }
# 
# # posteriors
# dev.new()
# par(mfrow=c(2,3))
# for(i in 1:length(Tc.val)){
# 	namevar<-names(Tc.val)[i]
# 	get.estimate(c.vals[-(1:niterations/2),i],name=namevar,leg=TRUE,true.val=Tc.val[i])
# 	abline(v=bglm$coefficients[[namevar]],col="pink3")
# }
# 
# get.estimate(Kvs,name="Kv",leg=TRUE,true.val=TKv)



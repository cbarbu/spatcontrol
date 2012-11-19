
source("extrapol_field.R")
set.seed(777)
niterations<-10000
# generate parameters
N<-500
Tc.val<-c(1.2,0.21,-0.09,-6,0.05)
names(Tc.val)<-c("CU","PE","oanimal","I.NO","P.NO")
c.map<-mat.or.vec(N,length(Tc.val))
names(c.map)<-names(Tc.val)
for(i in 1:length(Tc.val)){
  c.map[sample(1:500,100),i]<-rep(1,100)
}
sdc<-sqrt(1/TKv)
values<-as.matrix(c.map)%*% Tc.val 
c.comp<-rep(0,N)
c.val<-rep(0,length(Tc.val))

c.vals<-mat.or.vec(niterations,length(Tc.val))

## fit 
sdc.val<-0.03
nit<-1
zNA<-c()
Kc<-0.01 # precision of vals prior
while(nit<=niterations){
	c.all<-mhsamplec(c.val,c.comp,c.map,sdc.val,Kc,0*values,values,zNA);
	c.val<-c.all[[1]]
	c.comp<-c.all[[2]]
	c.vals[nit,]<-c.val
	nit<-nit+1
}
par(mfrow=c(2,3))
for(i in 1:length(Tc.val)){
  est<-get.estimate(c.vals[-(1:niterations/2),i],name=names(Tc.val)[i],leg=TRUE,true.val=Tc.val[i])

}



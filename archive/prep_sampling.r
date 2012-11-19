source("convergence_test.r")
cat("init prep_sampling db:",exists("db"),"\n")
## set the name
if(use.generated){racgen<-"gen"}else{racgen<-"true"}
if(use.yprime){racy<-"yprime"}else{racy<-"ydirect"}
if(use.v){racv<-"v"}else{racv<-"uv"}
if(use.insp){racinsp<-"Insp"}else{racinsp<-"NoInsp"}

name<-paste(nameSimul,racgen,racy,racv,racinsp,"tr",threshold,nbsimul,sep="_")

## prep db 
if(use.generated){
	zpos <- which(z.r==1);
	zneg <- which(z.r==0);
	zNA <- which(z.r==9);
}else{
	z.r<-db$status
	graphics.off()
	## initialisation specific to true db
	zpos <- which(db$status==1);
	zneg <- which(db$status==0);
	zNA <- which(db$status==9);
	dev.new()
	select<-c(zneg,zpos)
	plot_reel(db$X[select],db$Y[select],2*db$status[select]-1,main="data")
}

if(use.intercept){
	intercept<-qnorm(mean(z.r[which(z.r!=9)]),0,1)
}

# priors and sampling settings
nparam<-3
if(use.v){
	nparam=nparam+1;
}
if(use.cofactors){
	nparam=nparam+1;
}
if(use.insp){
	nparam=nparam+1;
}
full.screen()
par(mfrow=c(1,nparam))
# # f with gamma law:
# mf<-0.05
# vf<- 0.01
# shf<-(2+mf^2/vf)+sqrt(mf^2/vf * (4+mf^2/vf))
# scf<-mf/(shf-1)
# xabs<-seq(0,0.2,0.0001)
# plot(xabs,dgamma(xabs,shape=shf,scale=scf),main="f prior",xlab="f",ylab="density")
# f with log normal prior

flnparam<-meansd2meansdlognorm(mean=fprior,sdlog=sdlf,sd=sdfprior)

# mf<-log(fprior)-sdlf^2/2;
mf<-flnparam$meanlog

xabs<-seq(0.001,5*fprior,0.1)
# plot(xabs,lik.f(exp(xabs),mf,sdlf),main="f prior",xlab="f",ylab="density")
plot(xabs,lik.f(xabs,mf,sdlf,log=FALSE),main="f prior",xlab="f",ylab="density",type="l")
abline(v=c(fprior-2*sdfprior,fprior+2*sdfprior))
logsdfprop<-0.1

# T
Tlnparam<-meansd2meansdlognorm(mean=fprior,sdlog=sdlf,sd=sdfprior)
mT<-Tlnparam$meanlog
xabs<-seq(0,50,0.1)
# plot(xabs,lik.T(exp(xabs),mT,sdlT),main="T prior",xlab="T (log)",ylab="density (log)")
plot(xabs,lik.T(xabs,mT,sdlT,log=FALSE),main="T prior",xlab="T",ylab="density",type="l")
logsdTprop<-0.1

if(use.cofactors){
# c
	sdc.val<-0.1/nbfact.gen;
}

# Ku
# Kushape <- 0.0001; Kuscale <- 0.0001; # approximation of flat inverse prior
# Kushape <- 2; Kuscale <- 1; # not to big variances
xabs<-seq(0.0001,10,0.1)
plot(xabs,dgamma(xabs,shape=Kushape,scale=Kuscale),main="Ku prior",xlab="Ku",ylab="Density",type="l")

if(use.v){
# Kv
# Kvshape<-2 ; Kvscale<-1 # not to big variances
# Kvshape <- 0.0001; Kvscale <- 0.0001; # approximation of flat inverse prior
xabs<-seq(0.0001,10,0.1)
plot(xabs,dgamma(xabs,shape=Kvshape,scale=Kvscale),main="Kv prior",xlab="Kv",ylab="Density",type="l")
K.hyper <- c(Kushape,Kuscale,Kvshape,Kvscale); 
}

if(use.cofactors){
  # Kc
  Kcshape<-2
  Kcscale<-1
  xabs<-seq(0.0001,10,0.1)
  plot(xabs,dgamma(xabs,shape=Kcshape,scale=Kcscale),main="Kc prior",xlab="Kc",ylab="Density",type="l")
}

# b inspectors found rate when infested
## high prior
# abeta <- 35;
# bbeta <- 4;
## uniform prior
if(use.insp){
  xabs<-seq(0,1,0.01)
  plot(xabs,dbeta(xabs,abeta,bbeta),main="Betas prior",xlab="Beta",ylab="Density",type="l")
}

priorint <- c(meant,Kt); # intercept of spatial component
k <- c(Kushape,Kuscale);  # hyperparameters of Ku
b <- c(abeta,bbeta) # hyperparameters of the detection quality
hyper <- list(t=priorint, k=k, b=b,f=c(mf,sdlf),T=c(mT,sdlT));


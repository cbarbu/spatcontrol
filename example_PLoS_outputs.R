
### 
if(visu.progression){
  printdev(device=png,file="finalmaps.pdf",width=1000,height=500)
}
### fit structure

# ("visualize_sampler.r",local=TRUE)
### traces
dev.new(file="traces.pdf")
par(mfcol=c(2,6))
plot(c(1,nbiterations),log10(c(min(sampled[,1],T.r),max(sampled[,1],T.r))),xlab="T",ylab="log10(T)",type="n")
lines(log10(sampled[,1]),type="l")
abline(h=log10(T.r))
plot(sampled[,2],xlab="LLHT",type="l")
plot(c(1,nbiterations),log10(c(min(sampled[,3],f.r),max(sampled[,3],f.r))),xlab="f",ylab="log10(f)",type="n")
lines(log10(sampled[,3]),xlab="f",type="l")
if(use.generated){
  abline(h=log10(f.r))
}
# plot(sampled[,4],xlab="LLHf",type="l")
plot(sampled[,5],xlab="Ku",type="l")
if(use.generated){
  abline(h=Ku.r)
}
plot(sampled[,7],xlab="LLH z given w",type="l")
if(use.v){
  plot(sampled[,10],xlab="Kv",type="l")
  if(use.generated){
    abline(h=Kv.r)
  }
}
if(use.cofactors){
  c.vals<-read.table("cofactors.txt",header=TRUE)
  if(use.generated){
    plot(c(1,length(c.vals[,1])),c(min(c.vals,c.val.r),max(c.vals,c.val.r)),type="n")
    for(i in 1:nbfact.gen){
      abline(h=c.val.r[i],col=i);	
    }
  }else{
    plot(c(1,nbiterations),c(min(c.vals),max(c.vals)),type="n")
  }
  for(i in 1:nbfact.gen){
    lines(c.vals[,i],col=i)
  }
  legend("topright",cofs,col=(1:nbfact.gen),lty=1)

  # plot(sampled[,12],xlab="Kc",type="l")
  # if(use.generated){
  # abline(h=Kc.r)
  # }
  c.val<-drop(t(c.vals[dim(c.vals)[1],]))
  c.comp<-c.map%*%c.val
  est.c.val<-apply(c.vals,2,mean)
  est.c.comp<-c.map%*%est.c.val

  v<-drop(w)-drop(u)-drop(c.comp)
}
if(use.insp){
  betas<-read.table("betasamples.txt");
  tracemeanbeta<-apply(betas,1,mean)
  infbetas<-apply(betas,1,quantile,0.025)
  supbetas<-apply(betas,1,quantile,0.975)
  beta<-t(betas[length(betas[,1]),])
  bivect<-inspector%*%beta
  plot(tracemeanbeta,type="l",main="inspectors mean quality\n(with 5% quantiles)",ylim=c(0,1))
  lines(infbetas,lty=1)
  lines(supbetas,lty=1)
}

plot(sampled[,6],xlab="LLH z given y",type="l")
plot(sampled[,11],xlab="mean(u)",type="l")
abline(h=mu.r)
plot(sampled[,8],xlab="LLH y given w",type="l")

printdev(device=png,paste("traces.png",sep=""),width=800,height=400)

### posteriors

dev.new(file="prior_post.pdf")
par(mfrow=c(2,nparam))
# xabs<-seq(,3,0.05)
# plot(xabs,lik.f(10^(xabs),mf,sdlf,log=FALSE),main="f prior",xlab="f (log10)",ylab="density",type="l")
xabs<-seq(1,500,1)
plot(xabs,lik.f(xabs,mf,sdlf,log=FALSE),main="f prior",xlab="f",ylab="density",type="l")
xabs<-seq(0,50,0.1)
plot(xabs,lik.T(xabs,mT,sdlT,log=FALSE),main="T prior",xlab="T",ylab="density",type="l")
xabs<-seq(0,30,0.1)
plot(xabs,dgamma(xabs,shape=Kushape,scale=Kuscale),main="Ku prior",xlab="Ku",ylab="Density",type="l")
if(use.cofactors){
  plot(xabs,dgamma(xabs,shape=Kcshape,scale=Kcscale),main="Kc prior",xlab="Kc",ylab="Density",type="l")
}
if(use.v){
  plot(xabs,dgamma(xabs,shape=Kvshape,scale=Kvscale),main="Kv prior",xlab="Kv",ylab="Density",type="l")
}
if(use.insp){
  xabs<-seq(0,1,0.01)
  plot(xabs,dbeta(xabs,abeta,bbeta),main="Betas prior",xlab="Beta",ylab="Density",type="l")
}

hist(sampled[(beginEstimate):nbiterations,3],main=paste("f posterior (mean=",signif(mean(sampled[(beginEstimate):(nbiterations),3]),4),")",sep=""))
abline(v=mean(sampled[(beginEstimate):(nbiterations),3]))
if(use.generated){
  abline(v=f.r,col=4)
}
hist(sampled[(beginEstimate):(nbiterations),1],main=paste("T posterior (mean=",signif(mean(sampled[(beginEstimate):(nbiterations),1]),4),")",sep=""))
abline(v=mean(sampled[(beginEstimate):(nbiterations),1]))
if(use.generated){
  abline(v=T.r,col=4)
}
hist(sampled[(beginEstimate):(nbiterations),5],main=paste("Ku posterior (mean=",signif(mean(sampled[(beginEstimate):(nbiterations),5]),4),")",sep=""))
abline(v=mean(sampled[(beginEstimate):(nbiterations),5]))
if(use.generated){
  abline(v=Ku.r,col=4)
}
if(use.v){
  hist(sampled[(beginEstimate):(nbiterations),10],main=paste("Kv posterior (mean=",signif(mean(sampled[(beginEstimate):(nbiterations),10]),4),")",sep=""))
  abline(v=mean(sampled[(beginEstimate):(nbiterations),10]))
  if(use.generated){
    abline(v=Kv.r,col=4)
  }
}
if(use.insp){
  nbbetas<-length(betas[,1])
  meanbetas<-apply(betas[(nbbetas/2):(nbbetas),],2,mean)
  sdbetas<-apply(betas[(nbbetas/2):(nbbetas),],2,sd)
  limbetas<-apply(betas[(nbbetas/2):(nbbetas),],2,quantile,probs=c(0.025,0.975))
  if(use.generated){
    plot(meanbetas[order(beta.r)]~beta.r[order(beta.r)],ylim=c(0,1),xlab="inspector detection rate")
    lines(limbetas[1,order(beta.r)]~beta.r[order(beta.r)],type="l",pch=3)
    lines(limbetas[2,order(beta.r)]~beta.r[order(beta.r)],type="l",pch=3)
  }else{
    plot(meanbetas[order(meanbetas)],meanbetas[order(meanbetas)],ylim=c(0,1),xlab="inspector detection rate")
    lines(limbetas[1,order(meanbetas)]~meanbetas[order(meanbetas)],type="l",pch=3)
    lines(limbetas[2,order(meanbetas)]~meanbetas[order(meanbetas)],type="l",pch=3)
  }
  # hist(as.vector(as.matrix(betas[(nbbetas/2):nbbetas,])),main=paste("Betas posterior (mean=",signif(mean(meanbetas),4),")",sep=""),xlim=c(0,1))
  # abline(v=mean(meanbetas))
  # if(use.generated){
  # 	abline(v=beta.r,col=4)
  # }
}
printdev(device=pdf,paste("prior_post.pdf",sep=""))

if(use.cofactors){
  smallsampledcof<-c.vals[(beginEstimate):nbiterations,]
  smallsampledcof<-smallsampledcof[seq(1,length(smallsampledcof[,1]),freqsave),]
  smallsampledcof<-as.data.frame(smallsampledcof)
  nbcol<-floor(sqrt(nbfact.gen))
  dev.new(file="post_cofactors.pdf")
  par(mfrow=c(3,3))
  for(i in 1:nbfact.gen){
    hist(smallsampledcof[,i],main=cofs[i])
  }
  printdev(device=pdf,paste("post_cofactors.pdf",sep=""))
}

## final fields
dev.new(file="maps.pdf")
ncolvisu=3;
if(use.insp){
  ncolvisu=ncolvisu+1;
}
par(mfrow=c(2,ncolvisu))
visudata<-z.r
visudata[visudata==9]<-0
plot_reel(db$easting,db$northing,2*visudata-1,main="data")

visupseudodata<-2*generate_z(y,bivect,zNA)-1;
visupseudodata[is.na(db$status)]<-0
plot_reel(db$easting,db$northing,visupseudodata,main="generated z final")
if(use.insp){
  plot_reel(db$easting,db$northing,bivect*2-1,main="beta final")
}
plot_reel(db$easting,db$northing,y,main="y final")
plot_reel(db$easting,db$northing,w,main="w final")
plot_reel(db$easting,db$northing,w-u,main="c final")
plot_reel(db$easting,db$northing,u,main="u final")

printdev(device=pdf,paste("maps.pdf",sep=""))

# ## hist LLH
# dev.set(4)
# hist(sampled[(beginEstimate):(nbiterations),7],main="LLH z given u")
# abline(v=mean(sampled[(beginEstimate):(nbiterations),7]))
# hist(sampled[(beginEstimate):(nbiterations),4],main="LLHf (u given Q)")

### correlations
smallsampled<-sampled[(beginEstimate):nbiterations,]
# autoc.fact<-30
# smallsampled<-smallsampled[seq(1,length(smallsampled[,1]),autoc.fact),]
## select only arround a 1000 regularly sampled
spacing<-max(1,floor((nbiterations-beginEstimate)/1000))
smallsampled<-smallsampled[seq(1,length(smallsampled[,1]),spacing),]

smallsampled<-as.data.frame(smallsampled)

dev.new(file="correlations.pdf")
par(mfcol=c(3,nparam))

colnames(smallsampled)<-namesSampled

plot(smallsampled$LLH~smallsampled$f)
if(use.generated){
  abline(v=f.r,col=4)
}
plot(smallsampled$LLH~smallsampled$T)
if(use.generated){
  abline(v=Delta.r,col=4)
}
plot(smallsampled$LLH~smallsampled$Ku)
if(use.generated){
  abline(v=Ku.r,col=4)
}

plot(smallsampled$LLHfu~smallsampled$f)
if(use.generated){
  abline(v=f.r,col=4)
}
plot(smallsampled$LLHTu~smallsampled$T)
if(use.generated){
  abline(v=Delta.r,col=4)
}
plot(smallsampled$LLH~smallsampled$Ku)
if(use.generated){
  abline(v=Ku.r,col=4)
}

plot(smallsampled$T~smallsampled$f)
if(use.generated){
  abline(v=f.r,col=4)
  abline(h=T.r,col=4)
}
plot(smallsampled$T~smallsampled$Ku)
if(use.generated){
  abline(h=T.r,col=4)
  abline(v=Ku.r,col=4)
}
plot(smallsampled$f~smallsampled$Ku)
if(use.generated){
  abline(h=f.r,col=4)
  abline(v=Ku.r,col=4)
}
if(use.cofactors){
  plot(smallsampled$T~smallsampled$Kc)
  if(use.generated){
    abline(v=Kc.r,col=4)
    abline(h=T.r,col=4)
  }
  plot(smallsampled$f~smallsampled$Kc)
  if(use.generated){
    abline(h=f.r,col=4)
    abline(v=Kc.r,col=4)
  }
  plot(smallsampled$Ku~smallsampled$Kc)
  if(use.generated){
    abline(h=Ku.r,col=4)
    abline(v=Kc.r,col=4)
  }
}
if(use.v){
  plot(smallsampled$T~smallsampled$Kv)
  if(use.generated){
    abline(v=Kv.r,col=4)
    abline(h=T.r,col=4)
  }
  plot(smallsampled$f~smallsampled$Kv)
  if(use.generated){
    abline(h=f.r,col=4)
    abline(v=Kv.r,col=4)
  }
  plot(smallsampled$Ku~smallsampled$Kv)
  if(use.generated){
    abline(h=Ku.r,col=4)
    abline(v=Kv.r,col=4)
  }
}
printdev(device=png,paste("LLH_and_cov.png",sep=""),width=600,height=600)

## mean kernel
meanT<-mean(smallsampled$T)
meanf<-mean(smallsampled$f)
xabs<-seq(0,threshold);
startSample<-(beginEstimate)
minT<-quantile(sampled[startSample:nbiterations,1],0.025)
maxT<-quantile(sampled[startSample:nbiterations,1],0.975)
minf<-quantile(sampled[startSample:nbiterations,3],0.025)
maxf<-quantile(sampled[startSample:nbiterations,3],0.975)
mindistAS<-apply_by_row_not_null.spam(dist_mat*AS,min)
SBnoSelf<-SB
diag(SBnoSelf)<-0
SBnoSelf<-as.spam(SBnoSelf)
mindistSB<-apply_by_row_not_null.spam(dist_mat*SBnoSelf,min)
## difference at first neighbourg
dev.new(file="first_neigh_influence.pdf")
par(mfrow=c(1,2))
# weight of first neighbourg for estimated f and T
kernelT<-kern(meanT,mindistAS,meanf)
kernelNoT<-kern(1,mindistSB,meanf)
weightFirstNeigh<-data.frame(t(c(mean(kernelNoT,na.rm=TRUE),mean(kernelT,na.rm=TRUE))))
names(weightFirstNeigh)<-c("SB","NB")
wFN<-barplot(as.matrix(weightFirstNeigh),ylab="weight of the first neighbourg for estimated f and T",ylim=c(0,drop(as.matrix(kern(1,mean(mindistSB,na.rm=TRUE),maxf)))),main="spatial variations")
segments(wFN,c(quantile(kernelNoT,na.rm=TRUE,0.025),quantile(kernelT,0.025,na.rm=TRUE)),wFN,c(quantile(kernelNoT,0.975,na.rm=TRUE),quantile(kernelT,0.975,na.rm=TRUE)))

# weight of the mean first neighbourg
kernelT<-kern(sampled[startSample:nbiterations,1],mean(mindistAS),sampled[startSample:nbiterations,3])
weightFirstNeigh<-data.frame(t(c(as.matrix(kern(1,mean(mindistSB),meanf)),as.matrix(kern(meanT,mean(mindistAS),meanf)))))
names(weightFirstNeigh)<-c("SB","NB")
wFN<-barplot(as.matrix(weightFirstNeigh),ylab="weight for the mean distance of the first neighbourg",ylim=c(0,drop(as.matrix(kern(1,mean(mindistSB,na.rm=TRUE),maxf)))),main="variations across simulations")
segments(wFN,c(as.matrix(kern(1,mean(mindistSB,na.rm=TRUE),minf)),quantile(kernelT,0.025,na.rm=TRUE)),wFN,c(as.matrix(kern(1,mean(mindistSB,na.rm=TRUE),maxf)),quantile(kernelT,na.rm=TRUE,0.975)))

printdev(device=pdf,"first_neigh_influence.pdf")

#### first neighbor distances and kernels
par(mfrow=c(2,2))
hist(mindistSB)
hist(mindistAS)
plot(xabs,Ku*kern(1,xabs,meanf),ylim=c(0,1))
plot(xabs,Ku*kern(meanT,xabs,meanf),ylim=c(0,1))

dev.new(file="mean_kernel.pdf")
quantile(sampled[(beginEstimate):nbiterations,1],0.025)
# this gives falsely the impression that the two kernels are similar
plot(xabs,kern(1,xabs,meanf),type="l",ylim=c(0,1),col=4,main=paste("mean kernel (T: ",meanT,", f: ",meanf,") \n blue:intrablocks ; red: interblocks", seq=""))
lines(xabs,kern(1,xabs,minf),col=4,lty=2)
lines(xabs,kern(1,xabs,maxf),col=4,lty=2)
lines(xabs,kern(meanT,xabs,meanf),col=2)
lines(xabs,kern(minT,xabs,minf),col=2,lty=2)
lines(xabs,kern(maxT,xabs,maxf),col=2,lty=2)

printdev(device=pdf,"mean_kernel.pdf")

# plot of the difference between the two kernels in function of the distance
# as the min/max change with the distance, we use as reference the distance of first neighbourg across streets
dev.new(file="diff_rap_kernels.pdf")
par(mfrow=c(1,2))
diffKernel<-kern((1-sampled[startSample:nbiterations,1]),mean(mindistAS,na.rm=TRUE),sampled[startSample:nbiterations,3])
hist(diffKernel,main="difference between the kernels")
hist(sampled[startSample:nbiterations,1],main="rapport between the kernels",xlab="T")
# respectively above 0 and under 1
printdev(device=pdf,"diff_rap_kernels.pdf")

# at the end what we need is to show that within block matter more than across streets
# that is to say that the weight of everything SB is heigher than the weight of everything AS
# for all simulations
dev.new(file="general_influence_SB.pdf")

nbssampled<-length(smallsampled[,1])
GroupWeight<-mat.or.vec(nbssampled,2)
for(i in 1:nbssampled){
  cat(i," ")
  kernSB<-kern(1,SBdist,smallsampled[i,3])
  kernAS<-kern(smallsampled[i,1],ASdist,smallsampled[i,3])
  GroupWeight[i,1]<-mean(apply_by_row_not_null.spam(kernSB,sum,na.rm=TRUE),na.rm=TRUE)
  GroupWeight[i,2]<-mean(apply_by_row_not_null.spam(kernAS,sum,na.rm=TRUE),na.rm=TRUE)
}
GroupWeight[,2]<-GroupWeight[,1]+GroupWeight[,2]
GroupWeight[,1]<-GroupWeight[,1]/GroupWeight[,2]
GroupWeight[,2]<-1-GroupWeight[,1] 

# variance on the map for the estimated f and T
kernSB<-kern(1,SBdist,meanf)
kernAS<-kern(meanT,ASdist,meanf)
SBweights<-apply_by_row_not_null.spam(kernSB,sum)
ASweights<-apply_by_row_not_null.spam(kernAS,sum)
SBweights<-SBweights/(SBweights+ASweights)
ASweights<-1-SBweights

par(mfrow=c(2,2))
hist(SBweights,main="influence of the same block\nspatial variations")
hist(GroupWeight[,1],main="influence of the same block\nvariations across simulations")
spatVar<-locfit(~SBweights,xlim=c(0,1),alpha=0.5)
plot(spatVar,xlim=c(0,1),xlab="spatial disparities")
varMean<-locfit(~GroupWeight[,1],xlim=c(0,1),alpha=0.5)
plot(varMean,xlim=c(0,1),xlab="posterior of the mean")
printdev(device=pdf,"general_influence_SB.pdf")

if(use.insp){
  ## inspectors classification
  inspectors_class<-data.frame(inspectors,meanbetas)[order(meanbetas),]
  print(inspectors_class)
}

## print main results
cat(file="in_brief.txt","mean T:",meanT,"\n");
cat(file="in_brief.txt","mean f:",meanf,"\n",append=TRUE);
est.Ku<-mean(smallsampled$Ku)
cat(file="in_brief.txt","mean Ku:",est.Ku,"\n",append=TRUE);
if(use.v){
  est.Kv<-mean(smallsampled$Kv)
  cat(file="in_brief.txt","mean Kv:",est.Kv,"\n",append=TRUE);
}
est.mu<-mean(est.u)
cat(file="in_brief.txt","mean final u:",est.mu,"\n",append=TRUE);
if(use.v){
  est.mv<-mean(est.v)
  cat(file="in_brief.txt","mean final v:",est.mv,"\n",append=TRUE);
}
if(exists("est.y")){
  cat(file="in_brief.txt","mean final y:",mean(est.y),"\n",append=TRUE);
}

## save main results for generation/prediction
est.mean.beta<-mean(inspector%*%est.beta)
dump(c("meanf","meanT","est.Ku","est.Kv","est.mu","est.mv","est.mean.beta"),file="estimated.txt",append=TRUE)

p.i<-est.yp # base is the estimate by the model
# post spray corrections
pstaypositive<-0.00705 # probability of positive if sprayed positive
p.i[db$status==1]<-pstaypositive
# probability of newly positive 
pnewinfestation6months<-0.001101 
# probability of positive if not observed positive
p.i[db$status==0]<-p.i[db$status==0]*pstaypositive+pnewinfestation6months

# probability of positive if not sprayed
p.i[db$status==9]<-p.i[db$status==9]+pnewinfestation6months

## output a map 
db$observed<-db$status!=9
db$positive<-db$status==1
probamap<-with(db,cbind(X,Y,unicode_gps,GroupNum,FR_D,FR_M,FR_A,observed,positive,p.i))
write.csv(probamap,file="priormap.csv",row.names=FALSE)

# save.image(file="EndVisu.img")
# cat("end visu ok\n")

# ("getDIC.r",local=TRUE)
## DIC calculus (see Spiegelhalter2002)
if(!INTERMEDIARY){
  try(source("estimated.txt"),silent=TRUE)
}
# ("make_mean_field.r")
# make the mean of the fields to allow a generation of data according to it

getField<-function(filename){
  nblinesString<-system(paste("wc -l ",filename,sep=""),intern=TRUE)
  nblinesTable<-strsplit(nblinesString,split=" +")[[1]]
  nonVoid<-which(nzchar(nblinesTable))[1]
  nblines<-as.integer(nblinesTable[nonVoid[1]])
  cat("initially",nblines,"lines\n")

  skiped<-floor(nblines/2);

  nblinesread<-min(1000,nblines-skiped);
  # nblinesread<-nblines-skiped;
  cat("going to import",nblinesread,"lines beginning at",skiped,"...")
  fields<-scan(filename,skip=skiped,sep="\t",nlines=nblinesread); # can take some time (15s of for 750 lignes skipping 750 lines
  cat("Imported\n")
  fields<-matrix(fields,nrow=nblinesread,byrow=TRUE)

  return(fields)
}
if(exists("INTERMEDIARY")){
  if(INTERMEDIARY){
    # importOk<-try(source("estimated.txt"),silent=TRUE)
    # if(class(importOk)=="try-error"){
    us<-getField("usamples.txt")
    par(mfrow=c(1,2))
    est.u<-apply(us,2,mean)
    sd.u<-apply(us,2,sd)
    plot_reel(db$easting,db$northing,est.u,base=quantile(est.u,prob=c(0.05)),top=quantile(est.u,prob=c(0.95)),main="est.u")
    plot_reel(db$easting,db$northing,sd.u,base=0,top=quantile(sd.u,prob=c(0.95)),main="sd.u")
    rm(us)
    ws<-getField("wsamples.txt")
    est.w<-apply(ws,2,mean)
    w<-ws[dim(ws)[1],]
    rm(ws)
    est.comp<-c.map%*%est.c.val
    est.v<-est.w-est.u-est.comp
  }
  # }
}

par(mfrow=c(2,2))
plot_reel(db$easting,db$northing,visudata*2-1,main="data")
CIest.u<-as.vector(quantile(est.u,probs=c(0.025,0.975)))
plot_reel(db$easting,db$northing,est.u,main="mean u",base=CIest.u[1],top=CIest.u[2]);
plot_reel(db$easting,db$northing,est.c.comp,main="mean c",base=CIest.u[1],top=CIest.u[2])
plot_reel(db$easting,db$northing,est.w,main="mean w",base=CIest.u[1],top=CIest.u[2])
# dump("est.u",file="meanu.r")

# ## for film
# dev.new()
# for(i in 1:nblinesread){
# 	plot_reel(db$easting,db$northing,mus[i,])
# 	# Sys.sleep(1)
# }
if(INTERMEDIARY){
  # if(!exists("betas")){
  betas<-read.table("betasamples.txt");
  # }
  est.beta<-apply(betas,2,mean)
}

# we are in fact only interested by the influence of f/T but filtering out everything else
# may introduce a bad bias so we calculate with everything
Q<-QfromfT(Dmat,AS,SB,meanf,meanT,kern=kern);
meanKu<-mean(smallsampled[,5]);
# Dmean<-llh.ugivQ(dimension,est.u,Q,meanKu)+llh.zgivw(est.w,zpos,zneg,inspector%*%est.beta)+sum(dnorm(v,0,Kv,log=TRUE))
# meanD<-mean(sampled[,4]+sampled[,7]+sampled[,13]+sampled[,14]+sampled[,15]);
# ptheta<-meanD-Dmean;
# DIC <- ptheta+meanD;

# as the variables of interest here are the variables describing the spatial relationship, 
# we use a partial likelihood in the calculus of the DIC, restricted to the description of the data given the mean spatial parameters: f,T,Ku,u
# see 
# Dmean<- -2*(llh.ugivQ(dimension,est.u,Q,meanKu)+llh.zgivw(est.w,zpos,zneg,inspector%*%est.beta))
# meanD<-mean(-2*(smallsampled[,4]+smallsampled[,7]));
# if we condiser that D(theta)=-2*log(y|theta)+c, c being a constant that disappear in the comparisons
# and that y depends directly on w and beta then we can simplify the expression of D:
Dmean<- -2*(llh.zgivw(est.w,zpos,zneg,inspector%*%est.beta))
meanD<-mean(-2*(smallsampled[,7]));
ptheta<-meanD-Dmean;
DIC <- ptheta+meanD;
cat("\nDIC (partial):",DIC,"ptheta:",ptheta,"meanD:",meanD,"Dmean:",Dmean,"\n");


# cat(file="in_brief.txt","\nDIC (partial):",DIC,"ptheta:",ptheta,"meanD:",meanD,"Dmean:",Dmean,"\n",append=TRUE);
## T
par(mfcol=c(2,3))
if(use.streets){
  #pdf("Tfit.pdf",width=4,height=4)
  #par(mar=c(4,4,0.5,0.5))
  Tfit<-locfit(~sampled[(beginEstimate):(nbiterations),1],xlim=c(0,15),alpha=0.5)
  plot(Tfit,main="autocorrelation ratio accross streets\n versus inside block",xlab="T")
  xabs<-seq(0.01,15,0.01)
  lines(xabs,lik.T(xabs,mT,sdlT,log=FALSE),type="l",lty=2)
  # dev.off()
}


#pdf("ffit.pdf",width=4,height=4)
#	par(mar=c(4,4,0.5,0.5))
ffit<-locfit(~sampled[(beginEstimate):(nbiterations),3],xlim=c(1,150),alpha=1)
plot(ffit,main="kernel slope",xlab="f")
xabs<-seq(0.01,500,0.01)
lines(xabs,lik.f(xabs,mf,sdlf,log=FALSE),type="l",lty=2)
# dev.off()

Kufit<-locfit(~sampled[(beginEstimate):(nbiterations),5],xlim=c(0,150),alpha=1)
plot(Kufit,main="Spatial precision",xlab="Ku",xlim=c(0,10))
xabs<-seq(0.01,500,0.01)
lines(xabs,dgamma(xabs,shape=Kushape,scale=Kuscale),type="l",lty=2)
# lines(xabs,dgamma(xabs,shape=1.1,scale=0.5),type="l",lty=2)

Kvfit<-locfit(~sampled[(beginEstimate):(nbiterations),10],xlim=c(0,150),alpha=1)
plot(Kvfit,main="Non-Spatial precision",xlab="Kv",xlim=c(0,10))
xabs<-seq(0.01,500,0.01)
lines(xabs,dgamma(xabs,shape=Kvshape,scale=Kvscale),type="l",lty=2)
# lines(xabs,dgamma(xabs,shape=1.1,scale=0.5),type="l",lty=2)

if(use.insp){
  nbbetas<-length(betas[,1])

  plot(c(0,1),c(0,15),main="inspectors detection rates",xlab="beta",type="n")
  for(insp in 1:length(betas[1,])){
    bfit<-locfit(~betas[beginEstimate:nbbetas,insp],xlim=c(0,1),alpha=0.5)
    lines(bfit)
  }
  xabs<-seq(0,1,0.01)
  lines(xabs,dbeta(xabs,abeta,bbeta),type="l",lty=2,col=3)
}

if(use.cofactors){
  xabs<-seq(-2,2,0.01)
  plot(xabs,dnorm(xabs,0,1/sqrt(Kc)),xlab="c.val",ylab="Density",type="l",lty=2,ylim=c(0,6),main="cofactors")
  for(cof in 1:length(c.vals[1,])){
    cfit<-locfit(~c.vals[beginEstimate:nbiterations,cof],xlim=c(-15,15),alpha=1)
    lines(cfit)
  }
}
printdev(device=pdf,"prior_post_sup.pdf")



par(mfrow=c(4,2))
minx<-min(c(u,v,w))
maxx<-max(c(u,v,w))

#### all splited
## basic histogram
hist(v,xlim=c(minx,maxx),col="grey",main="risks induced for one simulation")
hist(u,add=TRUE,col="black",density=15,angle=-45)
hist(c.comp,add=TRUE,col="black",density=50,angle=-45)

hist(est.v,col="grey",xlim=c(minx,maxx),main="estimated risks induced")
hist(est.u,add=TRUE,col="black",density=15,angle=-45)
hist(c.map%*%est.c.val,add=TRUE,col="black",density=50,angle=-45)

## locfit
vfit<-locfit(~v,xlim=c(minx,maxx),alpha=0.5)
plot(vfit,main="comparison of u,v,c for specific values",xlab="probit risk",xlim=c(-5,5),lty=2,ylim=c(0,1))

ufit<-locfit(~u,xlim=c(minx,maxx),alpha=0.5)
lines(ufit,lty=1)

cfit<-locfit(~c.comp,xlim=c(minx,maxx),alpha=0.5)
lines(cfit,lty=3)

vfit<-locfit(~est.v,xlim=c(minx,maxx),alpha=0.5)
plot(vfit,main="comparison of u,v,c for specific values",xlab="probit risk",xlim=c(-5,5),lty=2,ylim=c(0,1))

ufit<-locfit(~est.u,xlim=c(minx,maxx),alpha=0.5)
lines(ufit,lty=1)

cfit<-locfit(~(c.map%*%est.c.val),xlim=c(minx,maxx),alpha=0.5)
lines(cfit,lty=3)


#### spatial versus non-spatial
## basic histogram
hist(v+c.comp,xlim=c(minx,maxx),col="grey",main="risks induced for one simulation")
hist(u,add=TRUE,col="black",density=15,angle=-45)

hist(est.v+c.map%*%est.c.val,col="grey",xlim=c(minx,maxx),main="estimated risks induced")
hist(est.u,add=TRUE,col="black",density=15,angle=-45)

## locfit
totnonspat<-(v+c.comp)
vfit<-locfit(~totnonspat,xlim=c(minx,maxx),alpha=0.5)
plot(vfit,main="comparison of u,v a specific value",xlab="probit risk",xlim=c(minx,maxx),lty=2,ylim=c(0,1))

ufit<-locfit(~u,xlim=c(minx,maxx),alpha=0.5)
lines(ufit,lty=1)

totnonspat<-(est.v+c.map%*%est.c.val)
vfit<-locfit(~totnonspat,xlim=c(minx,maxx),alpha=0.5)
plot(vfit,main="comparison of u,v estimated",xlab="probit risk",xlim=c(minx,maxx),lty=2)

ufit<-locfit(~est.u,xlim=c(minx,maxx),alpha=0.5)
lines(ufit,lty=1)

printdev(device=pdf,"comparison_risks.pdf")
plot.prob.to.observed<-function(Prob,Obs,nb_classes=10,xlab="Predicted",ylab="Observed",...){
  nb_classes<-10
  mid_class<-rep(0,nb_classes)
  prop_pos<-rep(0,nb_classes)
  tot_nb_pos<-rep(0,nb_classes)
  size_class<-rep(0,nb_classes)
  for(i in 1:nb_classes){
    min_prob<-(i-1)/nb_classes
    max_prob<-i/nb_classes
    mid_class[i]<-(min_prob+max_prob)/2
    cat("class",min_prob,"to",max_prob,":")
    risk_group<-which(Prob<max_prob & Prob>=min_prob)
    nb_pos<-sum(Obs[risk_group])
    nb_total<-length(Obs[risk_group])
    prop_pos[i]<-nb_pos/nb_total
    tot_nb_pos[i]<-nb_pos;
    cat(prop_pos[i],"of",nb_total,"\n")
    size_class[i]<-nb_total
  }
  plot(mid_class,prop_pos,xlab=xlab,ylab=ylab,...)
  abline(a=0,b=1)

  return(list(mid_class,prop_pos,tot_nb_pos,size_class))
}
par(mfrow=c(2,2))
Qmean<-QfromfT(Dmat,AS,SB,f=meanf,T=meanT,kern=kern)
spatPrec<-diag(Qmean)
u_pred<-pnorm(est.u,0,1+sqrt(1/(meanKu/spatPrec)))
plot.prob.to.observed(u_pred,visudata)
hist(u_pred)

# adjusting for inspectors
est.detection<-inspector%*%est.beta
plot.prob.to.observed(u_pred*est.detection,visudata)
hist(u_pred*est.detection)
printdev(device=pdf,"fit_by_spat.pdf")

c_pred<-pnorm(est.c.val,0,1)
w_pred<-pnorm(est.w,0,1)
final_proba<-w_pred*est.detection
par(mfrow=c(1,4))
plot_reel(db$easting, db$northing,visudata,base=0,top=1,main="data")
plot_reel(db$easting, db$northing,u_pred,base=-0,top=1,main="spatial proba")
plot_reel(db$easting, db$northing,c.map%*%c_pred,base=-0,top=1,main="cofactors proba")
plot_reel(db$easting, db$northing,final_proba,base=0,top=1,main="final proba")
printdev(device=pdf,"maps_prob.pdf")

# save.image(file="BeforePredictCheck.img")
## use estimates from visualize_sampler.r and other analyses before in full_sampler.r
## to evaluate the quality of the prediction


#### parameters
max.nb.simul<-100000
PercentHousesRemoved<-5 # % of houses removed
nbRepeat<-20 # nb of test of PercentHousesRemoved houses
randomRemoved<-FALSE # resample PercentHousesRemoved or avoid previously removed (then limit nbRepeat to 100/PercentHousesRemoved)
same.map<-FALSE
subsetFromNew<-FALSE
perBlockSample<-FALSE # Nota: imply randomRemoved if nbRepeat>1
if(perBlockSample){
  randomRemoved<-TRUE
}

#### init the kernel_fit_space.r
## general before any simulation
if(!same.map){
  ## parameters
  source("parameters_sampler.r")
  GLOBALSETPARAMETERS<-FALSE
  use.cofactors<-FALSE
  # changes in the map to use:
  period<-"nofall.2007"
  # get general parameters from file estimated.txt
  source("estimated.txt")
  Ku<-est.Ku
  Kv<-est.Kv
  Kc<-1
  K<-c(Ku,Kv,Kc);
  f<-meanf
  T<-meanT
  # technique initialisation
  # getsource("pseudo_data_generation.r")
  Q.est<-QfromfT(dist_mat,AS,SB,f,T,kern=kern);
  est.mean.beta<-mean(est.beta)
  est.detection<-rep(est.mean.beta,dimension)
  use.autostop<-TRUE
  adaptOK<-TRUE
}else{
  source("estimated.txt")
  Ku<-est.Ku
  Kv<-est.Kv
  Kc<-1
  K<-c(Ku,Kv,Kc);
  c.val<-est.c.val
  f<-meanf
  T<-meanT 
  Q.est<-QfromfT(dist_mat,AS,SB,f,T,kern=kern);
  est.detection<-est.beta
}
# getsource("prep_sampling.r")

#### initialisation before each simulation
## choice of the data to be analyzed
# # choice of only one block
# block<-2008 
# sel<-which(db$GroupNum==block)
# choice of everything

# actual subsetting
if(subsetFromNew){
  ## prior map

  ## new map
  newdata<-read.csv("knownSubset.csv")
  newknowns<-newdata[which(newdata$observed == 1),]
  db$X<-db$easting
  db$Y<-db$northing
  block_data<-subsetAround(db,newknowns$unicode,threshold)
}else{
  sel<-(1:length(db$status))
  block_data<-db[sel,]
}

## technical initialisation
QualityResult<-list()
known<-which(block_data$status<9)# choice only in known status
block_data$TrueStatus<-block_data$status # ref when choice n% of the houses to be set as unknown
nbHousesPerSlice<-ceiling(PercentHousesRemoved*length(known)/100)
nbRepeat<-min(ceiling(length(known)/nbHousesPerSlice),nbRepeat)
if(!randomRemoved){
  TotalToGuess<-sample(known,min(round(nbRepeat*PercentHousesRemoved*length(known)/100),length(known)))
  cat("Will remove",nbRepeat,"times",nbHousesPerSlice,"for a total guessed of",length(TotalToGuess),"among",length(known),"known\n")
}

## ploting the chosen data
par(mfrow=c(2,2))
block_visudata<-block_data$pos
block_visudata[block_data$insp==0]<-0.5
plot_reel(block_data$easting,block_data$northing,block_visudata,base=0,top=1)

# general initialisation according to the chosen dataset
# Nota: the variables don't really need to be reinitialized
dimension<-length(sel)
if(same.map){
  w<-est.w[sel]
  u<-est.u[sel]
}else{
  w<-rep(0,dimension)
  u<-rep(0,dimension)
}
bivect<-est.detection[sel]
y<-as.integer(rnorm(length(w),mean=w,sd=1)>0)
Q<-Q.est[sel,sel]
R <- makeRuv(dimension,Q,K);
cholR <- chol.spam(R,memory=list(nnzcolindices=300000));
cholQ<-chol.spam(Q)
if(use.cofactors){
  c.comp<-drop(c.map[sel,]%*%c.val)
}
grid.stab<-seq(1,length(w),ceiling(length(w)/5))# values of the field tested for stability, keep 5 values
## calculus
est.yp.b.total<-rep(0,length(block_data$status))
est.w.b.total<-rep(0,length(block_data$status))
est.u.b.total<-rep(0,length(block_data$status))
est.y.b.total<-rep(0,length(block_data$status))
est.sd.w.total<-rep(0,length(block_data$status))

for(numRepeat in 1:nbRepeat){
  if(randomRemoved){
    if(perBlockSample){
      keptAsKnown<-c()
      for(idBlock in levels(as.factor(block_data$GroupNum))){
	KnownInBlock<-intersect(which(block_data$GroupNum==idBlock),known)
	nbKnownInBlock<-length(KnownInBlock)
	if(nbKnownInBlock>0){
	  nbMaxToKeep<-max(round((100-PercentHousesRemoved)/100*nbKnownInBlock),1)
	  linesHouses<-sample(KnownInBlock,nbMaxToKeep)
	  keptAsKnown<-c(keptAsKnown,linesHouses)
	  cat("keep:",nbMaxToKeep,"out of",nbKnownInBlock,"in block",idBlock,"\n")
	}
      }
      ToGuess<-not.in(keptAsKnown,dimension)
      TotalToGuess<-ToGuess

      cat("will keep as known",length(keptAsKnown),"out of",length(known),"known in",dimension,"initially, leaving",length(ToGuess),"to predict\n")
    }else{
      ToGuess<-sample(known,round(PercentHousesRemoved*length(known)/100))
    }
  }else{
    lowBound<-1+(numRepeat-1)*nbHousesPerSlice
    upBound<-min(lowBound+nbHousesPerSlice-1,length(TotalToGuess))
    ToGuess<-TotalToGuess[lowBound:upBound]
    cat("\ntesting",lowBound,"to",upBound,"of",length(TotalToGuess),"\n")
  }
  block_data$status<-block_data$TrueStatus
  block_data$status[ToGuess]<-rep(9,length(ToGuess))
  plot(block_data$easting,block_data$northing,col=(block_data$GroupNum%%5+1))
  sel<-which(block_data$status!=9)
  lines(block_data$easting[sel],block_data$northing[sel],col="blue",pch=13,type="p")

  # technical inititialisation 
  starter<-1
  zposb<-which(block_data$status==1)
  znegb<-which(block_data$status==0)
  zNAb<-which( block_data$status==9)

  ItTestNum<- gibbsit(NULL,NminOnly=TRUE);
  beginEstimate<-1
  AdaptOK<-TRUE

  nbiterations<-ItTestNum+beginEstimate
  nbtraced<-2*(2+length(grid.stab))+4
  spacer<-(2+length(grid.stab))

  sampledb<-as.matrix(mat.or.vec(nbiterations+1,nbtraced));
  sampledb[1,1]<-mean(u)
  sampledb[1,2]<-sd(u)
  sampledb[1,3:spacer]<-u[grid.stab]
  sampledb[1,spacer+1]<-mean(w)
  sampledb[1,spacer+2]<-sd(u)
  sampledb[1,(spacer+3):(2*spacer)]<-w[grid.stab]
  LLHu<-llh.ugivQ(dimension,u,Q,K[1])
  sampledb[1,(2*spacer)+1]<-llh.ugivQ(dimension,u,Q,K[1])
  sampledb[1,(2*spacer)+2]<-llh.ygivw(y,w);
  sampledb[1,(2*spacer)+3]<-llh.zgivy(y,zposb,znegb,bivect);
  sampledb[1,(2*spacer)+4]<-llh.zgivw(w,zposb,znegb,bivect);

  # Rprof()
  source("kernel_fit_space.r")
  # Rprof(NULL)
  ## checking on one block
  # # par(mfrow=c(2,1))
  # # plot_reel(block_data$easting,block_data$northing,block_visudata,base=0,top=1)
  # u_pred.b<-pnorm(est.u.b,0,1)
  # plot_reel(db$easting[sel],db$northing[sel],pnorm(c.comp,0,1))
  # text(block_data$easting,block_data$northing,labels=round(pnorm(c.comp,0,1),2),pos=4)
  # plot_reel(db$easting[sel],db$northing[sel],u_pred.b)
  # text(block_data$easting,block_data$northing,labels=round(u_pred.b,2),pos=4)

  # ## general analysis
  # par(mfrow=c(2,2))
  # hist(est.w)
  # hist(est.w[db$status==9],add=T,col=4)
  # 
  # hist(pnorm(est.w))
  # hist(pnorm(est.w[db$status==9]),add=T,col=4)
  # 
  # hist(est.w.b)
  # hist(est.w.b[db$status==9],add=T,col=4)
  # 
  # hist(pnorm(est.w))
  # hist(pnorm(est.w[db$status==9]),add=T,col=4)

  # est.detection.b<-inspector[sel,]%*%est.beta
  # plot.prob.to.observed((pnorm(est.w.b,0,1)*est.detection.b)[ToGuess],block_data$TrueStatus[ToGuess])
  # plot(est.w.b[ToGuess],(est.w[sel])[ToGuess])
  # abline(a=0,b=1)
  # hist(est.yp.b[ToGuess])
  # plot.prob.to.observed(u_pred*est.detection,visudata)
  # plot.prob.to.observed(u_pred*est.detection,visudata)
  # ## quality of prediction for to be guessed
  # QualityResult[[numRepeat]]<-plot.prob.to.observed(est.yp.b[ToGuess],block_data$TrueStatus[ToGuess],xlim=c(0,1),ylim=c(0,1))
  est.yp.b.total[ToGuess]<-est.yp.b[ToGuess]
  est.yp.b.total[ToGuess]<-est.yp.b[ToGuess]
  est.u.b.total[ToGuess]<-est.u.b[ToGuess]
  est.w.b.total[ToGuess]<-est.w.b[ToGuess]
  est.sd.w.total[ToGuess]<-est.sd.w[ToGuess]
  est.y.b.total[ToGuess]<-est.y.b[ToGuess]
}
save.image("EndPredictLoop.img")
# dump(c("QualityResult","est.yp.b.total","est.u.b.total","est.w.b.total","est.y.b.total"),file="QualityResults.r")
dump(c("est.yp.b.total","est.u.b.total","est.w.b.total","est.y.b.total"),file="QualityResults.r")

#### analysis of the results
# ## transform the initial big list in digestible chunks
# nbRepeat<-length(QualityResult)
# QualityClasses<-QualityResult[[1]][[1]]
# asMatQualityResult<-mat.or.vec(nbRepeat,length(QualityResult[[1]][[1]]))
# asMatCountPos<-mat.or.vec(nbRepeat,length(QualityResult[[1]][[1]]))
# asMatCountTot<-mat.or.vec(nbRepeat,length(QualityResult[[1]][[1]]))
# for(numRepeat in 1:nbRepeat){
# 	asMatQualityResult[numRepeat,]<-QualityResult[[numRepeat]][[2]]
# 	asMatCountPos[numRepeat,]<-QualityResult[[numRepeat]][[3]]
# 	asMatCountTot[numRepeat,]<-QualityResult[[numRepeat]][[4]]
# }
# CountPos<-apply(asMatCountPos,2,sum)
# CountTot<-apply(asMatCountTot,2,sum)
# ObservedRatePos<-CountPos/CountTot
# 
# ## plot
# # base plot
# par(mfrow=c(1,2))
# plot(QualityClasses,ObservedRatePos,xlim=c(0,1),ylim=c(0,1))
# abline(a=0,b=1)
# # confidence interval
# if(randomRemoved){
# 	# simply use the standard deviation
# 	sdQualityResult<-rep(0,length(QualityResult[[1]][[1]]))
# 	for(numPoint in 1:length(QualityResult[[1]][[1]])){
# 		sdQualityResult[numPoint]<-sd(asMatQualityResult[,numPoint],na.rm=TRUE)
# 	}
# 	errbar(QualityClasses,ObservedRatePos,ObservedRatePos+sdQualityResult,ObservedRatePos-sdQualityResult,add=TRUE)
# }else{
# 	# use a binomial confidence interval as we look at the confidence interval of
# 	# the probability underlying CountPos 1 among CountTot draws
# 	BinomAnalysis<-binom.confint(x=CountPos,n=CountTot,conf.level=0.95,methods="exact")
# 	errbar(QualityClasses,BinomAnalysis$mean,BinomAnalysis$upper,BinomAnalysis$lower,add=T)
# }

# making the same for the omitted data

# making the same for the omitted data, taking into account only the same block

#### new better method, using directly est.prob.b.pos (even if equivalent)
# estimated prob for each house 
est.prob.b.pos<-est.yp.b.total*est.detection
nbClasses<-10
meanQualityByClass<-c()
meanProbInf<-c()
CountPos2<-c()
CountTot2<-c()

for(numClass in 1:nbClasses){
  selClass<-intersect(which(est.prob.b.pos<numClass/nbClasses & est.prob.b.pos>= (numClass-1)/nbClasses),TotalToGuess)

  meanQualityByClass[numClass]<-mean(est.prob.b.pos[selClass])
  CountPos2[numClass]<-length(which(block_data$TrueStatus[selClass]==1))
  TruelyObserved<-intersect(which(block_data$TrueStatus!=9),selClass)
  CountTot2[numClass]<-length(TruelyObserved)
  meanProbInf[numClass]<-mean(block_data$TrueStatus[TruelyObserved])
}

dev.new(width=3.5,height=3.5)
op <- par(mar = c(4,4,0.5,0.5))
plot(meanQualityByClass,meanProbInf,xlim=c(0,1),ylim=c(0,1),xlab="Predicted",ylab="Observed",asp=1,xaxs="i",yaxs="i")
abline(a=0,b=1)
BinomAnalysis<-binom.confint(x=CountPos2,n=CountTot2,conf.level=0.95,methods="exact")
errbar(meanQualityByClass,BinomAnalysis$mean,BinomAnalysis$upper,BinomAnalysis$lower,add=T)

dev.new(width=4.5,height=4.5)
op <- par(mar = c(4,4,1,1))
plot_reel(block_data$easting,block_data$northing,block_visudata,base=0,top=1)

# dev.print(device=pdf,"predict_quality.pdf")

# nice plotting of the predicted probability surface
zNoNA<-c(zneg,zpos)
plot_reel(block_data$easting,block_data$northing,est.prob.b.pos,base=0,top=0.7)
block_data_visu_used<-block_data$TrueStatus
block_data_visu_used[TotalToGuess]<- 8
sel<-block_data_visu_used==1 | block_data_visu_used==0
lines(block_data$easting[sel],block_data$northing[sel],col="blue",pch=13,type="p")

out<-grid.from.kernel(block_data$easting,block_data$northing,est.prob.b.pos,kern,T=T,f,steps=150,tr=threshold)
dist.weight<-matrix(as.vector(as.matrix(out$z)),nrow=length(out$xs))

par(mfrow=c(1,2))
image(x=out$xs,y=out$ys,z=dist.weight,asp=1,col=heat.colors(100),xlab="x (m)",ylab="y (m)")
pmat<-persp(out$xs,out$ys,dist.weight,asp=1,zlim=c(0,1),
	    col=make.col.persp(dist.weight,color.function=heat.colors),
	    phi=20,theta=-30,
	    border=NULL, # set cell borders: NULL->border ; NA-> no limits
	    xlab="x",ylab="y",zlab="influence")

#### plotting data/prediction/quality of prediction
block_visudata[block_data$TrueStatus==9]<-0
dev.new(width=12,height=4)
op <- par(mar = c(5,5,3,3),mfrow=c(1,3),cex.lab=1.5)
plot_reel(block_data$easting,block_data$northing,block_visudata,base=0,top=1,main="Infested households")
plot_reel(block_data$easting,block_data$northing,est.prob.b.pos,base=0,top=0.7,main="Prediction of infestation")
plot(meanQualityByClass,meanProbInf,xlim=c(0,1),ylim=c(0,1),xlab="Predicted observation prevalence",ylab="Observed infestation prevalence",asp=1,xaxs="i",yaxs="i",main="Quality of prediction")
abline(a=0,b=1)
errbar(meanQualityByClass,BinomAnalysis$mean,BinomAnalysis$upper,BinomAnalysis$lower,add=T)
printdev(device=pdf,"QualityPrediction.pdf")

## plot map of the predicted quality of the prediction
dev.new(width=12,height=4)
par(mfrow=c(1,3))
plot_reel(block_data$easting,block_data$northing,block_visudata,base=0,top=1,main="Infested households")
plot_reel(block_data$easting[TotalToGuess],block_data$northing[TotalToGuess],est.prob.b.pos[TotalToGuess],base=0,top=est.mean.beta,main="Prediction of infestation")
plot_reel(block_data$easting[TotalToGuess],block_data$northing[TotalToGuess],est.sd.w.total[TotalToGuess],base=-5,top=5,main="estimated quality of the prediction")

# Credible interval for the estimate probability of finding something:
est.prob.inf.fromw<-pnorm(0,est.w.b.total,1,lower.tail=FALSE)
est.prob.inf.fromw.min<-pnorm(0,est.w.b.total-est.sd.w.total,1,lower.tail=FALSE)
est.prob.inf.fromw.max<-pnorm(0,est.w.b.total+est.sd.w.total,1,lower.tail=FALSE)

plot(est.prob.inf.fromw[TotalToGuess]~est.yp.b.total[TotalToGuess])
lines(est.prob.inf.fromw.max[TotalToGuess]~est.yp.b.total[TotalToGuess],type="p")
lines(est.prob.inf.fromw.min[TotalToGuess]~est.yp.b.total[TotalToGuess],type="p")
abline(a=0,b=1)

# geometric mean of probability of observation
predicted.rates<-est.prob.b.pos[TotalToGuess]
observed<-block_visudata[TotalToGuess]
prob.obs<-rep(-1,length(observed))
prob.obs[observed==1]<-predicted.rates[observed==1]
prob.obs[observed==0]<-(1-predicted.rates[observed==0])

geom.mean.prob.obs<-exp(mean(log(prob.obs)))

# geometric mean only for guessed as good prediction
sd.prediction<-est.sd.w.total[TotalToGuess]
good.pred<-which(sd.prediction<quantile(sd.prediction,prob=0.95))
good.predicted.rates<-predicted.rates[good.pred]
good.observed<-observed[good.pred]
good.prob.obs<-rep(-1,length(good.observed))
good.prob.obs[good.observed==1]<-good.predicted.rates[good.observed==1]
good.prob.obs[good.observed==0]<-(1-good.predicted.rates[good.observed==0])

good.geom.mean.prob.obs<-exp(mean(log(good.prob.obs)))
cat("Prob of good prediction:",geom.mean.prob.obs,"only for small variances:",good.geom.mean.prob.obs,"\n")

# % of pos catched for small fractions of the population
cut.rate<-10/100
risky<-which(predicted.rates>cut.rate)
nb.pos.caught<-sum(observed[risky])
nb.pos.total<-sum(observed)
cat("by checking",length(risky)/length(predicted.rates)*100,"% of the houses, we would catch",nb.pos.caught,"out of",nb.pos.total,"positives (",100*nb.pos.caught/nb.pos.total,"%)\n")

# R2
rsquared<-1-sum((observed-predicted.rates)^2)/sum((observed-mean(observed))^2)
# MacFadden
prob.obs.null.model<-rep(-1,length(observed))
prob.obs.null.model[observed==1]<-mean(observed)
prob.obs.null.model[observed==0]<-(1-mean(observed))

MacFadden<-sum(log(prob.obs))/sum(log(prob.obs.null.model))

# Adjusted count
num.correct<-length(which(predicted.rates>0.5 & observed==1))+length(which(predicted.rates<=0.5 & observed==0))
num.most.frequent.outcome<-max(length(which(observed==1)),length(which(observed==0)))
adjusted.count<-(num.correct-num.most.frequent.outcome)/(length(observed)-num.most.frequent.outcome)

cat("R2:",rsquared,"; MacFadden:",MacFadden,"; adjusted.count:",adjusted.count,"\n")

# R2 for the predicted vs observed by class:
not.nan<-which(!is.nan(meanProbInf))
rsquared<-1-sum((meanQualityByClass[not.nan]-meanProbInf[not.nan])^2)/sum((meanProbInf[not.nan]-mean(meanProbInf[not.nan]))^2)


# #### quality of prediction by city block
# # get the mean of infestation by block in predicted households
# for(num.block in levels(as.factor(block_data$GroupNum))){
# 	predicted.in.this.block<-intersect(which(block_data$block.num),TotalToGuess)
# 	mean.this.block<-mean(block_data$TrueStatus[predicted.in.this.block])
# 
# 	}
#### extrapolation
dev.new()
par(mfrow=c(2,2))
if(!fit.spatstruct){
  plot(sampled[,1],main="mean u")
  plot(sampled[,2],main="sd u")
  plot(sampled[,spacer+1],main="mean w")
  plot(sampled[,spacer+2],main="sd w")
}


#### prediction of infestation post spraying
dev.new()
par(mfrow=c(2,3))
plot_reel(block_data$easting[zNA],block_data$northing[zNA],est.yp.b[zNA],base=0,top=1)

out<-grid.from.kernel(block_data$easting[zNA],block_data$northing[zNA],est.yp.b[zNA],kern,T=T,f,steps=150,tr=threshold)
dist.weight<-matrix(as.vector(as.matrix(out$z)),nrow=length(out$xs))
image.plot(x=out$xs,y=out$ys,z=dist.weight,asp=1,col=heat.colors(100),xlab="x (m)",ylab="y (m)")

## example of prediction at a house level
draw.positive<-rbinom(length(est.yp.b[zNA]),1,prob=est.yp.b[zNA])
plot_reel(block_data$easting[zNA],block_data$northing[zNA],draw.positive,base=0,top=1)

# =>1119 positive houses post spraying on Paucarpata: need to do that on cyclo2 and then integrate the 
# probability to open the door we may go down quite a bit but still 1119 households + after spraying is huge
# nevertheless if the opendoor model is write, the probability to open would be almost 100% when positive so may go much much down
# Nota: in sprayed house, we can count that even if the guys didn't find something they sprayed so the population has been shot down
# save.image(file="EndPredictCheck.img")


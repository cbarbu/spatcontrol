library("spam")
### traces
full.screen(file="traces.pdf")
par(mfcol=c(2,6))
plot(c(1,nbsimul),log10(c(min(sampled[,1],T.r),max(sampled[,1],T.r))),xlab="T",ylab="log10(T)",type="n")
lines(log10(sampled[,1]),type="l")
abline(h=log10(T.r))
plot(sampled[,2],xlab="LLHT",type="l")
plot(c(1,nbsimul),log10(c(min(sampled[,3],f.r),max(sampled[,3],f.r))),xlab="f",ylab="log10(f)",type="n")
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
	c.vals<-read.table("cofactors.txt")
	if(use.generated){
		plot(c(1,length(c.vals[,1])),c(min(c.vals,c.val.r),max(c.vals,c.val.r)),type="n")
		for(i in 1:nbfact.gen){
			abline(h=c.val.r[i],col=i);	
		}
	}else{
		plot(c(1,nbsimul),c(min(c.vals),max(c.vals)),type="n")
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
library(locfit)

full.screen(file="prior_post.pdf")
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

hist(sampled[(beginEstimate):nbsimul,3],main=paste("f posterior (mean=",signif(mean(sampled[(beginEstimate):(nbsimul),3]),4),")",sep=""))
abline(v=mean(sampled[(beginEstimate):(nbsimul),3]))
if(use.generated){
abline(v=f.r,col=4)
}
hist(sampled[(beginEstimate):(nbsimul),1],main=paste("T posterior (mean=",signif(mean(sampled[(beginEstimate):(nbsimul),1]),4),")",sep=""))
abline(v=mean(sampled[(beginEstimate):(nbsimul),1]))
if(use.generated){
abline(v=T.r,col=4)
}
hist(sampled[(beginEstimate):(nbsimul),5],main=paste("Ku posterior (mean=",signif(mean(sampled[(beginEstimate):(nbsimul),5]),4),")",sep=""))
abline(v=mean(sampled[(beginEstimate):(nbsimul),5]))
if(use.generated){
abline(v=Ku.r,col=4)
}
if(use.v){
	hist(sampled[(beginEstimate):(nbsimul),10],main=paste("Kv posterior (mean=",signif(mean(sampled[(beginEstimate):(nbsimul),10]),4),")",sep=""))
	abline(v=mean(sampled[(beginEstimate):(nbsimul),10]))
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
	smallsampledcof<-c.vals[(beginEstimate):nbsimul,]
	smallsampledcof<-smallsampledcof[seq(1,length(smallsampledcof[,1]),freqsave),]
	smallsampledcof<-as.data.frame(smallsampledcof)
	nbcol<-floor(sqrt(nbfact.gen))
	full.screen(file="post_cofactors.pdf")
	par(mfrow=c(3,3))
	for(i in 1:nbfact.gen){
		hist(smallsampledcof[,i],main=cofs[i])
	}
printdev(device=pdf,paste("post_cofactors.pdf",sep=""))
}

## final fields
full.screen(file="maps.pdf")
ncolvisu=3;
if(use.insp){
	ncolvisu=ncolvisu+1;
}
par(mfrow=c(2,ncolvisu))
visudata<-z.r
visudata[visudata==9]<-0
plot_reel(db$easting,db$northing,2*visudata-1,main="data")
if(use.NA){
dataPaleNA<-z.r
dataPaleNA[zNA]<-0.5
plot_reel(db$easting,db$northing,2*dataPaleNA-1,main="data with pale NA")
}
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
# hist(sampled[(beginEstimate):(nbsimul),7],main="LLH z given u")
# abline(v=mean(sampled[(beginEstimate):(nbsimul),7]))
# hist(sampled[(beginEstimate):(nbsimul),4],main="LLHf (u given Q)")

### correlations
smallsampled<-sampled[(beginEstimate):nbsimul,]
# autoc.fact<-30
# smallsampled<-smallsampled[seq(1,length(smallsampled[,1]),autoc.fact),]
## select only arround a 1000 regularly sampled
spacing<-max(1,floor((nbsimul-beginEstimate)/1000))
smallsampled<-smallsampled[seq(1,length(smallsampled[,1]),spacing),]

smallsampled<-as.data.frame(smallsampled)
namessmallsampled<-c("T","LLHTu","f","LLHfu","Ku","LLHy","LLH","LLHyw","i","Kv","mu","Kc","LLHv","LLHc","LLHb")

full.screen(file="correlations.pdf")
par(mfcol=c(3,nparam))

colnames(smallsampled)<-namessmallsampled

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
minT<-quantile(sampled[startSample:nbsimul,1],0.025)
maxT<-quantile(sampled[startSample:nbsimul,1],0.975)
minf<-quantile(sampled[startSample:nbsimul,3],0.025)
maxf<-quantile(sampled[startSample:nbsimul,3],0.975)
mindistAS<-apply_by_row_not_null.spam(dist_mat*AS,min)
SBnoSelf<-SB
diag(SBnoSelf)<-0
SBnoSelf<-as.spam(SBnoSelf)
mindistSB<-apply_by_row_not_null.spam(dist_mat*SBnoSelf,min)
## difference at first neighbourg
dev.new(file="first_neigh_influence.pdf")
par(mfrow=c(1,2))
# weight of first neighbourg for estimated f and T
kernelT<-Kernel(meanT,mindistAS,meanf)
kernelNoT<-Kernel(1,mindistSB,meanf)
weightFirstNeigh<-data.frame(t(c(mean(kernelNoT,na.rm=TRUE),mean(kernelT,na.rm=TRUE))))
names(weightFirstNeigh)<-c("SB","NB")
wFN<-barplot(as.matrix(weightFirstNeigh),ylab="weight of the first neighbourg for estimated f and T",ylim=c(0,drop(as.matrix(Kernel(1,mean(mindistSB,na.rm=TRUE),maxf)))),main="spatial variations")
segments(wFN,c(quantile(kernelNoT,na.rm=TRUE,0.025),quantile(kernelT,0.025,na.rm=TRUE)),wFN,c(quantile(kernelNoT,0.975,na.rm=TRUE),quantile(kernelT,0.975,na.rm=TRUE)))

# weight of the mean first neighbourg
kernelT<-Kernel(sampled[startSample:nbsimul,1],mean(mindistAS),sampled[startSample:nbsimul,3])
weightFirstNeigh<-data.frame(t(c(as.matrix(Kernel(1,mean(mindistSB),meanf)),as.matrix(Kernel(meanT,mean(mindistAS),meanf)))))
names(weightFirstNeigh)<-c("SB","NB")
wFN<-barplot(as.matrix(weightFirstNeigh),ylab="weight for the mean distance of the first neighbourg",ylim=c(0,drop(as.matrix(Kernel(1,mean(mindistSB,na.rm=TRUE),maxf)))),main="variations across simulations")
segments(wFN,c(as.matrix(Kernel(1,mean(mindistSB,na.rm=TRUE),minf)),quantile(kernelT,0.025,na.rm=TRUE)),wFN,c(as.matrix(Kernel(1,mean(mindistSB,na.rm=TRUE),maxf)),quantile(kernelT,na.rm=TRUE,0.975)))

printdev(device=pdf,"first_neigh_influence.pdf")

#### first neighbor distances and kernels
par(mfrow=c(2,2))
hist(mindistSB)
hist(mindistAS)
plot(xabs,Ku*Kernel(1,xabs,meanf),ylim=c(0,1))
plot(xabs,Ku*Kernel(meanT,xabs,meanf),ylim=c(0,1))

dev.new(file="mean_kernel.pdf")
quantile(sampled[(beginEstimate):nbsimul,1],0.025)
# this gives falsely the impression that the two kernels are similar
plot(xabs,Kernel(1,xabs,meanf),type="l",ylim=c(0,1),col=4,main=paste("mean kernel (T: ",meanT,", f: ",meanf,") \n blue:intrablocks ; red: interblocks", seq=""))
lines(xabs,Kernel(1,xabs,minf),col=4,lty=2)
lines(xabs,Kernel(1,xabs,maxf),col=4,lty=2)
lines(xabs,Kernel(meanT,xabs,meanf),col=2)
lines(xabs,Kernel(minT,xabs,minf),col=2,lty=2)
lines(xabs,Kernel(maxT,xabs,maxf),col=2,lty=2)

printdev(device=pdf,"mean_kernel.pdf")

# plot of the difference between the two kernels in function of the distance
# as the min/max change with the distance, we use as reference the distance of first neighbourg across streets
dev.new(file="diff_rap_kernels.pdf")
par(mfrow=c(1,2))
diffKernel<-Kernel((1-sampled[startSample:nbsimul,1]),mean(mindistAS,na.rm=TRUE),sampled[startSample:nbsimul,3])
hist(diffKernel,main="difference between the kernels")
hist(sampled[startSample:nbsimul,1],main="rapport between the kernels",xlab="T")
# respectively above 0 and under 1
printdev(device=pdf,"diff_rap_kernels.pdf")

# at the end what we need is to show that within block matter more than across streets
# that is to say that the weight of everything SB is heigher than the weight of everything AS
# for all simulations
dev.new(file="general_influence_SB.pdf")
ASdist<-as.spam(dist_mat*AS)
SBdist<-as.spam(dist_mat*SB)

nbssampled<-length(smallsampled[,1])
GroupWeight<-mat.or.vec(nbssampled,2)
for(i in 1:nbssampled){
	cat(i," ")
	kernSB<-Kernel(1,SBdist,smallsampled[i,3])
	kernAS<-Kernel(smallsampled[i,1],ASdist,smallsampled[i,3])
	GroupWeight[i,1]<-mean(apply_by_row_not_null.spam(kernSB,sum,na.rm=TRUE),na.rm=TRUE)
	GroupWeight[i,2]<-mean(apply_by_row_not_null.spam(kernAS,sum,na.rm=TRUE),na.rm=TRUE)
}
GroupWeight[,2]<-GroupWeight[,1]+GroupWeight[,2]
GroupWeight[,1]<-GroupWeight[,1]/GroupWeight[,2]
GroupWeight[,2]<-1-GroupWeight[,1] 

# variance on the map for the estimated f and T
kernSB<-Kernel(1,SBdist,meanf)
kernAS<-Kernel(meanT,ASdist,meanf)
SBweights<-apply_by_row_not_null.spam(kernSB,sum)
ASweights<-apply_by_row_not_null.spam(kernAS,sum)
SBweights<-SBweights/(SBweights+ASweights)
ASweights<-1-SBweights

par(mfrow=c(2,2))
hist(SBweights,main="influence of the same block\nspatial variations")
hist(GroupWeight[,1],main="influence of the same block\nvariations across simulations")
library("locfit")
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
pstayinfested<-0.00705 # probability of infested if sprayed infested
p.i[db$status==1]<-pstayinfested
# probability of newly infested 
pnewinfestation6months<-0.001101 
# probability of infested if not observed infested
p.i[db$status==0]<-p.i[db$status==0]*pstayinfested+pnewinfestation6months

# probability of infested if not sprayed
p.i[db$status==9]<-p.i[db$status==9]+pnewinfestation6months

## output a map 
db$opened<-db$status!=9
db$infested<-db$status==1
probamap<-with(db,cbind(X,Y,unicode_gps,block_num,FR_D,FR_M,FR_A,opened,infested,p.i))
write.csv(probamap,file="priormap.csv",row.names=FALSE)


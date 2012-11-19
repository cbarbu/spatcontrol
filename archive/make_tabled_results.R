load("EndSampleImage.img")

# sampling interval:
startSample<-round(length(sampled[,1])/2)

# T (lambda)
quantile(sampled[-(1:startSample),1],prob=c(0.025,0.975))
#=> give lower and upper bound of interval
mean(sampled[-(1:startSample),1])
#=> mean

# f (delta)
quantile(sampled[-(1:startSample),3],prob=c(0.025,0.975))
mean(sampled[-(1:startSample),3])
#=> mean

# Ku (tau u)
quantile(sampled[-(1:startSample),5],prob=c(0.025,0.975))
mean(sampled[-(1:startSample),5])
#=> mean

# Kv (tau v)
quantile(sampled[-(1:startSample),10],prob=c(0.025,0.975))
mean(sampled[-(1:startSample),10])
#=> mean

# mean(GroupWeight[,1])
# median(GroupWeight[,1])
# quantile(GroupWeight[,1],prob=c(0.025,0.975))


## full table from bsampled (binned by freqsave(20))
## keeping the last half of the chain (everybody may not have the same exact start
## so keep the same last half for everybody
## only DIC and SBI are taken on a smaller sample as they are much more computationally intensive

bsampled<-sampled[seq(startSample,dim(sampled)[1],freqsave),]
toBeKept<-dim(bsampled)[1]

paramsNames<-c("lambda","delta","tau_u","tau_v")
paramsNums<-c(1,3,5,10)
tabParam<-cbind(apply(bsampled[,paramsNums],2,mean),
t(apply(bsampled[,paramsNums],2,quantile,prob=c(0.025,0.975))))
rownames(tabParam)<-paramsNames

# tabParam<-rbind(tabParam,c(mean(GroupWeight[,1]),quantile(GroupWeight[,1],prob=c(0.025,0.975))))
# rownames(tabParam)[dim(tabParam)[1]]<-"SBI"

startSmallsampledcof<-dim(smallsampledcof)[1]-toBeKept
bsampledcof<-smallsampledcof[-(1:startSmallsampledcof),]
tablecof<-cbind(apply(bsampledcof,2,mean),t(apply(bsampledcof,2,quantile,prob=c(0.025,0.975))))
rownames(tablecof)<-cofs
tabParam<-rbind(tabParam,tablecof)

# mean sensitivity (beta bar)
betas$mean<-apply(betas,1,mean)
# betas$q0.975<-apply(betas,1,quantile,prob=0.975)
# betas$q0.025<-apply(betas,1,quantile,prob=0.025)
startSampleBeta<-length(betas[,1])-toBeKept
tabParam<-rbind(tabParam,c(
mean(betas[-(1:startSampleBeta),"mean"]),
quantile(betas[-(1:startSampleBeta),"mean"],prob=c(0.025,0.975))))
rownames(tabParam)[dim(tabParam)[1]]<-"mbeta"

# correctly ordered
tabParam[match(c("delta","lambda","SBI","tau_u",cofs,"tau_v","mbeta"),rownames(tabParam)),]->tabParam

roundTabParam<-signif(tabParam,4)

# plus DIC and significance of lambda
cat(paste(rownames(roundTabParam),"\t",roundTabParam[,1],"\t(",roundTabParam[,2]," - ",roundTabParam[,3],")",sep=""),
paste("DIC:",round(DIC,2),"\t lambda quantile 0.90/0.95/0.99/0.999:",paste(signif(quantile(bsampled[,1],prob=c(0.9,0.95,0.99,0.999)),3),collapse=" ")),
sep="\n")

# time per iteration
(mainLoopStop-mainLoopStart)/dim(sampled)[1]

signif(tabParam,3)


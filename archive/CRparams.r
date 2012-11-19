
nbIt<-length(sampled[,1])
first<-floor(nbIt/2)

CRmean<-function(vect,prob=c(0.025,0.975),discarded=1:first){
	value<-quantile(vect[-discarded],prob=prob)
	value["mean"]<-mean(vect[-discarded])
	value["sd"]<-sd(vect[-discarded])
	return(value)
}
CR<-c()
CR$T<-CRmean(sampled[,1],prob=c(0.05,0.95))
CR$f<-CRmean(sampled[,3])

names(c.vals)<-cofs
for(nameCof in cofs){
	CR[[nameCof]]<-CRmean(as.vector(c.vals[nameCof])[[1]],prob=c(0.05,0.95))
}

# CR of the cofactors influence
nbUsed<-1000
used<-seq(first,nbIt,round((nbIt-first)/nbUsed))
nbUsed<-length(used)
summaryCs<-mat.or.vec(nbUsed,4)
for( i in 1:length(used)){
	cs<-as.vector(c.map%*%as.vector(t(c.vals[used[i],])));
	extremesCs<-quantile(cs,prob=c(0.025,0.975))
	meanCs<-mean(cs)
	sdCs<-sd(cs)
	summaryCs[i,]<-c(extremesCs[[1]],extremesCs[[2]],meanCs,sdCs)
}
CR$meanCRc<-c(apply(summaryCs,2,mean),quantile(summaryCs[,4],prob=0.025),quantile(summaryCs[,4],prob=0.975))
names(CR$meanCRc)<c("2.5%","97.5%","mean","sd","sd 2.5%","sd 97.5%")
extremesUs<-apply(us,1,quantile,prob=c(0.025,0.975))
meanUs<-apply(us,1,mean)
sdUs<-apply(us,1,sd)
CR$meanCRu<-c(apply(extremesUs,1,mean),mean=mean(meanUs),sd=mean(sdUs),quantile(sdUs,prob=c(0.025,0.975)))

CR$est.u.sd<-sd(est.u)
CR$est.c.comp.sd<-sd(est.c.comp[,1])

CR$Beta<-c(min(est.beta),max(est.beta),mean(est.beta))
names(CR$Beta)<-c("min","max","mean")
cat("est.Beta order:",order(est.beta),"\n")

## print CR as line of grouped columns
nbfields<-c()
cat("\n",file="CR.csv")
for(name in names(CR)){
	nbfields[name]<-length(CR[[name]])
	cat(name,file="CR.csv",append=TRUE)
	for(i in 1:nbfields[name]){
		cat("\t",file="CR.csv",append=TRUE)
	}
}
cat("\n",file="CR.csv",append=TRUE)
for(name in names(CR)){
	cat(names(CR[[name]]),file="CR.csv",sep="\t",append=TRUE)
	cat("\t",file="CR.csv",append=TRUE)
}
cat("\n",file="CR.csv",append=TRUE)
for(name in names(CR)){
	cat(CR[[name]],file="CR.csv",sep="\t",append=TRUE)
	cat("\t",file="CR.csv",append=TRUE)
}


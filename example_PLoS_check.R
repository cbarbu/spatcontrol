library(locfit)
get.sampled<-function(samples=NULL,file="sampled.txt"){
  if(is.null(samples$sampled)){
    # for some reason this is much faster than read.table()
    cat("importing spatial autocorrelation values\n")
    noms<-read.table(file,nrow=1,header=TRUE)
    sampled<-as.data.frame(matrix(scan(skip=1,file=file,sep="\t"),ncol=dim(noms)[2],byrow=TRUE))
    names(sampled)<-names(noms)
  }else{
    sampled<-samples$sampled
  }

  return(sampled)
}
get.cofactors<-function(samples=NULL,file="cofactors.txt"){
  if(is.null(samples$c.vals)){
    cat("import cofactors\n")
    noms<-read.table(file,nrow=1,header=TRUE)
    c.vals<-as.data.frame(matrix(scan(skip=1,file=file,sep="\t"),ncol=dim(noms)[2],byrow=TRUE))
    names(c.vals)<-names(noms)
  }else{
    c.vals<-samples$c.vals
  }
  return(c.vals)
}
get.betas<-function(samples=NULL,file="betasamples.txt"){
  cat("import observers \n")
  noms<-read.table(file,nrow=1,header=FALSE)
  betas<-as.data.frame(matrix(scan(file=file,sep="\t"),ncol=dim(noms)[2],byrow=TRUE))
  return(betas)
}

trace.mcmc<-function(samples=NULL){
  # import autocorrelation parameters and likelihoods

  ## monitored variables/parameters
  sampled<-get.sampled(samples)
  niterations<-dim(sampled)[1]

  par(mfrow=c(2,3))
  plot(1/sampled$T,main="Barrier effect (1/T)",pch=".")
  plot(sampled$f,main="Characteristic distance",pch=".")
  plot(c(1,niterations),c(0,1),type="n",main="Mean SB share",xlab="Index")
  lines(sampled$meanSBshare,pch=".")
  lines(sampled$meanSBshareNoT,main="Mean SB share",pch=".",col="blue",type="p")
  legend("bottomright",c("Actual","If no barrier effect"),pch=20,col=c("black","blue"),pt.cex=0.5)
  plot(sampled$meanBeta,main="Mean observer quality",pch=".")
  plot(sampled$LLH,main="LLH",pch=".")

  ## cofactors
  c.vals<-get.cofactors(samples)
  
  # traces it
  dev.new()
  par(mfrow=c(3,2))
  for(cofactor in names(c.vals)){
    plot(c.vals[,cofactor],main=cofactor,pch=".")
  }

  ## observers
  betas<-get.betas(samples)
  dev.new()
  plot(c(1,dim(betas)[1]),c(0,1),type="n",main="Individual quality of observers",xlab="Iterations",ylab="Detection rate of positive")
  for(insp in 1: dim(betas)[2]){
    lines(betas[,insp],col=sample(colors(),1),pch=".",type="p")
  }

  return(invisible(list(sampled=sampled,c.vals=c.vals,betas=betas)))
}
get.estimate<-function(C,name="",visu=TRUE,leg=TRUE){
  estimate<-c(mean(C),quantile(C,probs=c(0.025,0.5,0.975)))
  names(estimate)[1]<-"mean"
  densfit<-locfit(~C)
  vals<-predict(densfit,estimate)
  if(visu){
    plot(densfit,xlab=name)
    lines(rep(estimate[1],2),c(0,vals[1]),col="black")
    for(q in 2:4){
      lines(rep(estimate[q],2),c(0,vals[q]),col="blue")
    }

    if(leg){
      # legend
      if(mean(densfit$box)>estimate[3]){
	loc<-"topright"
      }else{
	loc<-"topleft"
      }
      legend(loc,c("CrI/med.","Mean"),col=c("blue","black"),lty=1)
    }
  }
  attributes(estimate)$densfit<-densfit
  attributes(estimate)$vals<-vals

  return(estimate)
}
posteriors.mcmc<-function(samples=NULL){
  ### plot posteriors and get estimates

  estimates<-list()
  ## monitored parameters and variables
  sampled<-get.sampled(samples)

  dev.new()
  par(mfrow=c(2,2))
  estimates$T<-get.estimate(1/sampled$T,name="Barrier effect (1/sigma)",leg=FALSE)
  estimates$f<-get.estimate(sampled$f,name="Charac. Dist. ",leg=FALSE)
  estimates$meanBeta<-get.estimate(sampled$meanBeta,name="Mean rate obs.",leg=FALSE)
  estimates$meanSBshare<-get.estimate(sampled$meanSBshare,name="Mean share SG",leg=TRUE)

  ## cofactors
  c.vals<-get.cofactors(samples)
  maxdens<-0
  for(cofactor in names(c.vals)){
    estimates[[cofactor]]<-get.estimate(c.vals[,cofactor],visu=FALSE)
    maxdens<-max(maxdens,attributes(estimates[[cofactor]])$vals)
  }
  dev.new()
  plot(range(c.vals),c(0,maxdens),type="n",xlab="Cofactors' posteriors",ylab="Density")
  for(cofactor in names(c.vals)){
    lines(attributes(estimates[[cofactor]])$densfit,col=which(names(c.vals)==cofactor))
  }
  legend("topleft",names(c.vals),col=(1:length(names(c.vals))),lty=1)

  ## inspectors
  betas<-get.betas(samples)
  for(obs in 1:dim(betas)[2]){
    nameObs<-paste("obs",obs,sep="")
    estimates[[nameObs]]<-get.estimate(betas[,obs],visu=FALSE)
    maxdens<-max(maxdens,attributes(estimates[[nameObs]])$vals)
  }

  dev.new()
  plot(range(betas),c(0,maxdens),type="n",xlab="Observers' quality posteriors",ylab="Density")
  for(obs in 1:dim(betas)[2]){
    nameObs<-paste("obs",obs,sep="")
    lines(attributes(estimates[[nameObs]])$densfit,col=sample(colors(),1))
  }

  return(estimates)
}
graphics.off()
samples<-trace.mcmc()
posteriors.mcmc(samples)


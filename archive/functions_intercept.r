# from a natural mean and a standard deviation 
# get the mean of a corresponding log normal distribution

full.screen<-function(...){
	dev.new(width=13.5,height=9.5,...)
}
printdev<-function(...){
	# if(!exists("X11_plot")){
	if(interactive()){
		out<-try(X11())
		if(class(out)=="try-error"){
			X11_plot<<-FALSE
		}else{
			dev.off();
			X11_plot<<-TRUE
		}
		# }
	}else{
		X11_plot<<-FALSE
	}

	if(X11_plot){
		dev.print(...);
	}
	return(X11_plot);
}
## return the matrix A extended to the given dimensions
## keeping the values in
resized<-function(A,nr=nrow(A),nc=ncol(A)){
	B<-as.matrix(mat.or.vec(nr,nc));
	B[1:(dim(A)[1]),1:dim(A)[2]]<-A
	return(B);
}


meansd2meansdlognorm <-function(mean=mean,sd=NA,sdlog=NA){
	if(is.na(sdlog)){
		sdlog<-sqrt(log((sd/mean)^2+1))
	}
	meanlog<-log(mean)-sdlog^2/2
	return(list(meanlog=meanlog,sdlog=sdlog))
}
meansdlognorm2meansd <-function(meanlog=meanlog,sdlog=NA){
	mean<-exp(meanlog+sdlog^2/2)
	sd<-sqrt(exp(2*meanlog+sdlog^2)*(exp(sdlog^2)-1))
	return(list(mean=mean,sd=sd))
}
# # testing
# fprior<-100
# sdfprior<-50
# flnparam<-meansd2meansdlognorm(mean=fprior,sd=sdfprior)
# 
# # mf<-log(fprior)-sdlf^2/2;
# mf<-flnparam$meanlog
# sdlf<-flnparam$sdlog
# fparam<-meansdlognorm2meansd(meanlog=mf,sdlog=sdlf)
# 
# cat("fprior",fprior,"back",fparam$mean,"sdfprior:",sdfprior,"back:",fparam$sd)


# # should use intersect() in place of this one
# numbers.in.both<-function(a,b){
# 	if(is.null(a)||is.null(b)){
# 		return(NULL)
# 	}
# 	a<-a[order(a)]
# 	b<-b[order(b)]
# 	pa<-1
# 	pb<-1
# 	out<-c()
# 	while(pa <= length(a) && pb <= length(b)){
# 		if(a[pa]==b[pb]){
# 			out<-c(out,a[pa])
# 			pa<-pa+1
# 		}else if(b[pb]>a[pa]){
# 			pa<-pa+1
# 		}else{
# 			pb<-pb+1
# 		}
# 	}
# 	return(out)
# }
# # # testing
# # a<-c(1,2,4,5,9)
# # b<-c(2,3,9,10)
# # numbers.in.both(a,b)

not.in<-function(linesToRemove,originalLength){
	originalLines<-rep(0,originalLength)
	originalLines[linesToRemove]<-1

	linesToKeep<-which(originalLines==0)

	return(linesToKeep)
}
# testing
not.in(c(1,3,4),10)

# Kernel, NB: specific treatment of spam is 2 to 10 times more efficient

sampleku <- function(dim, hyper, Q, x) {
  ku.a <- hyper$k[1];
  ku.b <- hyper$k[2];
  u <- x[1:dim];
  pos.shape <- (0.5*(dim-1) + ku.a);
  pos.scale <- (0.5*as.numeric(u %*% (Q %*% u)) + ku.b^(-1))^(-1);
  ku <- rgamma(n=1, shape=pos.shape, scale=pos.scale);
  return(ku);
}

makeRuint <- function(dim, hyper, Q, ku) {
  t.pr <- hyper$t[2];
  R <- ku*Q;
  diag.spam(R) <- diag.spam(R) + 1;
  R <- cbind.spam(R, rep(1,dim));
  R <- rbind.spam(R, c(rep(1,dim), dim + t.pr));
  return(R);
}
makeRuv <- function(dim,Q,K) {
  Ku <- K[1];
  Kv <- K[2];
  R <- Ku*Q;
  diag.spam(R) <- diag.spam(R) + Kv;
  R <- cbind(R, diag.spam(-1*Kv,dim,dim));
  R <- rbind(R, cbind(diag.spam(-1*Kv,dim,dim), diag.spam((Kv+1), dim, dim)));
  return(R);
}

samplex <- function(dim, hyper, cholR, s, sy0) {
  t.mn <- hyper$t[1];
  t.pr <- hyper$t[2];
  x <- rnorm(n=(dim+1), mean=0, sd=1);
  center <- (s^(-1))*sy0;
  center <- c(center, sum(center) + (t.mn*t.pr));
  center <- backsolve(cholR, forwardsolve(cholR, center));
  x <- backsolve(cholR,x);
  x <- x + center;
  return(x);
}

krig <- function(dim, cholR, x) {
  dir <- c(rep(1, dim), 0);
  dir <- backsolve(cholR, forwardsolve(cholR, dir));
  scale <- (sum(dir[1:dim])^(-1))*sum(x[1:dim]);
  x <- x - scale*dir;
  return(x);
}

samples <- function(dim, cholR, sy0) {
  pos.shape <- (0.5*dim);
  pos.scale <- (0.5*sum(sy0*(sy0 - backsolve(cholR, forwardsolve(cholR, sy0)))));
  invs2 <- rgamma(n=1, shape=pos.shape, scale=pos.scale);
  s <- (sqrt(invs2))^(-1);
  return(s);
}

samplesy0 <- function(dim, s, x, y1) {
  center <- s*(x[1:dim] + x[(dim+1)]);
  lwbd <- rep(0,dim);
  lwbd[(y1==0)] <- -Inf;
  upbd <- rep(0,dim);
  upbd[(y1==1)] <- Inf;
  sy0 <- rtnorm(n=dim, mean=center, sd=s, lower=lwbd, upper=upbd);
  return(sy0);
}

sampley1 <- function(dim, I, z, x, beta) {
  center <- x[1:dim] + x[(dim+1)];
  y1 <- rep(0,dim);
  y1[z$pos] <- 1;
  prob_succes <- c(t(I) %*% beta);
  p <- (1-prob_succes[z$neg])*pnorm(q=center[z$neg], mean=0, sd=1, lower.tail=TRUE); # prob infested, not detected
  q <- pnorm(q=center[z$neg], mean=0, sd=1, lower.tail=FALSE); # prob not infested
  pneg <- p/(p + q); # prob obsevation negative
  pna <- pnorm(q=center[z$na], mean=0, sd=1);
  y1[z$neg] <- rbinom(n=length(z$neg), size=1, prob=pneg);
  y1[z$na] <- rbinom(n=length(z$na), size=1, prob=pna);
  y1 <- as.numeric(y1);
  return(y1);
}
sampley <- function(dim,u,yprime) {
  center <- u;
  lwbd <- rep(0,dim);
  lwbd[(yprime==0)] <- -Inf;
  upbd <- rep(0,dim);
  upbd[(yprime==1)] <- Inf;
  y <- rtnorm(n=dim, mean=center, sd=1, lower=lwbd, upper=upbd);
  return(y);
}

sampleyprime <- function(dim,w,betaprime,zpos,zneg,zNA) {
  yprime <- rep(0,dim);
  yprime[zpos] <- 1;
  p <- (1-betaprime[zneg])*(1-pnorm(q=0, mean=(w[zneg]), sd=1));
  q <- (pnorm(q=0, mean=(w[zneg]), sd=1));
  area<-p+q
  area[area==0]<-1; # avoid NA and associated warnings
  Pzneg <- p/(area);
  PzNA <- 1-pnorm(q=0, mean=(w[zNA]), sd=1);
  yprime[zneg] <- rbinom(length(zneg),1,Pzneg);
  yprime[zNA] <-  rbinom(length(zNA),1,PzNA); # w[zNA]  
  return(yprime);
}

sampleb <- function(dim, hyper, I, z, y1) {
  beta.a <- hyper$b[1];
  beta.b <- hyper$b[2];
  pos <- rep(0, dim);
  pos[z$pos] <- 1;
  total <- c(I %*% y1);
  ident <- c(I %*% (y1*pos));
  pos.a <- ident + beta.a; 
  pos.b <- total - ident + beta.b; 
  beta <- rbeta(n=nrow(I), shape1=pos.a, shape2=pos.b); 
  return(beta);
}
lik.f<-function(f,mf,sdlf,log=TRUE){
	LLH<-dlnorm(f,mf,sdlf,log=log)
	# LLH<-0*f # flat prior
	return(LLH);
}
acceptf<-{};
sample_f <- function(u,Ku,T,logsdfprop,f,mf,sdlf,Q,LLHu,AS,SB,cholQ=NULL,Dmat=NULL){
    # # sampling symmetrically arround f (may be buggy)	
	# flogmean<-meansd2meansdlognorm(mean=f,sdlog=logsdfprop)$meanlog
    # f_prop<-rlnorm(1,flogmean,logsdfprop);
    # flogmeanprop<-meansd2meansdlognorm(f_prop,sdlog=logsdfprop)$meanlog
    # hasting_term=dlnorm(f,flogmeanprop,logsdfprop,log=TRUE)-dlnorm(f_prop,flogmean,logsdfprop,log=TRUE)
    
    # sampling not symmetrically (but simple and correct) arround T
    f_prop<-rlnorm(1,log(f),logsdfprop);
    hasting_term=dlnorm(f,log(f_prop),logsdfprop,log=TRUE)-dlnorm(f_prop,log(f),logsdfprop,log=TRUE);
    
	Qprop<-QfromfT(Dmat,AS,SB,f=f_prop,T=T);
	# Qprop<-Qfromf(Dmat,f=f_prop);

	cholQprop<-get.cholMat(Qprop,cholQ)
	LLHuprop<-fast.llh.ugivQ(dimension,u,Qprop,Ku,cholQ=cholQprop);

	llhfprop<-lik.f(f_prop,mf,sdlf);
	llhf<-lik.f(f,mf,sdlf);

	LLHproposal <- LLHuprop+llhfprop; 
	LLH <-LLHu+llhf;
	lnr <- LLHproposal-LLH+hasting_term;

	cat("f:",f," LLH:",LLH,"(",LLHu,"+",llhf,") f prop:",f_prop," LLHproposal:",LLHproposal,"(",LLHuprop,"+",llhfprop,") h_t:",hasting_term," lnr",lnr,sep="");

	if(lnr>=log(runif(1))) {
		f <- f_prop;
		Q<-Qprop;
		cholQ<-cholQprop;
		LLHu<-LLHuprop;
		acceptf<<-c(acceptf,1);
		cat(" accept 1\n");
	}else{
		acceptf<<-c(acceptf,0);
		cat("accept 0\n");
	}
	return(list(f=f,Q=Q,LLHu=LLHu,cholQ=cholQ));
}
lik.T<-function(f,mT,sdlT,log=TRUE){
	LLH<-dlnorm(f,mT,sdlT,log=log)
	# LLH<-0*f # flat prior
	return(LLH);
}

acceptT<-{};
sample_T <- function(u,Ku,f,T,logsdTprop,mT,sdT,Q,LLHu,AS,SB,cholQ=NULL,Dmat=NULL){
    # # sampling symmetrically arround T (may be buggy)	
    # Tlogmean<-meansd2meansdlognorm(mean=T,sdlog=logsdTprop)$meanlog
    # T_prop<-rlnorm(1,Tlogmean,logsdTprop);
    # Tlogmeanprop<-meansd2meansdlognorm(mean=T_prop,sdlog=logsdTprop)$meanlog
	# hasting_term=dlnorm(T,Tlogmeanprop,logsdTprop,log=TRUE)-dlnorm(T_prop,Tlogmean,logsdTprop,log=TRUE);

    # sampling not symmetrically (but simple and correct) arround T	
    T_prop<-rlnorm(1,log(T),logsdTprop);
    hasting_term=dlnorm(T,log(T_prop),logsdTprop,log=TRUE)-dlnorm(T_prop,log(T),logsdTprop,log=TRUE);

    # calculate the LLH (common for both samplings)
	Qprop<-QfromfT(Dmat,AS,SB,f=f,T=T_prop);

	cholQprop<-get.cholMat(Qprop,cholQ)
	LLHuprop<-fast.llh.ugivQ(dimension,u,Qprop,Ku,cholQ=cholQprop);
	llhTprop<-lik.T(T_prop,mT,sdlT);
	llhT<-lik.T(T,mT,sdlT);
	LLHproposal <- LLHuprop+llhTprop; 
	LLH <- 	       LLHu+llhT;
	lnr <- LLHproposal-LLH+hasting_term;

	cat("T:",T," LLH:",LLH,"(",LLHu,"+",llhT,") T prop:",T_prop," LLHproposal:",LLHproposal,"(",LLHuprop,"+",llhTprop,") h_t:",hasting_term," lnr",lnr,sep="");

	if(lnr>=log(runif(1))) {
		T <- T_prop;
		Q<-Qprop;
		cholQ<-cholQprop;
		LLHu<-LLHuprop;
		acceptT<<-c(acceptT,1);
		cat(" accept 1\n");
	}else{
		acceptT<<-c(acceptT,0);
		cat("accept 0\n");
	}
	return(list(T=T,Q=Q,LLHu=LLHu,cholQ=cholQ));
}
acceptfT={}
sample_fT <- function(u,K,T,logsdTprop,mT,sdlT,f,logsdfprop,sdCoupleFact,mf,sdlf,Q,LLHu,AS,SB,cholQ=NULL,Dmat=NULL){
	f_prop<-rlnorm(1,log(f),logsdfprop*sdCoupleFact);
	T_prop<-rlnorm(1,log(T),logsdTprop*sdCoupleFact);

	Qprop<-QfromfT(Dmat,AS,SB,f=f_prop,T=T_prop);

	h1=dlnorm(f,log(f_prop),logsdfprop,log=TRUE)-dlnorm(f_prop,log(f),logsdfprop,log=TRUE);
	h2=dlnorm(T,log(T_prop),logsdTprop,log=TRUE)-dlnorm(T_prop,log(T),logsdTprop,log=TRUE);
	h<-h1+h2;

	cholQprop<-get.cholMat(Qprop,cholQ)
	LLHuprop<-fast.llh.ugivQ(dimension,u,Qprop,Ku,cholQ=cholQprop);

	llhfprop<-lik.f(f_prop,mf,sdlf);
	llhTprop<-lik.T(T_prop,mT,sdlT);
	llhf<-lik.f(f,mf,sdlf);
	llhT<-lik.T(T,mT,sdlT);

	LLHproposal <- LLHuprop+llhfprop+llhTprop; 
	LLH <- 	       LLHu+llhf+llhT;
	lnr <- LLHproposal-LLH+h;

	cat("f :",f," T :",T," LLH:",LLH,"(",LLHu,"+",llhf,"+",llhT,")\nfp:",f_prop," Tp:",T_prop," LLHproposal:",LLHproposal,"(",LLHuprop,"+",llhfprop,"+",llhTprop,") h_t:",h," lnr",lnr,sep="");

	if(lnr>=log(runif(1))) {
		f <- f_prop;
		T <- T_prop;
		Q<-Qprop;
		cholQ<-cholQprop;
		LLHu<-LLHuprop;
		acceptfT<<-c(acceptfT,1);
		cat(" accept 1\n");
	}else{
		acceptfT<<-c(acceptfT,0);
		cat("accept 0\n");
	}
	return(list(f,T,Q,LLHu,cholQ));
}


get.cholMat<-function(Mat,cholMat=NULL){ # compute cholMat, update if possible 
	if(is.null(cholMat)){
		cholMat<-chol.spam(Mat);
	}else{
		out<-try(cholMat<-update.spam.chol.NgPeyton(cholMat,Mat));
		if(class(out)=="try-error"){ 
			cholMatbad<<-cholMat;
			Matbad<<-Mat;
			dump(list("cholMatbad","Matbad"),file="badMat.r")
			cholMat<-chol.spam(Mat);
		}
	}
	return(cholMat)
}
# general function for likelihood of u given Q
llh.ugivQ<-function(dimension,u,Q,Ku,cholQ=NULL){
	# cholQ<-chol(Q);
	cholQ<-get.cholMat(Q,cholMat=cholQ)

	exp_part<- Ku*(t(u-mean(u))%*%Q%*%(u-mean(u))); # version with mean of u has no importance
	# exp_part<- Ku*(t(u)%*%Q%*%(u)); # prior mean is at 0
	det_part<- dimension*log(Ku)+2*determinant(cholQ)$modulus
	logpi<-dimension*log(2*pi);
	LLH= -1/2*(exp_part-det_part+logpi); 
	# cat("LLH",LLH,"exppart:",exp_part,"det_part",det_part,"logpi",logpi,"\n");
	return(LLH);
}
fast.llh.ugivQ<-function(dimension,u,Q,Ku,cholQ,uQu){ # when cholQ is for sure known
	exp_part<- Ku*(t(u-mean(u))%*%Q%*%(u-mean(u))); # version with mean of u has no importance
	# exp_part<- Ku*(t(u)%*%Q%*%(u)); # prior mean is at 0
	det_part<- dimension*log(Ku)+2*determinant(cholQ)$modulus
	logpi<-dimension*log(2*pi);
	LLH= -1/2*(exp_part-det_part+logpi); 
	return(LLH);
}
llh.vgivKv<-function(dimension,v,Kv){
	LLH<-sum(dnorm(v,mean=0,sd=sqrt(1/Kv),log=TRUE))

	return(LLH);
}
llh.ygivw<-function(y,w){
	## return the loglikelihood of y given w
	logpi<-length(y)*log(2*pi);
	LLH<- -1/2 * t(y-w)%*%(y-w) - logpi;
	return(LLH);
}
llh.zgivy<-function(y,zpos,zneg,bivect){
	probypos<-(y>0)
	LLHpos<-sum(log(probypos[zpos]*bivect[zpos]));
	LLHneg<-sum(log(probypos[zneg]*(1-bivect[zneg])+(1-probypos[zneg])));
	LLH<-LLHpos+LLHneg;
	return(LLH);
}
llh.zgivw<-function(w,zpos,zneg,bivect){
	probypos<-1-pnorm(0,w,1)
	LLHpos<-sum(log(probypos[zpos]*bivect[zpos]));
	LLHneg<-sum(log(probypos[zneg]*(1-bivect[zneg])+(1-probypos[zneg])));
	LLH<-LLHpos+LLHneg;
	return(LLH);
}

gibbs <- function(dim, hyper, I, z, s, state,MHsd,AutMats) {
  ku <- state$ku;
  x <- state$x;
  u <- x[1:dimension]
  sy0 <- state$sy0;
  y1 <- state$y1;
  beta <- state$beta;
  T<- state$T
  f<- state$f
  Q<-state$Q
  cholQ<-state$cholQ
  LLHu<-state$LLHu

  # update Q
  if(use.f){
	  out<-sample_f(u,ku,T,MHsd$f,f,hyper$f[1],hyper$f[2],Q,LLHu=LLHu,AutMats$AS,AutMats$SB,cholQ=cholQ,Dmat=AutMats$Dmat); 
	  f<-out$f
	  Q<-out$Q
	  cholQ<-out$cholQ
	  LLHu<-out$LLHu
  }
  if(use.streets){
	  out<-sample_T(u,ku,f,T,MHsd$T,hyper$T[1],hyper$T[2],Q,LLHu,AutMats$AS,AutMats$SB,cholQ=cholQ,Dmat=AutMats$Dmat);
	  T<-out$T;
	  Q<-out$Q;
	  cholQ<-out$cholQ;
	  LLHu<-out$LLHu;
  }

  R <- makeRuint(dim=dim, hyper=hyper, Q=Q, ku=ku);
  # cholR <- chol.spam(R, memory=list(nnzcolindices=4e6),);
  cholR <<- update.spam.chol.NgPeyton(cholR,R, memory=list(nnzcolindices=4e6));
  rm(R);
  x <- samplex(dim=dim, hyper=hyper, cholR=cholR, s=s, sy0=sy0);
  x <- krig(dim=dim, cholR=cholR, x=x);
  ku <- sampleku(dim=dim, hyper=hyper, Q=Q, x=x);
  sy0 <- samplesy0(dim=dim, s=s, x=x, y1=y1);
  y1 <- sampley1(dim=dim, I=I, z=z, x=x, beta=beta);
  beta <- sampleb(dim=dim, hyper=hyper, I=I, z=z, y1=y1);

  LLHu<-fast.llh.ugivQ(dimension,x[1:dimension],Q,ku,cholQ); # important for f/T sampling

  state <- list(ku=ku, x=x, sy0=sy0, y1=y1, beta=beta,T=T,f=f,Q=Q,cholQ=cholQ,LLHu=LLHu);
  return(state);
}

lik <- function(dim, hyper, Q, I, state) {
  t.mn <- hyper$t[1];
  t.pr <- hyper$t[2];
  ku.a <- hyper$k[1];
  ku.b <- hyper$k[2];
  beta.a <- hyper$b[1];
  beta.b <- hyper$b[2];
  ku <- state$ku;
  x <- state$x;
  y1 <- state$y1;
  beta <- state$beta;
  pos <- rep(0, dim);
  pos[z$pos] <- 1;
  total <- c(I %*% y1);
  ident <- c(I %*% (y1*pos));
  t <- x[(dim+1)];
  u <- x[1:dim];
  L1 <- (-0.5*t.pr)*((t - t.mn)^2);
  L2 <- (ku.a - 1)*log(ku) - (ku/ku.b);
  L3 <- (0.5*(dim-1))*log(ku) - (0.5*ku)*as.numeric(u %*% (Q %*% u)); # missing the determinant part
  L4 <- sum(y1*pnorm(q=(u+t), mean=0, sd=1, log.p=TRUE)) + sum((1 - y1)*pnorm(q=(u+t), mean=0, sd=1, lower.tail=FALSE, log.p=TRUE));
  L5 <- sum((ident + beta.a - 1)*log(beta)) + sum((total - ident + beta.b - 1)*log(1 - beta));
  L <- L1 + L2 + L3 + L4 + L5;
  return(L);
}
generator <- function(dim, hyper, Q, I, z, s, ku, t) {
  Q.temp <- (ku)*Q;
  diag.spam(Q.temp) <- diag.spam(Q.temp) + (ku)*(1e-6);
  cholQ.temp <- chol(Q.temp);
  u <- rnorm(n=dim, mean=0, sd=1);
  u <- backsolve(cholQ.temp,u);
  u <- u - mean(u);
  y0 <- rnorm(n=dimension, mean=(u + t), sd=1);
  sy0 <- s*y0;
  y1 <- as.numeric(y0 > 0);
  beta <- rbeta(n=nrow(I), shape1=hyper$b[1], shape2=hyper$b[2]);
  p <- c(t(I) %*% beta)*(y1);
  zsim <- rbinom(n=dim, size=1, prob=p);
  na <- z$na;
  zsim[na] <- 9;
  pos <- which(zsim==1);
  neg <- which(zsim==0);
  zsim <- list(pos=pos, neg=neg, na=na);
  rm(Q.temp,cholQ.temp);
  state <- list(ku=ku, x=c(u,t), sy0=sy0, y1=y1, beta=beta);
  return(list(state=state, zsim=zsim));
}

moran <- function(dim, W, y1) {
  num <- dim*as.numeric((y1 - mean(y1)) %*% (W %*% (y1 - mean(y1))));
  den <- sum(W)*as.numeric((y1 - mean(y1)) %*% (y1 - mean(y1)));
  return((num/den));
}

neighbor <- function(dim, easting, northing, i) {
  z <- rep(0, dimension);
  z[which(W[i,] != 0)] <- 5;
  z[i] <- 10;
  res <- data.frame(z);
  coordinates(object=res) <- cbind(easting, northing);
  spplot(res);
}
updateMHparam<-function(MHsd){
	adaptOK<-TRUE
	if(use.f){
		logsdfprop<-MHsd$f
		rateaccept<-mean(tail(acceptf,20))
		cat("accept rate f:",rateaccept);
		if(rateaccept<lowAcceptRate){
			logsdfprop<-0.9*logsdfprop;
			cat("update of logsdfprop to:",logsdfprop);
			adaptOK<-FALSE
		}else if(rateaccept>highAcceptRate){
			logsdfprop<-1.1*logsdfprop;
			cat("update of logsdfprop to:",logsdfprop);
			adaptOK<-FALSE
		}
		MHsd$f<-logsdfprop
	}
	if(use.streets){
		rateaccept<-mean(tail(acceptT,20))
		cat("accept rate T:",rateaccept);
		if(rateaccept<lowAcceptRate){
			logsdTprop<-0.9*logsdTprop;
			cat("update of logsdTprop to:",logsdTprop);
			adaptOK<-FALSE
		}else if(rateaccept>highAcceptRate){
			logsdTprop<-1.1*logsdTprop;
			cat("update of logsdTprop to:",logsdTprop);
			adaptOK<-FALSE
		}
	}

	MHsd$OK<-adaptOK

	return(MHsd)
}
sample_u <- function(dimension,Q,K,y,cholQ=NULL){
	R <- K[1]*Q+diag.spam(1,dimension);
	center <- y;
	u <- rmvnorm.canonical(n=1, b=center, Q=R,Rstruct=cholQ);
	
	return(drop(u));
}
source("sample_y_direct.r")

sampleK <- function(dim,Q,x,K.hyper) {
  Ku.a <- K.hyper[1];
  Ku.b <- K.hyper[2];
  Kv.a <- K.hyper[3];
  Kv.b <- K.hyper[4];
  u <- x[1:dim];
  v <- x[(dim + (1:dim))] - u;
  u.pshape <- (0.5*(dim-1) + Ku.a);
  # u.pscale <- (0.5*as.numeric(u %*% (Q %*% u)) + Ku.b^(-1))^(-1); # the bad thing you can do against convergence
  u.pscale <- (0.5*as.numeric(u-mean(u)) %*% (Q %*% (u-mean(u))) + Ku.b^(-1))^(-1); # centered u
  v.pshape <- (0.5*dim + Kv.a);
  v.pscale <- (0.5*as.numeric(v %*% v) + Kv.b^(-1))^(-1);
  Ku <- rgamma(n=1, shape=u.pshape, scale=u.pscale);
  Kv <- rgamma(n=1, shape=v.pshape, scale=v.pscale);
  K <- c(Ku,Kv);
  return(K);
}
acceptKv={}
sampleKvMHunifSigma<-function(Kv,v,logsdDraw){
	sig<-sqrt(1/Kv)

	sigProp<-rlnorm(1,mean=log(sig),sd=logsdDraw)
	
	LLHproposal<-sum(dnorm(v,mean=0,sd=sigProp,log=TRUE))
	LLH<-sum(dnorm(v,mean=0,sd=sig,log=TRUE))

	hasting_term=dlnorm(sig,log(sigProp),logsdDraw,log=TRUE)-dlnorm(sigProp,log(sig),logsdDraw,log=TRUE);

	lnr <- LLHproposal-LLH+hasting_term;

	cat("sigv:",sig," LLH:",LLH," sigv prop:",sigProp," LLHproposal:",LLHproposal," lnr",lnr,sep="");

	if(lnr>=log(runif(1))) {
		Kv <- 1/(sigProp^2);
		acceptKv<<-c(acceptKv,1);
		cat(" accept 1\n");
	}else{
		acceptKv<<-c(acceptKv,0);
		cat("accept 0\n");
	}
	return(list(Kv=Kv));
}
# # test sampleKvMHunifSigma
# ntest<-10000
# Kv<-0.1*Kv.r
# sampledKv<-rep(0,ntest)
# for(i in 1: ntest){
# Kv<-sampleKvMHunifSigma(Kv,v.r,0.05)$Kv
# sampledKv[i]<-Kv
# }
# par(mfrow=c(1,2))
# plot(sampledKv)
# abline(h=Kv.r,col=2)
# abline(h=1/sd(v.r)^2,col=4)
# hist(sampledKv[-(1:100)])
# abline(v=Kv.r,col=2)
# abline(v=1/sd(v.r)^2,col=4)
# mean(acceptKv)
# mean(sampledKv[-(1:100)])

acceptKu={}
sampleKuMHunifSigma<-function(Ku,u,logsdDraw,Q,LLHu=NULL,cholQ=NULL,ucQuc=NULL){
	sig<-sqrt(1/Ku)

	sigProp<-rlnorm(1,mean=log(sig),sd=logsdDraw)
	
	LLHuprop<-fast.llh.ugivQ(dim(Q)[1],u,Q,1/sigProp^2,cholQ)
	LLHproposal<-LLHuprop
	if(is.null(LLHu)){
		LLHu<-llh.ugivQ(dim(Q)[1],u,Q,Ku,cholQ=cholQ)
	}
	LLH<-LLHu

	hasting_term=dlnorm(sig,log(sigProp),logsdDraw,log=TRUE)-dlnorm(sigProp,log(sig),logsdDraw,log=TRUE);

	lnr <- LLHproposal-LLH+hasting_term;

	cat("sigu:",sig," LLH:",LLH," sigu prop:",sigProp," LLHproposal:",LLHproposal," lnr",lnr,sep="");

	if(lnr>=log(runif(1))) {
		Ku <- 1/(sigProp^2);
		LLHu<-LLHuprop;
		acceptKu<<-c(acceptKu,1);
		cat(" accept 1\n");
	}else{
		acceptKu<<-c(acceptKu,0);
		cat("accept 0\n");
	}
	return(list(Ku=Ku,LLHu=LLHu));
}
# # test sampleKuMHunifSigma
# ntest<-1000
# Ku<-0.1*Ku.r
# sampledKu<-rep(0,ntest)
# LLHu<-NULL
# for(i in 1: ntest){
# 	cat(i)
# 	out<-sampleKuMHunifSigma(Ku,u.r,0.05,Q.r,LLHu=LLHu,cholQ=cholQ.r)
# 	Ku<-out$Ku
# 	LLHu<-out$LLHu
# 	sampledKu[i]<-Ku
# }
# par(mfrow=c(1,2))
# plot(sampledKu)
# abline(h=Ku.r,col=2)
# abline(h=1/sd(u.r)^2,col=4)
# hist(sampledKu[-(1:100)])
# abline(v=Ku.r,col=2)
# abline(v=1/sd(u.r)^2,col=4)
# mean(acceptKu)
# mean(sampledKu[-(1:100)])

## ok for both Kv and Kc (and could be used for Ku with dim<-dimension-1)
sampleKg <- function(dim,Q,field,K.shape,K.scale,n=1) {
  pshape <- (0.5*dim + K.shape);
  pscale <- (0.5*as.numeric(t(field) %*%Q%*%field) + K.scale^(-1))^(-1);
  K <- rgamma(n=n, shape=pshape, scale=pscale);
  return(K);
}
sampleKu <- function(dim,Q,u,Kushape,Kuscale) {
	pos.shape <- (0.5*(dim-1) + Kushape);
	pos.scale <- (0.5*as.numeric(u %*% Q %*% u) + Kuscale^(-1))^(-1);
	Ku <- rgamma(n=1, shape=pos.shape, scale=pos.scale);
	return(Ku);
}

# makeRc<- function(Q,K,Qc){
# 	# make Rc using kronecker product
# 	# substantively faster
# 	Ku <- K[1];
# 	Kv <- K[2];
# 	Kc <- K[3];
# 	KQc<-Kc*Qc;
# 	R <- Ku*Q;
# 	dimension<-nrow(Qc)
# 	zero<-spam(0,dimension,dimension)
# 	R <- cbind.spam(R, zero);
# 	R <- rbind.spam(R, cbind.spam(zero, diag.spam(1,nrow(Qc))));
# 
# 	K<-kronecker.spam(matrix(c(1,-1,-1,1),2),KQc)
# 	R <-R+K;
# 	return(R);
# }
makeRc<- function(Q,K,Qc){
	Ku <- K[1];
	Kv <- K[2];
	Kc <- K[3];
	R <- Ku*Q;
	dimension<-nrow(Qc)
	zero<-spam(0,dimension,dimension)
	R <- cbind.spam(R, zero);
	R <- rbind.spam(R, cbind.spam(zero, diag.spam(1,nrow(Qc))));

	R <-R+Kc*kQc;
	return(R);
}
baseRc<- function(Q,K,Qc){
	Ku <- K[1];
	Kv <- K[2];
	Kc <- K[3];
	R <- Ku*Q;
	dimension<-nrow(Qc)
	zero<-spam(0,dimension,dimension)
	R <- cbind.spam(R, zero);
	R <- rbind.spam(R, cbind.spam(zero, zero));

	return(R);
}
basediag1<-function(dimension){
	zero<-spam(0,dimension,dimension)
	Rdiag1 <- cbind.spam(zero, zero);
	Rdiag1 <- rbind.spam(Rdiag1, cbind.spam(zero, diag.spam(1,nrow(Qc))));
	return(Rdiag1);
}

updateRc<- function(Q,RQc,QentriesInR,diag1entriesInR,K){
	# make Rc using kronecker product
	# substantially faster
	Ku <- K[1];
	Kc <- K[3];
	RKQc <- Kc*RQc;

	RKQc@entries[QentriesInR]<-Ku*Q@entries+RKQc@entries[QentriesInR];
	RKQc@entries[diag1entriesInR]<-RKQc@entries[diag1entriesInR]+1;
	return(RKQc);
}


samplexuvcof<-function(dimension,Q,K,y,Qc,cholR){
	x <- rnorm(n=(2*dimension), mean=0, sd=1);
	center <- c(rep(0,dimension), y);
	# R <- makeRc(Q,K,Qc);
	R <-updateRc(Q,RQc,QentriesInR,diag1entriesInR,K) ## consistently faster than makeRc
	if(is.null(cholR)){
		cholR <- chol.spam(R, memory=list(nnzcolindices=4e6),);
	}else{
		cholR <- update.spam.chol.NgPeyton(cholR,R);
	}

	center <- backsolve(cholR, forwardsolve(cholR, center));
	x <- backsolve(cholR,x);
	x <- x + center;
	cholR<<-cholR;
	return(x);
}
##  as we quickly run into memory troubles with the gibbs sampler
##  a metropolis hasting version is welcome
# pi(c | y,u,v)  \propto pi(y | c.val,u,v)pi(c.val)
#                \propto N(c.comp+u+v,1) N(c.val,Kc)
llh.c<-	function(c.val,c.comp,Kc,wnoc,y){
	w<-wnoc+c.comp;
	LLH<-llh.ygivw(y,w)+sum(dnorm(c.val,mean=0,sd=sqrt(1/Kc),log=TRUE));
	return(LLH);
}
acceptc.val<-{}
mhsamplec<-function(c.val,c.comp,c.map,sdc.val,Kc,wnoc,y){
	nbc<-length(c.val);
	c.valprop<-rnorm(nbc,c.val,sdc.val);
	# hasting_term=dnorm(c.val,log(c.valprop),logsdc.val,log=TRUE)-dnorm(c.valprop,log(c.val),logsdc.val,log=TRUE);
	c.compprop<-drop(c.map%*%c.valprop);
	LLHproposal <- llh.c(c.valprop,c.compprop,Kc,wnoc,y); 
	LLH <-llh.c(c.val,c.comp,Kc,wnoc,y);
	lnr <- LLHproposal-LLH# +hasting_term;
	cat("c.val:",c.val,"LLH:",LLH,"c.valprop:",c.valprop,"LLH proposal:",LLHproposal,"lnr",lnr);

	if(lnr>=log(runif(1))) {
		c.val <- c.valprop;
		c.comp<-c.compprop;
		acceptc.val<<-c(acceptc.val,1);
		cat(" accept 1\n");
	}else{
		acceptc.val<<-c(acceptc.val,0);
		cat("accept 0\n");
	}

	return(list(c.val,c.comp));
}
samplexuv <- function(dim,Q,K,y,cholR=NULL) {
  x <- rnorm(n=(2*dim), mean=0, sd=1);
  center <- c(rep(0,dim), y);
  R <- makeRuv(dim,Q,K);
  if(is.null(cholR)){
	  cholR <- chol.spam(R, memory=list(nnzcolindices=4e6),);
  }else{
	  cholR <- update.spam.chol.NgPeyton(cholR,R);
  }

  center <- backsolve(cholR, forwardsolve(cholR, center));
  x <- backsolve(cholR,x);
  x <- x + center;
  cholR<<-cholR;
  return(x);
}
acceptx.val<-{}
samplexMH<-function(dimension,x,Q,K,y,zpos,zneg,bivect,c.comp,sdx.val,cholQ=NULL) {
	
	u<-x[1:dimension]
	v<-x[dimension+(1:dimension)]-u;
	w<-x[(1:dimension)+dimension]+c.comp

	xprop<-rnorm(length(x),mean=x,sd=sdx.val)
	# xprop <- samplexuv(dimension,Q,K,y-c.comp,cholR);

	uprop<-xprop[(1:dimension)]
	wprop<-xprop[(1:dimension)+dimension]+c.comp
	vprop<-xprop[(1:dimension)+dimension]-uprop;

	llhzprop<-llh.zgivw(wprop,zpos,zneg,bivect);
	llhuprop<-llh.ugivQ(dimension,uprop,Q,K[1],cholQ=cholQ);
	llhvprop<-llh.vgivKv(dimension,vprop,K[2]);
	LLHproposal<-llhzprop+llhuprop+llhvprop;

	llhz<-llh.zgivw(w,zpos,zneg,bivect)
	llhu<-llh.ugivQ(dimension,u,Q,K[1],cholQ=cholQ)
	llhv<-llh.vgivKv(dimension,v,K[2])
	LLH<-llhz+llhu+llhv;
	cat("orig llhs: u(",llhu,") v(",llhv,") z(",llhz,")\n")
	cat("prop llhs: u(",llhuprop,") v(",llhvprop,") z(",llhzprop,")\n")

	lnr <- LLHproposal-LLH# +hasting_term;
	cat("u mean:",mean(u),"LLH:",LLH,"uprop mean",mean(uprop),"LLH prop:",LLHproposal,"lnr",lnr);

	if(lnr>=log(runif(1))) {
		x <- xprop;
		acceptx.val<<-c(acceptx.val,1);
		cat(" accept 1\n");
	}else{
		acceptx.val<<-c(acceptx.val,0);
		cat("accept 0\n");
	}

	return(x);
}
fastsamplexuv <- function(dim,cholR,y) {
  x <- rnorm(n=(2*dim), mean=0, sd=1);
  center <- c(rep(0,dim), y);
  center <- backsolve(cholR, forwardsolve(cholR, center));
  x <- backsolve(cholR,x);
  x <- x + center;
  return(x);
}

samplebeta <- function(zpos,zneg,matrix,yprime,a,b) {
	yp.positive <- yprime;
	yp.negative <- yprime; 
	yp.positive[-zpos] <- 0; # z- turn y+ into 0 a.pos being then the sum of the y+ well observed
	yp.negative[-zneg] <- 0; # z+ turn y+ into 0 a.neg being then the sum of the y+ badly observed 
	# the NA are not apparent here
	a.pos <- (as.vector(t(matrix) %*% yp.positive) + a); 
	b.pos <- (as.vector(t(matrix) %*% yp.negative) + b); 
	beta <- rbeta(n=ncol(matrix), shape1=a.pos, shape2=b.pos); 
	return(beta);
}
generate_z <-function(y,b,zNA){
	# y the probit continuous value
	# b the vector with per house the discovery rate of inspectors when y>0
	# if NA, z is set to NA
	b[zNA]=0;
	pre_z<-(y>0);
	zypos<-rbinom(sum(pre_z),1,b[pre_z]);
	pre_z[pre_z]<-zypos;
	pre_z[zNA]<-NA
	return(pre_z);
}
# allow to plot boxplots with custom values
boxplot.free<-function(x,breaks=c(0.025,0.25,0.5,0.75,0.975),rm.out=TRUE,...){
	out<-boxplot(x,plot=FALSE)
	out$stats[1,]<-apply(x,2,quantile, probs = breaks[1], na.rm = TRUE, names= FALSE)
	out$stats[2,]<-apply(x,2,quantile, probs = breaks[2], na.rm = TRUE, names= FALSE)
	out$stats[3,]<-apply(x,2,quantile, probs = breaks[3], na.rm = TRUE, names= FALSE)
	out$stats[4,]<-apply(x,2,quantile, probs = breaks[4], na.rm = TRUE, names= FALSE)
	out$stats[5,]<-apply(x,2,quantile, probs = breaks[5], na.rm = TRUE, names= FALSE)
	# cat("outstats:",out$stats,"\n")
	if(rm.out){
		out$out<-c()
		out$group<-c()
	}
	bxp(out,...)
	return(out)
}
adjust.lim<-function(limsmall,limbig,steps){
	# limsmall: c(min,max) for coord with smallest range
	# limbig: c(min,max) for coord with biggest range
	steps.small<-steps
	limsmall<-c(min(limsmall),max(limsmall))
	stepsize<-(limsmall[2]-limsmall[1])/steps.small

	limbig<-c(min(limbig),max(limbig))
	size.big.init<-(limbig[2]-limbig[1])
	steps.big<-ceiling(size.big.init/stepsize)
	shift.big<-(steps.big*stepsize-size.big.init)/2
	limbig<-c(limbig[1]-shift.big,limbig[2]+shift.big)

	return(list(limsmall=limsmall,limbig=limbig,stepsize=stepsize))
}
#test
# adjust.lim(c(0,4),c(10.1,0),10)

library(sp)
library(spam)
source("spam_complement.r")
make.col.persp<-function(z,nbcol=100,color.function=jet.colors){
	# make col vector for persp from base package to get colors according to z
	# nrz: number of columns in z
	nrz <- nrow(z)
	ncz <- ncol(z)

	nbcol <- 100
	color <- color.function(nbcol)
	# Compute the z-value at the facet centres
	zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
	# Recode facet z-values into color indices
	facetcol <- cut(zfacet, nbcol)
	zcol<-color[facetcol]
	return(zcol)
}
expKernel<-function(T,matdist,f){
	if(class(matdist)=="spam"){
		K<-matdist;
		K@entries<- T*exp(-matdist@entries/f);
	}else{
		K<- T*exp(-matdist/f);
	}
	# K<- T*exp(-log(2)*matdist/f); ## normalized to have 0.5 when dist=f
	return(K);
}
gaussianKernel<-function(T,matdist,f){
	if(class(matdist)=="spam"){
		K<-matdist;
		K@entries<- T*exp(-(matdist@entries*matdist@entries)/f^2);
	}else{
		K<- T*exp(-(matdist*matdist)/f^2);
	}
	# K<- T*exp(-log(2)*(matdist*matdist)/f^2);## normalized to have 0.5 when dist=f
	return(K);
}
# system.time(for(i in 1:10){A<-gaussianKernel(1,Dmat,20)})
cauchyKernel<-function(T,matdist,f){
	if(class(matdist)=="spam"){
		K<-matdist;
		K@entries<- T/(1+(matdist@entries*matdist@entries)/f^2);
	}else{
		K<- T/(1+(matdist*matdist)/f^2);
	}
	return(K);
}
geometricKernel<-function(T,matdist,f){
	if(class(matdist)=="spam"){
		K<-matdist;
		K@entries<- T/(1+matdist@entries/f);
	}else{
		K<- T/(1+matdist/f);
	}
	return(K);
}
adjust.lim<-function(limsmall,limbig,steps){
	# limsmall: c(min,max) for coord with smallest range
	# limbig: c(min,max) for coord with biggest range
	steps.small<-steps
	limsmall<-c(min(limsmall),max(limsmall))
	stepsize<-(limsmall[2]-limsmall[1])/steps.small

	limbig<-c(min(limbig),max(limbig))
	size.big.init<-(limbig[2]-limbig[1])
	steps.big<-ceiling(size.big.init/stepsize)
	shift.big<-(steps.big*stepsize-size.big.init)/2
	limbig<-c(limbig[1]-shift.big,limbig[2]+shift.big)

	return(list(limsmall=limsmall,limbig=limbig,stepsize=stepsize))
}

grid.from.kernel<-function(known.x,known.y,known.z,
	Kernel=expKernel,f=NULL,T=1,
	xlim=NULL,
	ylim=NULL,
	steps=NULL,
	tr=NULL
	){
	if(is.null(tr)){
		tr<-min(max(known.x)-min(known.x),max(known.y)-min(known.y))/(steps/2);
	}
	if(is.null(xlim)){
		xlim=c(min(known.x)-tr,max(known.x)+tr);
	}
	if(is.null(ylim)){
		ylim=c(min(known.y)-tr,max(known.y)+tr);
	}
	if(is.null(f)){
		f<-tr/4
	}


	# get ToGuess locations
	if(abs(xlim[2]-xlim[1])>abs(ylim[2]-ylim[1])){
		out<-adjust.lim(ylim,xlim,steps)
		ylim<-out$limsmall
		xlim<-out$limbig
		stepsize<-out$stepsize
	}else{
		out<-adjust.lim(xlim,ylim,steps)
		xlim<-out$limsmall
		ylim<-out$limbig
		stepsize<-out$stepsize
	}
	# vectors with the xs and ys of the grid
	xs<-seq(xlim[1],xlim[2],stepsize)
	ys<-seq(ylim[1],ylim[2],stepsize)

	# coordonates for all each point
	ToGuess.x<-rep(xs,length(ys))
	ToGuess.y<-as.vector(sapply(ys,rep,length(xs)))

	# get distance matrix
	matdist<-nearest.dist(x=cbind(ToGuess.x,ToGuess.y),y=cbind(known.x,known.y),method="euclidian",delta=tr,upper=NULL)
	
	weightsKnownInToGuessRaw<-Kernel(T,matdist,f); # raw weights

	# get normalized by ToGuess weights
	sumR<-drop(weightsKnownInToGuessRaw%*%rep(1,dim(weightsKnownInToGuessRaw)[2]))
	isolated<-which(sumR<0.01)
	sumRsimple<-sumR
	sumRsimple[isolated]<-1
	NormMat<-diag.spam(1/sumRsimple)
	weightsKnownInToGuess<-NormMat%*%weightsKnownInToGuessRaw

	ToGuess.z<-weightsKnownInToGuess%*%known.z
	ToGuess.z[isolated]<- NA

	return(list(x=ToGuess.x,y=ToGuess.y,z=ToGuess.z,xs=xs,ys=ys,
			raw.weights=weightsKnownInToGuessRaw,dists=matdist
			));
}

num.per.level<-function(v){
	return(aggregate(rep(1,length(v)),by=list(v),sum))
}

zgenHighLevel<-function(est.detection,zNA,est.Q=NULL,est.Kv=NULL,est.Ku=NULL,est.mu=NULL,est.c.comp=NULL,est.v=NULL,est.u=NULL,est.w=NULL,force.mu=FALSE){
	dimension<-length(est.detection)
	u.p<-rep(0,dimension)
	c.p<-rep(0,dimension)
	v.p<-rep(0,dimension)
	if(is.null(est.w)){
		if(is.null(est.u)){
			if(!is.null(est.Ku) && !is.null(est.Q) && !is.null(est.mu)){
				u.p <- drop(rmvnorm.prec.pseudo(n=1, mu=rep(est.mu,dimension), Q=est.Ku*est.Q));
				if(force.mu){
					u.p <-(u.p-mean(u.p)+est.mu)
				}
				cat("mean u.p",mean(u.p),"sd u.p",sd(u.p),"\n")
			}else{
				cat("Use no u component\n");
			}
		}else{
			u.p<-est.u
		}
		if(!is.null(est.c.comp)){
			c.p<-est.c.comp
			cat("mean c.p",mean(est.c.comp),"sd c.p",sd(drop(est.c.comp)),"\n")
		}else{
			cat("Use no c component\n");
		}
		if(is.null(est.v)){
			if(!is.null(est.Kv)){
				Qvp<-diag.spam(est.Kv,dimension)
				if(is.null(est.v)){
					est.v<-rep(0,dimension)
				}
				v.p <- drop(rmvnorm.prec.pseudo(n=1, mu=est.v, Q=Qvp));
				cat("length v.p",length(v.p),"(",dimension,") mean v.p",mean(v.p),"sd v.p",sd(v.p),"\n")
			}else{
				cat("Use no v component\n");
			}
		}else{
			v.p+est.v
		}
		w.p<-u.p+c.p+v.p;
	}else{
		w.p<-est.w
	}

	y.p <- rnorm(n=dimension, mean=w.p, sd=1);
	cat("mean y.p",mean(y.p),"sd y.p",sd(y.p),"\n")
	z.p <-generate_z(y.p,est.detection,zNA)
	cat("mean z.p",mean(z.p[!is.na(z.p)]),"\n")
	return(list(z=z.p,y=y.p,w=w.p,v=v.p,c=c.p,u=u.p))
}

generated.morans.struct<-function(distances,mats_neigh,nbRep,est.detection=NULL,est.Q=NULL,est.Ku=NULL,est.Kv=NULL,est.mu=NULL,est.c.comp=NULL,est.v=NULL,true.val=NULL,est.u=NULL,est.w=NULL,force.mu=FALSE,trueStatus=FALSE){
	if(is.null(est.detection)){
		dimension<-max(c(length(est.u),length(est.v),length(est.w),dim(est.Q)[1]))

		est.detection<-rep(1,dimension)
	}
	if(is.null(true.val)){
		sel<-1:length(est.detection)
		zNA<-c()
	}else{
		sel<-(true.val!=9)
		zNA<-which(true.val==9)
	}
	IdMoinsIs<-mat.or.vec(nbRep,length(distances)-1)
	MI1<-mat.or.vec(nbRep,length(distances)-1)
	MI2<-mat.or.vec(nbRep,length(distances)-1)
	MI3<-mat.or.vec(nbRep,length(distances)-1)
	for(i in 1:nbRep){
		cat(i,"th generation\n")
		generated<-zgenHighLevel(est.detection,zNA,est.w=est.w,est.u=est.u,est.v=est.v,est.Q=est.Q,est.Ku=est.Ku,est.Kv=est.Kv,est.mu=est.mu,est.c.comp=est.c.comp,force.mu=force.mu);
		if(trueStatus){
			y.p<-generated$y
			cat("l y.p:",length(y.p));
			mIrefPseudoObs<-structured.moransI(distances,mats_neigh,y.p,nb_rep_sign=0);
		}else{
			z.p<-generated$z
			mIrefPseudoObs<-structured.moransI(distances,mats_neigh,z.p[sel],nb_rep_sign=0);
		}
		plotIval<-plot.structured.moransI(distances,mIrefPseudoObs,plot=FALSE);
		IdMoinsIs[i,]<-plotIval$mI2-plotIval$mI3
		MI2[i,]<-plotIval$mI2
		MI3[i,]<-plotIval$mI3
		MI1[i,]<-plotIval$mI1
	}
	# meanIdMoinsIs<-apply(IdMoinsIs,2,mean)
	# limIdMoinsIs<-apply(IdMoinsIs,2,quantile,prob=c(0.025,0.975))

	boxplot.free(MI1,ylim=c(0,max(MI1)),xaxt="n")
	plotTrueIval<-NULL
	if(!is.null(true.val)){
		mIrefData<-structured.moransI(distances,mats_neigh,true.val[sel]);
		plotTrueIval<-plot.structured.moransI(distances,mIrefData,plot=FALSE);
		trueIdMoinsIs<-plotTrueIval$mI2-plotTrueIval$mI3
		lines(plotTrueIval$mI1,col=4)
	}
	get.med.position.axis(distances)
	boxplot.free(IdMoinsIs,ylim=c(0,max(IdMoinsIs)),xaxt="n")
	if(!is.null(true.val)){
		lines(trueIdMoinsIs,col=4)
	}
	get.med.position.axis(distances)
	return(list(MI1=MI1,MI2=MI2,MI3=MI3,ref=plotTrueIval))
}

replot.gen.struct<-function(distances,MIs,ylim1=NULL,ylim2=NULL){
	MI1<-MIs$MI1
	MI2<-MIs$MI2
	MI3<-MIs$MI3
	if(is.null(ylim1)){
		ylim1<-c(0,max(MI1))
	}
	boxplot.free(MI1,ylim=ylim1,xaxt="n")
	if(!is.null(MIs$ref)){
		lines(MIs$ref$mI1,col=1)
	}
	dumb<-get.med.position.axis(distances)
	IdMoinsIs<-MI2-MI3
	if(is.null(ylim2)){
		ylim2<-c(0,max(IdMoinsIs))
	}
	boxplot.free(IdMoinsIs,ylim=ylim2,xaxt="n")
	if(!is.null(MIs$ref)){
		lines(MIs$ref$mI2-MIs$ref$mI3,col=1)
	}
	dumb<-get.med.position.axis(distances)
	return()
}

MeanCovPairAtDist<-function(mask,CovMat){
	linkedCovMat<-mask*CovMat
	linkedCovMat<-as.spam(linkedCovMat)
	meanCovPair<-mean(linkedCovMat@entries)

	return(meanCovPair)
}

StructCorrel<-function(distances,Q=NULL,CovMat=NULL,mats_neigh=NULL){
	if(is.null(CovMat)){
		if(!is.null(Q)){
			CovMat<-solve(Q)
		}else{
			stop("Need Q or CovMat\n")
		}
	}
	if(length(mats_neigh[[2]])==3){
		include_streets_anal<-TRUE
	}else{
		include_streets_anal<-FALSE
	}

	SC <- data.frame() ; # list of moran test results general
	for (i in 2:(length(distances))){
		limiteinf=distances[i-1];
		limitesup=distances[i];
		cat("\nlimiteinf=",limiteinf,"limitesup=",limitesup);
		# NB: any function can be put in the glist argument of nb2listw

		label<-paste(limiteinf,limitesup,sep="-")
		dmtr<-mats_neigh[[i]][[1]]
		SC[i,1]<-MeanCovPairAtDist(dmtr,CovMat);
		if(include_streets_anal){
			SBr<-mats_neigh[[i]][[2]]
			ASr<-mats_neigh[[i]][[3]]
			SC[i,2]<-MeanCovPairAtDist(SBr,CovMat);
			SC[i,3]<-MeanCovPairAtDist(ASr,CovMat);
		}
		# print(SC)
	}

	return(SC)
}
StructNeigh<-function(distances,mats_neigh){
	
	nb_mat<-length(mats_neigh[[2]])
	nb_neigh<-as.data.frame(mat.or.vec(length(distances),nb_mat+1))
	if(nb_mat==3){
		include_streets_anal<-TRUE
		names(nb_neigh)<-c("nbNeighTotal","med_position","nbNeighSB","nbNeighAS")
	}else{
		include_streets_anal<-FALSE
		names(nb_neigh)<-c("nbNeighTotal","med_position")
	}

	SC <- data.frame() ; # list of moran test results general
	for (i in 2:(length(distances))){
		limiteinf=distances[i-1];
		limitesup=distances[i];
		cat("\nlimiteinf=",limiteinf,"limitesup=",limitesup);
		# NB: any function can be put in the glist argument of nb2listw

		label<-paste(limiteinf,limitesup,sep="-")
		row.names(nb_neigh)[i]<-label
		dmtr<-mats_neigh[[i]][[1]]
		nb_neigh[i,1]<-sum(dmtr/2) # div by 2 as the matrix is sym
		nb_neigh[i,2]<- (limiteinf+limitesup)/2
		if(include_streets_anal){
			SBr<-mats_neigh[[i]][[2]]
			ASr<-mats_neigh[[i]][[3]]
			nb_neigh[i,3]<-sum(SBr/2)
			nb_neigh[i,4]<-sum(ASr/2)
		}
		# print(SC)
	}

	return(nb_neigh[-1,])
}

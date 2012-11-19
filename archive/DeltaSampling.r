#### Sampling Delta
# the effort induced by streets in terms of pseudo distance
# 1 step of metropolis in the gibbs sampler
# using a gaussian proposal and a gaussian prior for Delta

source("morans_functions.r")
# source("/home/cbarbu/Documents/outils/programmation/R/morans_functions.r")

# pseudo log likelihood for a multivariate normal vector r given it's center b and it's precision matrix Q
dmvnorm.canonical<-function(r,b,Q,logout=TRUE,cholQ=NULL,...){
	if(is.null(cholQ)){
		cholQ <- chol.spam(Q,...);
	}
	# log.det<- -2*determinant(cholQ)$modulus;
	log.det<- determinant(cholQ)$modulus*2
	mu_r <- drop(solve.spam(cholQ, b)) 
	# print("mu_r:")
	# print(mu_r)
	X<-(r-mu_r);
	# print("X:")
	# print(X)
	# log.det<-sum(log(diag(cholQ)))
	in_exp<-drop(t(X)%*%Q%*%X);

	d <- ncol(Q)
	logpi<-d*log(2*pi);

	cat("in_exp:",in_exp,"log.det",log.det,"logpi",logpi,"\n");
	lnP <- -1/2 * (in_exp - log.det+ logpi);

	if (logout==TRUE){
		return(lnP);
	}else{
		return(exp(lnP));
	}
}
dmvnorm.prec<-function(r,mu,Q,logout=TRUE,...){
	# recicling chol and tweaking it (see help(chol)) could still save a lot of time
	# WARNING return the log of the likelihood
	cholQ <- chol.spam(Q,...);
	# log.det<- -2*determinant(cholQ)$modulus;
	# log.det<- determinant.spam(Q,Rstruct=cholQ)$modulus
	log.det<- determinant(cholQ)$modulus*2
	X<-(r-mu);
	# log.det<-sum(log(diag(cholQ)))
	in_exp<-drop(t(X)%*%Q%*%X);

	d <- ncol(Q)
	logpi<- d*log(2*pi);

	# cat("in_exp:",in_exp,"log.det",log.det,"logpi",logpi,"\n");
	lnP <- -1/2 * (in_exp - log.det+ logpi);

	if (logout==TRUE){
		return(lnP);
	}else{
		return(exp(lnP));
	}
}
pseudo_inv<-function(A){
	A.svd<-svd(A);
	Ddiag<-as.spam(diag(A.svd$d));
	A.plus<-A.svd$v%*%(1/Ddiag)%*%t(A.svd$u)
	return(A.plus);
}

importOk<-try(dyn.load("filter_spam.so"),silent=TRUE)

# change all entries in A bigger than maxtobesetnull to 0
if(class(importOk)!="try-error"){
filter.spam<-function(A,maxtobesetnull){
	out <- .C("filter_spam",
		  entries = as.numeric(A@entries),
		  nbent = as.integer(length(A@entries)),
		  limit = as.numeric(maxtobesetnull))
	A@entries<-out$entries
	return(as.spam(A));
}
}else{
	cat("\nWARNING filter_spam.so not available, will use pure R code, may be slower.\n")
	filter.spam<-function(A,maxtobesetnull){
		A@entries[A@entries>maxtobesetnull]<- 0
		return(A);
	}
}

# multiply by f the entries of dist_mat that correspond to an entry in AS
if(class(importOk)!="try-error"){
SpecificMultiply.spam<-function(f,dist_mat,AS){
	out <- .C("specific_multiply",
		  nbent_dm = as.integer(dist_mat@dimension[1]),
		  entries_dm = as.numeric(dist_mat@entries),
		  nbent_dm = as.integer(length(dist_mat@entries)),
		  rpoint_dm = as.integer(dist_mat@rowpointers),
		  colind_dm = as.integer(dist_mat@colindices),
		  entries_AS = as.numeric(AS@entries),
		  nbent_AS = as.integer(length(AS@entries)),
		  rpoint_AS = as.integer(AS@rowpointers),
		  colind_AS = as.integer(AS@colindices),
		  f = as.numeric(f))
	Dmat<-dist_mat;
	Dmat@entries<-out$entries_dm;

	return(Dmat);
}
}else{
	SpecificMultiply.spam<-function(f,dist_mat,AS){
	#SpecificMultiply<-function(f,dist_mat,AS){
		return(dist_mat+AS*dist_mat*(f-1))
	}
}

Qfromf<-function(Dmat,f=1,addeps=epsilon,origin=NULL,kern=default_kern,T=T.r){
	## for an inverse relationship
	if(kern=="inv"){
		if(!is.null(origin)){
			Dmat@entries <- Dmat@entries+origin
		}
		# to have an 1/(D)^(1/f) relationship and the house itself to have a non null weight
		Dmat@entries <- (Dmat@entries)^(f)
		if(is.null(origin)){
			diag.spam(Dmat)<-1
		}
		Q<- -1/Dmat;#the precision matrix of u 
		if(is.null(origin)){
			diag.spam(Q) <- 0;
		}
	}else if(kern=="exp"){
		Q<- -exp(-f*Dmat);
		Q<- (T*AS+SB)*Q;#the precision matrix of u 
	}
	rsumQ <- -Q %*% rep(1,dimension); # the precision for each yi
	
	rsumQ<-rsumQ+addeps 
	norm_fact<- mean(1/rsu+addepsmQ) # the global variance of y
	cat("norm_fact (mean)",norm_fact,"\n");
	diag.spam(Q) <- rsumQ;

	# diag.spam(Q) <- rsumQ+addeps/norm_fact;
	Q<-Q*norm_fact;

	return(Q);
}
QfromfT<-function(Dmat,AS,SB,f=f.r,T=T.r,addeps=epsilon,origin=NULL,kern=Kernel,dimension=dim(Dmat)[1]){
	if(is.spam(SB)){
		Q<- SpecificMultiply.spam(1/T,-kern(T,Dmat,f),SB);#the precision matrix of u 
	}else{
		if(SB==1){
			Q<- -kern(1,Dmat,f)
		}else{
			stop("Bad SB in QfromfT\n")
		}
	}
	# if(length(diagDmat)==dim(AS)[1]){
	# 	Q@entries[diagDmat] <- rep(0,length(diagDmat)) # diagonal starts at 0
	# }else{
		diag(Q)<-0 # diagonal starts at 0
	# }

	rsumQ <- -Q %*% rep(1,dimension); # the precision for each yi
	rsumQ<-rsumQ+addeps # NB we want addeps to be in the normalisation so that we don't have too heavy outliers, addeps is the variance of isolated houses and logically is included in the normalization of the isolation component
	if(use.NormQ == "mean"){
		norm_fact<- mean(1/rsumQ) # the global variance of y
		cat("norm_fact (mean)",norm_fact,"\n");
	}else if(use.NormQ=="med"){
		norm_fact<- median(1/rsumQ) # the global variance of y
		cat("norm_fact (med)",norm_fact,"\n");
	}else{
		norm_fact<-1
	}
	# if(length(diagDmat)==dim(AS)[1]){
	# 	Q@entries[diagDmat] <- rsumQ;
	# }else{
		diag(Q) <- rsumQ;
	# }
	Q<-Q*norm_fact;

	return(Q);
}

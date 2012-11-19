source("pseudo_data_generation.r")
## conclusion: f is identifiable for small values of f (<=0.5), higher it is at risk of tending to +Inf
## increasing Ku.r and Kv.r helps to get higher values of f to be identifiable but still difficult
## b in this model should therefore no be identified 

## u being identifiable we test the simplest model
# w=u
# u~N(0,Q)
# so y~N(u,1)
# and u|y ~ N( (I+Q)^-1 . y,(I+Q)) 
# we want to be sure that for a generated y, the u got with various f 
# are not equivalent in terms of LLH
# so we :
# - with Delta.r -> u.r -> y.r
# - sample u from y.r with various Q(Delta)
# - check Delta.r is at the max

sample_us <- function(n,dimension,Q,K,y,cholR=NULL){
  center <- y;
  R <- Q+diag.spam(1,dimension);
  us <- rmvnorm.canonical(n=n, b=center, Q=R, Rstruct=cholR); # exploit the stability of the structure of R, at least 10 times faster.
  # cholR has to be computed at the initiation with the same structure
  return(us);
}

scan()

tx_passs<-c(0.001,0.005,0.01,0.05,0.1,0.5,1)
nbrepet<-5;
nb_us<-5;

i.r=1
cholQ<-list();
cholR<-list();
det_part<-list();
Q<-list();
for(tx_pass.r in tx_passs){
	nb.tests<-length(tx_passs);

	par(mfrow=c(1,1))
	LLH<-mat.or.vec(nbrepet,nb.tests);
	LLH_us<-mat.or.vec(nb_us,nb.tests);
	for(repet in 1:nbrepet){
		Q.r<-QfromDelta(Delta.r,dist_mat,AS,f=f.r,addeps=epsilon,tx_pass=tx_pass.r);
		if(repet==1){
			cholQ.r<-chol(Ku.r*Q.r);
		}
		u.r <- drop(rmvnorm.prec(n=1, mu=rep(0,dimension), Q=Q.r));
		u.r<-u.r-mean(u.r)+mu.r
		y.r<-drop(rmvnorm.prec(n=1,mu=u.r,Q=diag.spam(1,dimension)))
		# plot_reel(data$easting,data$northing,u.r,main=paste("Delta.r",Delta.r))

		# scan()

		i=1;
		for(tx_pass in tx_passs){
			Delta=Delta.r;
			if(repet==1 && i.r==1){
				Q[i]<-QfromDelta(Delta,dist_mat,AS,f=f.r,addeps=epsilon,tx_pass=tx_pass);
				cholQ[i]<-chol(Ku.r*Q[[i]]);
				det_part[[i]]<- 2*determinant(cholQ[[i]])$modulus
				cholR[i]<-chol(Q[[i]]+diag.spam(1,dimension));
				# this makes a huge difference, specially for big sizes Q
				# in the sampler we should use as much as possible saved values of Q/detQ
				# grid sampling ? 
			}else{
				cat("use saved det_part and Q\n")
			}
			# sample several u using the gibbs sampler
			us<-sample_us(nb_us,dimension,Q[[i]],K,y.r,cholR=cholR[[i]]);

			# get the LLH for each u
			cat("exp_parts:")
			sum_exp_part=0;
			for(k in 1:nb_us){
				u<-us[k,];
				exp_part<- t(u)%*%Q[[i]]%*%u; 
				LLH_us[k,i]= -1/2*(exp_part-det_part[[i]]); 
				sum_exp_part=sum_exp_part+exp_part;
				cat(exp_part,"; ");
			}
			exp_part=sum_exp_part/nb_us;
			cat("\n");
			LLH[repet,i]= -1/2*(exp_part-det_part[[i]]); 
			# as wrong as it can seem, this almost works for threshold = 150

			cat("kernel:",default_kern,"f:",f.r,"exp_part",exp_part,"det_part",det_part[[i]],"LLH:",LLH[repet,i],"\n");
			i=i+1;
		}

		if(repet==1){
			plot(log(tx_passs),LLH[repet,],type="l",main=paste("LLH of u.r fn of tx_pass for tx_pass.r",tx_pass.r," for f.r=",f.r,sep=""),xlab="log(tx_pass)",ylab="Log Likelihood")
			abline(v=log(tx_pass.r))
		}else{
			lines(log(tx_passs),LLH[repet,])
		}
	}
	par(mfrow=c(1,2))
	plot(c(min(log(tx_passs)),max(log(tx_passs))),c(min(LLH),max(LLH)),type="n",main=paste("LLH of u fn of tx_pass for various y.r, tx_pass.r=",tx_pass.r,"f.r=",f.r,sep=""),xlab="log(tx_pass)",ylab="Log Likelihood")
	abline(v=log(tx_pass.r))
	for(repet in 1:nbrepet){
		lines(log(tx_passs),LLH[repet,])
	}
	plot(c(min(log(tx_passs)),max(log(tx_passs))),c(max(LLH)-500,max(LLH)),type="n",main=paste("LLH of u fn of tx_pass for various y.r, f.r=",f.r,sep=""),xlab="log(tx_pass)",ylab="Log Likelihood")
	abline(v=log(tx_pass.r))
	for(num_u in 1:nbrepet){
		lines(log(tx_passs),LLH[num_u,])
	}
	dev.print(device=pdf,file=paste("identif_tx_pass_in_y_tr_",threshold,"_",default_kern,"_f.r_",f.r,"tx",tx_pass.r,".pdf",sep=""))
	i.r=i.r+1;
}




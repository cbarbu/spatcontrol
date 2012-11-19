source("pseudo_data_generation.r")

## u being identifiable we test the simplest model
# w=u
# u~N(0,Q)
# so y~N(u,1)
# and u|y ~ N( (I+Q)^-1 . y,(I+Q)) 
# we want to be sure that for a generated y, the u got with various Delta 
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

Deltas<-seq(0,threshold,10)
nbrepet<-5;
nb_us<-5;

i.r=1
cholQ<-list();
cholR<-list();
det_part<-list();
Q<-list();
for(Delta.r in seq(0,threshold,20)){
	par(mfrow=c(1,1))
	LLH<-mat.or.vec(nbrepet,length(Deltas));
	LLH_us<-mat.or.vec(nb_us,length(Deltas));
	for(repet in 1:nbrepet){
		Q.r<-Ku.r*QfromDelta(Delta.r,dist_mat,AS,f=f.r,addeps=epsilon);

		if(repet==1){
			cholQ.r<-chol(Q.r);
		}
		u.r <- drop(rmvnorm.prec(n=1, mu=rep(0,dimension), Q=Q.r));
		u.r<-u.r-mean(u.r)+mu.r
		y.r<-drop(rmvnorm.prec(n=1,mu=u.r,Q=diag.spam(1,dimension)))
		# plot_reel(data$easting,data$northing,u.r,main=paste("Delta.r",Delta.r))

		# scan()

		i=1;
		for(Delta in Deltas){
			if(repet==1 && i.r==1){
				Q[i]<-Ku.r*QfromDelta(Delta,dist_mat,AS,f=f.r,addeps=epsilon);
				cholQ[i]<-chol(Q[[i]]);
				det_part[[i]]<- 2*determinant(cholQ[[i]])$modulus
				cholR[i]<-chol(Q[[i]]+diag.spam(1,dimension));
				# this makes a huge difference, specially for big sizes Q
				# in the sampler we should use as much as possible saved values of Q/detQ
				# grid sampling ? 
			}else{
				# cat("use saved det_part and Q\n")
			}
			# sample several u using the gibbs sampler
			us<-sample_us(nb_us,dimension,Q[[i]],K,y.r,cholR=cholR[[i]]);

			# get the LLH for each u
			sum_exp_part=0;
			for(k in 1:nb_us){
				u<-us[k,];
				exp_part<- t(u)%*%Q[[i]]%*%u; 
				LLH_us[k,i]= -1/2*(exp_part-det_part[[i]]); 
				sum_exp_part=sum_exp_part+exp_part;
				# cat(exp_part,"; ");
			}
			exp_part=sum_exp_part/nb_us;
			# cat("\n");
			LLH[repet,i]= -1/2*(exp_part-det_part[[i]]); 
			# as wrong as it can seem, this almost works for threshold = 150

			cat("Delta:",Delta,"exp_part",exp_part,"det_part",det_part[[i]],"LLH:",LLH[repet,i],"\n");
			i=i+1;
		}

		if(repet==1){
			plot(Deltas,LLH[repet,],type="l",main=paste("LLH of u.r fn of Delta for Delta.r=",Delta.r,sep=""),xlab="Delta (m)",ylab="Log Likelihood")
			abline(v=Delta.r)
		}else{
			lines(Deltas,LLH[repet,])
		}
	}
	par(mfrow=c(1,2))
	plot(c(min(Deltas),max(Deltas)),c(min(LLH),max(LLH)),type="n",main=paste("LLH of u fn of Delta for various y.r, Delta.r=",Delta.r,sep=""),xlab="Delta (m)",ylab="Log Likelihood")
	abline(v=Delta.r)
	for(repet in 1:nbrepet){
		lines(Deltas,LLH[repet,])
	}
	plot(c(min(Deltas),max(Deltas)),c(max(LLH)-500,max(LLH)),type="n",main=paste("LLH of u fn of Delta same y.r, Delta.r=",Delta.r,sep=""),xlab="Delta (m)",ylab="Log Likelihood")
	abline(v=Delta.r)
	for(repet in 1:nbrepet){
		lines(Deltas,LLH[repet,])
	}
	dev.print(device=pdf,file=paste("identif_y_tr_",threshold,"_",default_kern,"_f_",f.r,"_Delta.r_",Delta.r,".pdf",sep=""))
	i.r=i.r+1;
}




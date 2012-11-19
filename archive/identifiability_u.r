source("pseudo_data_generation.r")
## first test of llh of u.r direclty with various delta

Deltas<-seq(0,threshold,10)
scan()

i.r=1
cholQ<-list();
det_part<-list();
Q<-list();
epsilon=0.1 ## the highest the epsilon the lowest the right part of the llh
#for(Delta.r in seq(0,threshold,20)){
	par(mfrow=c(1,1))
	nbrepet<-15;
	LLH<-mat.or.vec(nbrepet,length(Deltas));
	for(repet in 1:nbrepet){
		Q.r<-QfromDelta(Delta.r,dist_mat,AS,f=f.r,addeps=epsilon);
		if(repet==1){
			cholQ.r<-chol(Ku.r*Q.r);
		}
		u.r <- drop(rmvnorm.prec.pseudo(n=1, mu=rep(mu.r,dimension), Q=Ku.r*Q.r,Rstruct=cholQ.r));
		u.r<-u.r-mean(u.r)+mu.r
		par(mfrow=c(1,2))
		plot_reel(data$easting,data$northing,u.r,main=paste("Delta.r",Delta.r))
		plot(u.r);

		i=1;
		for(Delta in Deltas){
			if(repet==1 && i.r==1){
				Q[i]<-QfromDelta(Delta,dist_mat,AS,f=f.r,addeps=epsilon);
				cholQ[i]<-chol(Q[[i]]);
				det_part[[i]]<- 2*determinant(cholQ[[i]])$modulus
				# this makes a huge difference, specially for big sizes Q
				# in the sampler we should use as much as possible saved values of Q/detQ
				# grid sampling ? 
			}else{
				cat("use saved det_part and Q\n")
			}
			exp_part<- t(u.r)%*%Q[[i]]%*%u.r; 
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
	plot(c(min(Deltas),max(Deltas)),c(min(LLH),max(LLH)),type="n",main=paste("LLH of u.r fn of Delta for Delta.r=",Delta.r,sep=""),xlab="Delta (m)",ylab="Log Likelihood")
	abline(v=Delta.r)
	for(repet in 1:nbrepet){
		lines(Deltas,LLH[repet,])
	}
	dev.print(device=pdf,file=paste("identif_u_tr_",threshold,"_Delta.r_",Delta.r,".pdf",sep=""))
	i.r=i.r+1;
# }


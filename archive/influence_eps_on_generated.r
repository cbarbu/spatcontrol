source("pseudo_data_generation.r")

Delta.r=100;
par(mfrow=c(2,3));
for(epsilon in c(0.00001,0.001,0.01,0.1,1,10)){
	Q.r<-QfromDelta(Delta.r,dist_mat,AS,f=f.r,addeps=epsilon);
	if(repet==1){
		cholQ.r<-chol(Ku.r*Q.r);
	}
	u.r <- drop(rmvnorm.prec(n=1, mu=rep(0,dimension), Q=Q.r));
	u.r<-u.r-mean(u.r)+mu.r
	# y.r<-drop(rmvnorm.prec(n=1,mu=u.r,Q=diag.spam(1,dimension)))
	plot_reel(data$easting,data$northing,u.r,main=paste("epsilon:",epsilon))
}



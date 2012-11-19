
plot.prob.to.observed<-function(Prob,Obs,nb_classes=10,xlab="Predicted",ylab="Observed",...){
	nb_classes<-10
	mid_class<-rep(0,nb_classes)
	prop_pos<-rep(0,nb_classes)
	tot_nb_pos<-rep(0,nb_classes)
	size_class<-rep(0,nb_classes)
	for(i in 1:nb_classes){
		min_prob<-(i-1)/nb_classes
		max_prob<-i/nb_classes
		mid_class[i]<-(min_prob+max_prob)/2
		cat("class",min_prob,"to",max_prob,":")
		risk_group<-which(Prob<max_prob & Prob>=min_prob)
		nb_pos<-sum(Obs[risk_group])
		nb_total<-length(Obs[risk_group])
		prop_pos[i]<-nb_pos/nb_total
		tot_nb_pos[i]<-nb_pos;
		cat(prop_pos[i],"of",nb_total,"\n")
		size_class[i]<-nb_total
	}
	plot(mid_class,prop_pos,xlab=xlab,ylab=ylab,...)
	abline(a=0,b=1)

	return(list(mid_class,prop_pos,tot_nb_pos,size_class))
}
par(mfrow=c(2,2))
Qmean<-QfromfT(Dmat,AS,SB,f=meanf,T=meanT)
spatPrec<-diag(Qmean)
u_pred<-pnorm(est.u,0,1+sqrt(1/(meanKu/spatPrec)))
plot.prob.to.observed(u_pred,visudata)
hist(u_pred)

# adjusting for inspectors
est.detection<-inspector%*%est.beta
plot.prob.to.observed(u_pred*est.detection,visudata)
hist(u_pred*est.detection)
printdev(device=pdf,"fit_by_spat.pdf")

c_pred<-pnorm(est.c.val,0,1)
w_pred<-pnorm(est.w,0,1)
final_proba<-w_pred*est.detection
par(mfrow=c(1,4))
plot_reel(db$easting, db$northing,visudata,base=0,top=1,main="data")
plot_reel(db$easting, db$northing,u_pred,base=-0,top=1,main="spatial proba")
plot_reel(db$easting, db$northing,c.map%*%c_pred,base=-0,top=1,main="cofactors proba")
plot_reel(db$easting, db$northing,final_proba,base=0,top=1,main="final proba")
printdev(device=pdf,"maps_prob.pdf")


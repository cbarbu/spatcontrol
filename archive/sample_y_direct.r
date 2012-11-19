# enable to sample a mixture of two truncated normal with upper limit 
# of one being the lower limit of the other and this limit being 0
# this allow specifically to sample y, the "continuous reality" in the probit model
# directly according to the data

library(msm)
library("truncnorm")

composite_ptnorm <- function(y,x,bi){
	area_inf<-pnorm(0,mean=x,sd=1)
	A<-area_inf;
	B<-(1-area_inf)*(1-bi);

	ret<-rep(0,length(y));
	ret[y<0]<-pnorm(y[y<0],mean=x,sd=1)/(A+B);
	ret[y>=0]<-(pnorm(y[y>=0],mean=x,sd=1)*(1-bi)+area_inf*bi)/(A+B);

	return(ret);
}

# to pick something according to the previous curve
sample_composite_ptnorm_vect <- function(xvect,bivect){
	# sample y given that it's density is a normalized sum of 
	# dnorm(xvect,1,0,+Inf)*(1-beta)+dnorm(xvect,1,-Inf,0)
	# xvect: probability of being 1 vs. 0
	# bivect: probability of being observed as 1 vs. 0
	l<-length(xvect);
	# cat("l",l);
	tunifvect<-runif(l);
	A<-pnorm(0,mean=xvect,sd=1); 	# P(y-,z-|w)=P(y-|w) as P(y-,z+)=0
	B<-(1-A)*(1-bivect);		# P(y+,z-|w)=P(y+,obs-|w)=P(obs-|y+,w)*P(y+|w)=(1-beta)*(1-P(y-|w))
	area<-A+B			# P(z-)=P(y-|w)+P(y+,obs-|w)
	samp<-B/area;			# P(y+|z-,w)=P(y+,z-|w)/P(z-)
	samp[area==0]<-1 ; 		# avoid errors when bivect=1 and xvect>0

	tfinal<-mat.or.vec(l,1);
	ypos<-which(tunifvect<=samp);
	yneg<-which(tunifvect>samp);
	# tfinal[ypos]<- rtnorm(length(ypos),mean=xvect[ypos],sd=1,lower=0,upper=Inf); # yprime=1, determine a corresponding y
	# tfinal[yneg]<- rtnorm(length(yneg),mean=xvect[yneg],sd=1,lower=-Inf,upper=0); # yprime=0, determine a corresponpding y
	# cat("\n",length(tfinal[tunifvect<=samp]),length(which(tunifvect<=samp)),length(xvect[which(tunifvect<=samp)]));
	
	if(length(ypos)>0){
		tfinal[ypos]<- rtruncnorm(1,mean=xvect[ypos],sd=1,a=0,b=Inf); # yprime=1, determine a corresponding y
	}

	if(length(yneg)>0){
		tfinal[yneg]<- rtruncnorm(1,mean=xvect[yneg],sd=1,a=-Inf,b=0); # yprime=0, determine a corresponpding y
	}
	
	return(tfinal)
}

# from there we can sample y which is sampled differently according to the observation or not of bugs in the data
sample_y_direct <-function(w,zpos,zneg,zNA,bivect){
	# return the continuous variable result of the probit model
	# it is the augmented variable of yprime
	# y prime can simply be obtained by using
	# yprime<- (y>0)
	# w is the general risk in the probit model
	# zpos the vector of references of positive z
	# zneg --------------------------- negative z
	# zNA ---------------------------- unknown z
	# bivect the vector of per house probability to find infestation if there is

	# y[zpos]<- rtnorm(length(zpos),mean=w[zpos],sd=1,lower=0,upper=Inf);
	y<-0*w
	if(length(zpos)>0){
		y[zpos]<- rtruncnorm(1,mean=w[zpos],sd=1,a=0,b=Inf);
	}

	if(length(zneg)>0){
		y[zneg]<- sample_composite_ptnorm_vect(w[zneg],bivect[zneg]);
		## draw first the inspectors, stupid
		# lookWell<-which(rbinom(length(bivect[zneg]),1,prob=bivect[zneg])>0.5)
		# wneg<-w[zneg]
		# yneg<-0*wneg
		# if(length(yneg[-lookWell])>0){
		# 	yneg[-lookWell]<- rtruncnorm(1,mean=wneg[-lookWell],sd=1,a=0,b=Inf);
		# }
		# if(length(lookWell)>0){
		# 	yneg[lookWell]<- rtruncnorm(1,mean=wneg[lookWell],sd=1,a=-Inf,b=0);
		# }
		# cat("look Well:",length(lookWell),"out of",length(zneg),"\n")
		# y[zneg]<-yneg
	}

	# y[zNA] <- rnorm(length(zNA),mean=w[zNA],sd=1) # w[zNA];
	if(length(zNA)>0){
		y[zNA] <- rnorm(length(zNA),mean=w[zNA],sd=1) # w[zNA];
	}

	return(y);
}

dmtnorm<-function(x,mu=0,std=1,shift=0,w1=1,w2=1,logout=F){
	# x the value(s) to examin
	# shift the values at wich we shift from the first to the second trunc normal
	# w1/w2 the weight of the left/right part 
	# output the probability as a log
	A<-pnorm(shift,mean=mu,sd=std)*w1;
	B<-(1-pnorm(shift,mean=mu,sd=std))*w2;
	x1<-x;
	x1[x>=shift]<-Inf;
	x2<-x;
	x2[x<shift]<--Inf;
	px=dnorm(x1,mean=mu,sd=std)*w1/(A+B)+dnorm(x2,mean=mu,sd=std)*w2/(A+B);
	# cat("A",A,"B",B,"px",px,"\n");
	if(logout==T){
		px<-log(px);
	}
	return(px);
}
area.perso<-function(fun,lower,upper,eps=1e-05,...){
	# numerical integral of fun between lower and upper
	x<-seq(lower,upper,eps);
	dens<-fun(x,...);
	a<-(sum(dens)+(dens[2]-dens[length(dens)])/2)*eps;
	return(a);
}

# # testing:
# bi= 0.7 # find rate of inspectors
# x=-0.3
# A<-pnorm(0,mean=x,sd=1);
# B<-(1-pnorm(0,mean=x,sd=1))*(1-bi);
# xabs<-seq(-5,5,0.01)
# xabs1<-seq(-5,0,0.01)
# xabs2<-seq(0,5,0.01)
# 
# par(mfrow=c(1,6))
# plot(c(xabs1,xabs2),dnorm(c(xabs1,xabs2),mean=x,sd=1),type="n")
# lines(xabs1,dnorm(xabs1,mean=x,sd=1))
# lines(xabs2,dnorm(xabs2,mean=x,sd=1)*(1-bi))
# 
# plot(xabs,dmtnorm(xabs,mu=x,std=1,shift=0,w2=(1-bi)))
# 
# plot(c(xabs1,xabs2),pnorm(c(xabs1,xabs2),mean=x,sd=1),type="n")
# lines(xabs1,pnorm(xabs1,mean=x,sd=1))
# lines(xabs2,pnorm(xabs2,mean=x,sd=1)*(1-bi))
# 
# plot(c(xabs1,xabs2),pnorm(c(xabs1,xabs2),mean=x,sd=1),type="n")
# lines(xabs1,pnorm(xabs1,mean=x,sd=1)/(A+B))
# lines(xabs2,(pnorm(xabs2,mean=x,sd=1)*(1-bi)+pnorm(0,mean=x,sd=1)*bi)/(A+B))
# 
# plot(xabs,composite_ptnorm(xabs,x,bi))
# 
# xvect<-rep(x,10000)
# bivect<-rep(bi,10000)
# system.time(y<-sample_composite_ptnorm_vect(xvect,bivect))
# 
# hist(y,freq=FALSE)
# lines(xabs1,dnorm(xabs1,mean=x,sd=1)/(A+B))
# lines(xabs2,dnorm(xabs2,mean=x,sd=1)*(1-bi)/(A+B))
# 
# # very fast (less than 2times rtnorm for 10000) and exact

# test sample_y_direct

# nbItem<-10000
# w<-rnorm(nbItem,0,3)
# bivect<-rep(0.99,nbItem)
# zpos<-sample(1:nbItem,round(20*nbItem/100))
# zNoPos<-(1:nbItem)[-zpos]
# zNApos<-sample(1:length(zNoPos),round(20*nbItem/100))
# zNA<-zNoPos[zNApos]
# zneg<-zNoPos[-zNApos]
# 
# y<-sample_y_direct(w,zpos,zneg,zNA,bivect)
# par(mfrow=c(1,3))
# plot(y[zpos],w[zpos])
# abline(a=0,b=1)
# plot(y[zNA],w[zNA])
# abline(a=0,b=1)
# plot(y[zneg],w[zneg])
# abline(a=0,b=1)

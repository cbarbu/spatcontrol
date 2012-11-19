## parameters
Delta.r <- 100; # <- rtnorm(1,mean=muDelta,sd=sdDelta,lower=0,upper=Inf); # the supplementary distance due to streets in meters
threshold <- 150; # max number of meters to consider neighbourgs
f.r<- 1 # the strength of the relationship on long distances

## get necessary functions
source("/home/cbarbu/Documents/outils/programmation/R/zoom.r")
source("spam_complement.r")
library("spdep");
calculate_accross_streets_inside_blocks_correlograms<-function(distances,pres_abs,nb_sb,nb_as){
	#### calculate the spatial autocorrelogram
	# pres_abs:  presence absence data 
	# nb_sb: list of neighbourgs of the same block
	# nb_as: list of neighbourgs accross streets
	mI2 <- list() ; # list of moran test results in blocks
	mI3 <- list() ; # list of moran test results across streets
	for (d in 2:(length(distances))){
		print("d:");
		print(d);
		limiteinf=distances[d-1];
		limitesup=distances[d];
		cat("limiteinf=",limiteinf,"limitesup=",limitesup,"\n");

		dnb2<-nb_sb[[d-1]]
		dnb3<-nb_as[[d-1]]
		# nb alternatively glist can be used to give a customized list of weights
		temp=system.time(lw2 <- nb2listw(dnb2,zero.policy=TRUE));
		print(temp)
		temp=system.time(lw3 <- nb2listw(dnb3,zero.policy=TRUE));
		print(temp)
		cat("2 listw ok temp=",temp,"\n");
		temp=system.time(mt2 <- moran.test(pres_abs, lw2, zero.policy=TRUE,adjust.n=FALSE));
		print(temp)
		print(mt2)
		temp=system.time(mt3 <- moran.test(pres_abs, lw3, zero.policy=TRUE,adjust.n=FALSE));
		print(temp)
		print(mt3)
		cat("moran ok temp=",temp,"\n");
		mI2[[paste(limiteinf,limitesup,sep="-")]]<-mt2 # save it using as key the limite as a string
		mI3[[paste(limiteinf,limitesup,sep="-")]]<-mt3 # save it using as key the limite as a string
	}
	return(list(mI2,mI3));
}
dist_mat_for_moran_perso<-function(data,threshold){
	data_non_NA<-data[data$status!=9,]
	dim_mp<-length(data_non_NA[,1])

	dmtmp <-nearest.dist(x=data_non_NA[,c("easting","northing")], y=NULL, method="euclidian", delta=threshold, upper=NULL);          
	dmtmp@entries<-rep(1,length(dmtmp@entries))# [dmt@entries!=0]<-1 # 1 only when dist_mat not 0
	SBmp<- nearest.dist(x=cbind(data_non_NA$block_num,rep(0,dim_mp)), method="euclidian", upper=NULL,delta=0.1)
	ASmp<-dmtmp-SBmp; # get 1 whereever the distances matrix is defined(under threshold) and not same block
	ASmp<-as.spam(ASmp)
}
moran_perso<-function(dmt,values){
	# dmt a spam matrix with the weights between houses
	# values the values for each location
	# easily 20 and probably up to 100 times faster than moran.test

	# correction of dmt for correct calculus of I
	dim_mp<-length(data_non_NA[,1])
	one_vect<-rep(1,dim_mp)
	#need 0 on the diagonal of weight
	diag(dmt)<-0
	#need weight normalized by row
	rsum<-dmt%*%one_vect
	rsum[rsum==0]<-1;
	rsum_mat<-dmt*0;
	diag(rsum_mat)<-1/rsum
	dmt<-rsum_mat%*%dmt

	# matrix of residuals
	res_z<- values-mean(values)
	res_z_mat<-diag.spam(1,length(values));
	diag.spam(res_z_mat)<-res_z

	### Moran's I
	# covariance term
	int<-dmt%*%res_z_mat;
	int<-res_z_mat%*%int

	MI<-length(res_z_mat)*one_vect%*%int%*%one_vect/(sum(dmt@entries)*res_z%*%res_z)
	return(drop(MI))
}
plot_correlogram_from_mI2mI3<-function(distances,name,mI2mI3){
	#### extract values ok for plotting

	mI2<-mI2mI3[[1]];
	mI3<-mI2mI3[[2]];
	morans_I2 = {};
	morans_I2_ref = {}
	morans_I2_error = {}
	morans_I3 = {}
	morans_I3_ref = {}
	morans_I3_error = {}
	med_position = {}
	legend_position = {};
	for (i in 1:length(mI2)){
		med_position=c(med_position,distances[i]-eval(parse(text=(names(mI2[i]))))/2);
		legend_position=c(legend_position,names(mI2[i]));
		morans_I2 = c(morans_I2,mI2[[i]]$estimate[[1]])
		morans_I2_ref = c(morans_I2_ref,mI2[[i]]$estimate[[2]])
		morans_I2_error = c(morans_I2_error,sqrt(mI2[[i]]$estimate[[3]])*1.96);	
		morans_I3 = c(morans_I3,mI3[[i]]$estimate[[1]])
		morans_I3_ref = c(morans_I3_ref,mI3[[i]]$estimate[[2]])
		morans_I3_error = c(morans_I3_error,sqrt(mI3[[i]]$estimate[[3]])*1.96);	
	}
	cat("ready to plot");
	# #### display the spatial autocorrelogram
	dev.new()
	plot(c(distances[1],max(distances)),c(0,1),type='n',xaxt='n',xlab="distances (m)",ylab="Morans's I",zs=0,main=name)
	errbar(med_position,morans_I2,morans_I2-morans_I2_error,morans_I2+morans_I2_error,add=TRUE,col=2);
	errbar(med_position,morans_I3,morans_I3-morans_I3_error,morans_I3+morans_I3_error,add=TRUE,col=4);
	axis(1,at=med_position,labels=legend_position);
}
sample_composite_ptnorm_vect <- function(xvect,bivect){
	# sample y given that it's density is a normalized sum of 
	# dnorm(xvect,1,0,+Inf)*(1-beta)+dnorm(x,1,-Inf,0)
	l<-length(xvect);
	tunifvect<-runif(l);
	A<-pnorm(0,mean=xvect,sd=1)
	B<-(1-A)*(1-bivect);
	area<-A+B
	samp<-A/area;

	tfinal<-mat.or.vec(l,1);
	tfinal[tunifvect>=samp]<-rtnorm(length(which(tunifvect>=samp)),mean=xvect,sd=1,lower=0,upper=Inf);
	tfinal[tunifvect<samp]<- rtnorm(length(which(tunifvect<samp)),mean=xvect,sd=1,lower=-Inf,upper=0);
	
	return(tfinal)
}
sample_y_direct <-function(x,zpos,zneg,zNA,bivect){
	# return the continuous variable result of the probit model
	# it is the augmented variable of yprime
	# y prime can simply be obtained by using
	# yprime<- (y>0)

	y[zpos]<- rtnorm(length(zpos),mean=x,sd=1,lower=0,upper=Inf);
	y[zneg]<- sample_composite_ptnorm_vect(x[zneg],bivect[zneg]);
	y[zNA] <- rnorm(length(zNA),mean=x,sd=1)

	return(y);
}

## get clean data
data.base <- read.csv("DB_mm_blocks_05May2011.csv",header=T);
# data<-data.base[data.base$northing>8185000&data.base$northing<8185500&data.base$easting>233250&data.base$easting<233750,]
data<-data.base[data.base$northing&data.base$northing&data.base$easting&data.base$easting,]
library("spdep");
dnb <- dnearneigh(as.matrix(data[,c("easting","northing")]), 0, threshold)
isolated<-which(card(dnb)<1)

if(length(isolated)>0){
	data<-data[-isolated,];
}
# prep data for autocorrelation 
## the problem is that we have a lot of NA in the status for true DATA
## we should then compute the autocorrelation only on houses on which we have true data 
## to be able to compare directly the things while we apply the model to all the houses

houses_XY2<- data[c("easting","northing")]; # not unicode and number of bugs
houses_XY2<- houses_XY2[data$status != 9,]
blocks=data$block_num;
blocks=blocks[data$status != 9]
truedata<-data$status[data$status != 9]

# get the list of neighbourgs, same block or accross streets
nb_sb<-list()
nb_as<-list()
distances = c(0,20,40,60,80,100,120)
nb_houses<-length(blocks)
for (d in 2:(length(distances))){
	limiteinf=distances[d-1];
	limitesup=distances[d];
	cat("limiteinf=",limiteinf,"limitesup=",limitesup,"\n");
	temp=system.time(dnb <- dnearneigh(as.matrix(houses_XY2), limiteinf, limitesup));
	cat("near neigh ok temp=",temp,"\n");
	with_neigh<-1:nb_houses
	with_neigh<-with_neigh[card(dnb)>0];
	# NB: any function can be put in the glist argument of nb2listw
	dnb2 <-dnb # list of neighbourghs in the same block
	dnb3 <-dnb # list of neighbourghs not in the same block
	for(i in with_neigh){
		block=blocks[i];
		blocks_neigh=blocks[dnb[[i]]];

		val_voisins2<-dnb[[i]][which(blocks_neigh==block)];
		val_voisins3<-dnb[[i]][which(blocks_neigh!=block)];

		dnb2[[i]] <- as.integer(val_voisins2)
		dnb3[[i]] <- as.integer(val_voisins3)
	}
	nb_sb[[d-1]]<-dnb2;
	nb_as[[d-1]]<-dnb3;
}
graphics.off()
#### calculate autocorrelogram with true data
mI2mI3<-calculate_accross_streets_inside_blocks_correlograms(distances,truedata,nb_sb,nb_as);
plot_correlogram_from_mI2mI3(distances,paste("real data"),mI2mI3);

#### get matrix of distances
## matrix of basic distances
spam.options(nearestdistnnz=c(9058076,400))
dist_mat <-nearest.dist(x=data[,c("easting","northing")], y=NULL, method="euclidian", delta=threshold, upper=NULL);          
cat("got dist_mat");

#### generation of a spatial mean:
dimension<-length(data$status)
Id<-diag.spam(1,dimension) # matrix identity
ratio_pos<-length(which(data$status==1))/length(which(data$status!=9))
mu_data<-qnorm(ratio_pos,mean=0,sd=1)# the mean 

#### get spatial variance in the data
## generate pseudo y corresponding to the data: yd

yd<-drop(sample_y_direct(rep(mu_data,dimension),which(data$status==1),which(data$status==0),which(data$status==9),rep(1,dimension)));
wd<-drop(rmvnorm.prec(1,mu=yb,Q=Id))

#### generation of the base random vectors
y0<-rnorm(dimension,mu_data,sd=1)
w0<-drop(rmvnorm.prec(1,mu=y0,Q=Id))
norm_res<-rnorm(dimension,0,sd=1)


## prep to get matrix of distances corrected by delta
source("DeltaSampling.r") # allow to get AS correct
cat("got AS");

## for various Delta get the generated data and corresponding autocorrelogram
library(Hmisc);
for(Delta in c(0,50,100)){
	cat("Delta:",Delta);
	####  Qw : the matrices with weights non normalized (spatial sum without normalization for the spatial mean)
	## inverse relationship
	Qw <- - QfromDelta(Delta,dist_mat,AS,f=f.r,addeps=0,origin=NULL);
	diag(Qw)<-0; # the value for the very house is not taken into account 
	# NB: this is forced by the fact that we use a kernel in inverse of the distance but could be changed later

	#### Qn the matrix that get the spatial mean
	# matSum : the matrix with the sum by row on the diagonal to be able to get the spatial mean
	sbyrow<-rsum(Qw)
	isolated<-which(sbyrow==0)

	# to allow calculus of inverse after
	sbyrow[isolated]<- 1 
	 
	matSum<-0*Qw;
	diag(matSum)<-sbyrow
	# matNorm : the matrix normalizing Qw, exist as we took care to remove isolated houses
	matNorm<-0*Qw
	matNorm<-1/matSum
	Qm<-matNorm%*%Qw
	# Qn : the matrix given the spatial mean, so that mean(Qn%*%w) = mean(w)
	Qn<-matNorm%*%Qw
	diag(Qn)[isolated]<-1 # add isolated houses into the spatial term so that the variance 
	# that need to be added can be correctly evaluated

	## get spatial term
	u<-Qn%*%w0
	# that can be multiplied to get the observed spatial part in the data with this matrix
	ud<-Qn%*%wd
	uf<-drop(sd(ud)/sd(u)*u)

	## get precision matrix of residuals
	Vr<-drop(var(wd)-var(uf))
	norm_fact<-sum(diag(matSum))/dimension
	norm_cov_mat_res<-matSum/norm_fact;
	res<-sqrt(Vr)*norm_cov_mat_res%*%norm_res
	
	#### get generated data
	w<-uf+res; 
	# z<-as.real(u>0); # for simplist model 
	## with this model, u is too low with not enough dispersion, no positives
	y<-rmvnorm.spam(1,mu=w,Sigma=Id) 
	z<-as.real(w>0); # for base model
	# with this model, the lack of dispersion for u makes not possible to have any impact
	# of the streets on the output, at least not significative

	#### plot generated data
	dev.new()
	par(mfrow=c(1,3))
	plot_reel(data$easting,data$northing,2*uf,main=paste("u for Delta=",Delta))
	plot_reel(data$easting,data$northing,2*y,main=paste("y for Delta=",Delta))
	plot_reel(data$easting,data$northing,2*(z-0.5),main=paste("z for Delta=",Delta))

	
	mI2mI3<-calculate_accross_streets_inside_blocks_correlograms(distances,z[data$status!=9],nb_sb,nb_as);
	cat("mI2mI3 generated");
	
	plot_correlogram_from_mI2mI3(distances,paste("Delta=",Delta),mI2mI3);
}

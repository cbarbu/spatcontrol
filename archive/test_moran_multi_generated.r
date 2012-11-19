# generate and test for Moran's correlogram multiple maps
# the difficulty is to generate on the full map but test the moran's I only on the nonNA

graphics.off()
## parameters
source("parameters_sampler.r") # default to the values for the simulation
set.seed(200)
Ku.r <- 0.418# rgamma(n=1, shape=K.hyper[1], scale=K.hyper[2]);
Kv.r <- 212.97; # rgamma(n=1, shape=K.hyper[3], scale=K.hyper[4]);
threshold <- 50; # max number of meters to consider neighbourgs
f.r<- 14.98 # the strength of the relationship on long distances 0.5 favorise extremely 
# or simply permit the identification of Delta
T.r=0.264 # taux d'association accross streets over association within blocks
epsilon=0.01 # the value added to the diagonal of Q in order to allow for LLH calculations and data generation
use.NormQ<-"mean"
# need to be small enough to allow for the identication of Delta
mu.r<- -1.48
mv.r<- 0
beta.r<-0.7; # in testing, should use the real beta or at least something well distributed
default_kern="exp"
use.v.gen<-TRUE;

period<- "nofall.2007" # keep the all map, for the test

cat("u.r for mu.r",mu.r,", Ku.r",Ku.r,", tr",threshold,", f.r:",f.r,"\n")

## parameters 
# distances = c(0,20,40,60,80,100,120,140) # delimitation of the rings
# distances = seq(5,605,20)
distances = seq(5,125,20)
include_streets_anal<- TRUE # avoid when large distances
check_significance<-FALSE # check the significance of the diff AS/SB (takes time)
nb_rep_sign<-0
nbsimul=3;

# not to be changed
use.map.gen<-FALSE

## functions
library("spam")
# library("spdep");
source("spam_complement.r")
source("morans_functions.r")
source("zoom.r")
source("functions_intercept.r")

## data
GLOBALSETPARAMETERS<-FALSE
source("pseudo_data_generation.r")
zNA <- which(data$status==9);
nbNA<-length(zNA)

if(check_significance){
	dyn.load("useC2.so")
}

med_position = {}

dimension <- nrow(data);

K.r<-c(Ku.r,Kv.r)

Qv.r<-diag.spam(Kv.r,nrow(Q.r))
cholQv.r<-chol(Qv.r);
bivect<-rep(beta.r,dimension) # should be improved

## prep data for autocorrelation study
housestemp<-data
houses_XY <- housestemp[,c("easting","northing")]; 
values <- housestemp$status # get presence absence data 
values[values==9]<-NA
blocks <- housestemp$block_num;

## general plot of the data
dev.new()
par(mfrow=c(1,1))
sizepoint<-0.5
sel<- which(values==0)
plot(housestemp$easting[sel],housestemp$northing[sel],cex=sizepoint,pch=16,asp=1,col=3,xlab="easting (m)",ylab="northing (m)")
sel<- which(is.na(values))
lines(housestemp$easting[sel],housestemp$northing[sel],type="p",cex=sizepoint,pch=16,asp=1,col=4)
sel<- which(values==1)
lines(housestemp$easting[sel],housestemp$northing[sel],type="p",cex=sizepoint,pch=16,asp=1,col=2)

#### prep neighbourhoods
## define rings of neighbours general, same block and accross streets
cat("get neighbourgs\n");
# need to recalculate distances, for houses above threshold in Q

noNA<-which(!is.na(values))
mats_neigh<-gen.mats.neigh(distances,data$easting[noNA],data$northing[noNA],data$block_num[noNA])

#### calculus and plot autocorrelation in the original data
cat("true data autocorrelation\n")
mIref<-structured.moransI(distances,mats_neigh,values[noNA],nb_rep_sign=nb_rep_sign);
dev.new()
par(mfrow=c(1,2))
simple_mIref<-plot.structured.moransI(distances,mIref);

#### main loop: generation of "nbsimul" data and autocorrelation calculation
dev.new()
par(mfrow=c(1,2))
IdMoinsIsFull<-generated.morans.struct(distances,mats_neigh,nbsimul,est.Q=Q.r,est.Ku=Ku.r,est.mu=mu.r,est.Kv=Kv.r,bivect,true.val=data$status)
dev.new()
table.mI<-{}
saved.seeds<-round(runif(nbsimul,min=0,max=16000))
cat("start loop")
for(num_simul in 1:nbsimul){
	cat("num simul:",num_simul,"out of:",nbsimul,"\n")
	# set.seed(saved.seeds[num_simul])
	## generation of the pseudo data
	z.r<-zgenHighLevel(bivect,zNA,est.Q=Q.r,est.Ku=Ku.r,est.Kv=Kv.r,est.mu=mu.r,est.v=rep(0,dimension),est.c.comp=rep(0,dimension))
	values<-z.r
	# u.r <- drop(rmvnorm.prec.pseudo(n=1, mu=rep(mu.r,dimension), Q=Ku.r*Q.r,Rstruct=cholQ.r));
	# u.r<-u.r-mean(u.r)+mu.r

	# w.r<-2*u.r
	# if(use.v.gen){
	# 	v.r <- drop(rmvnorm.prec.pseudo(n=1, mu=rep(mv.r,dimension), Q=Qv.r,Rstruct=cholQv.r));
	# 	# v.r <- drop(rmvnorm.canonical(n=1, b=rep(mv.r,dimension), Q=Qv.r));
	# 	w.r<-w.r+v.r
	# }
	# y.r <- rnorm(n=dimension, mean=w.r, sd=1);

	# z.r<-generate_z(y.r,bivect,zNA); 
	# values<-z.r
	# values[z.r==9]<- NA

	## plot of the generated data
	par(mfrow=c(2,2))
	# nameplot <- paste("u.r for mu.r",mu.r,", Ku.r",Ku.r,", tr",threshold,", f.r:",f.r,")
	# plot_reel(data$easting,data$northing,u.r,main=nameplot)
	# if(use.v.gen){
	# 	name <- paste("v.r for Ku.r",Ku.r,"Kv.r",Kv.r,"mv.r",mv.r)
	# 	plot_reel(data$easting,data$northing,v.r,main=name)
	# }
	# nameplot <- paste("y.r for Ku.r",Ku.r,"Kv.r",Kv.r,"mu.r",mu.r,"mv.r",mv.r)
	# plot_reel(data$easting,data$northing,y.r,main=nameplot)
	# nameplot <- paste("z.r for Ku.r",Ku.r,"Kv.r",Kv.r,"mu.r",mu.r,"mv.r",mv.r)
	mI<-structured.moransI(distances,mats_neigh,values[noNA]);

	plot_reel(data$easting[noNA],data$northing[noNA],2*(data$status[noNA]-1),main="data")
	plot_reel(data$easting[noNA],data$northing[noNA],2*(z.r[noNA]-1),main="generated data")

	# dev.set(dev.prev())
	plot.structured.moransI(distances,mIref);
	simple_mI<-plot.structured.moransI(distances,mI);
	table.mI<-rbind(table.mI,c(simple_mI$mI1,simple_mI$mI2,simple_mI$mI3))
	# dev.set(dev.next())

	## plot the difference between mI2/mI3
	if(check_significance){
		difft<-mI[[4]]
		## plot diff AS/SB
		diffr<-morans_I2-morans_I3
		plot(c(min(med_position),max(med_position)),c(min(difft),max(diffr)),type="n",xlab="distance(m)",ylab="Isb - Ias")
		lines(med_position,diffr,col=4,type="b")
		## get the 95% confidence lines 
		ordered_difft<-mat.or.vec(length(distances),nb_rep_sign)
		for(i in 2:length(distances)){
			ordered_difft[i,]<-difft[i,order(difft[i,])]
		}
		dump("ordered_difft",file="ordered_difft.txt")
		boxplot(med_position,t(difft[-1,]),add=TRUE)
	}
}

#### plotting of the results
## plot nb pair neighbours
# count
nb_neigh<-mI[[5]]
med_position<-(distances[2:length(distances)]+distances[1:(length(distances)-1)])/2
plot(c(min(med_position),max(med_position)),c(min(nb_neigh),max(nb_neigh)),type="n")
lines(med_position,nb_neigh[-1,1],type="l",col=1)
lines(med_position,nb_neigh[-1,2],col=4)
lines(med_position,nb_neigh[-1,3],col=2)

# proportions
plot(c(min(med_position),max(med_position)),c(0,1),type="n",xlab="distance",ylab="rate of pairs of neighbours inside same block")
lines(med_position,nb_neigh[-1,2]/nb_neigh[-1,1])

# summary of autocorrelation
dev.new()
nbclasses<-length(distances)-1
par(mfcol=c(2,3))
simple_mIref<-plot.structured.moransI(distances,mIref);
morans_I2ref<-simple_mIref$mI2;
morans_I3ref<-simple_mIref$mI3;
plot(med_position,morans_I2ref-morans_I3ref,type="l");

mean.table.mI<-apply(table.mI,2,mean);
morans_I1mean<-mean.table.mI[1:nbclasses];
morans_I2mean<-mean.table.mI[(1:nbclasses)+nbclasses]
morans_I3mean<-mean.table.mI[(1:nbclasses)+2*nbclasses]

diff.mI<-table.mI[,(1:nbclasses)+nbclasses]-table.mI[,(1:nbclasses)+2*nbclasses]
mean.diff.mI<-apply(diff.mI,2,mean);
plot(c(distances[1],distances[(length(distances))]),c(0.8*min(morans_I1mean,morans_I2mean,morans_I3mean),1.1*max(morans_I1mean,morans_I2mean,morans_I3mean)),type='n',xaxt='n',xlab="distances (m)",ylab="Morans's I")
axis(1,at=med_position,labels=legend_position)
lines(med_position,morans_I1mean,col=1) # black general
if(include_streets_anal){
	lines(med_position,morans_I2mean,col=4) # blue within blocks
	lines(med_position,morans_I3mean,col=2) # red inter blocks
}
plot(mean.diff.mI,type="l",main="mean diff mI inter - intra blocks")
# dev.print(device=pdf,"autocorrelation_compared_to_mean_gen.pdf")

### plot one of the top of the other
# autocorrelation
simple_mIref<-plot.structured.moransI(distances,mIref);
lines(med_position,morans_I1mean,col=1,lty=2) # black general
if(include_streets_anal){
	lines(med_position,morans_I2mean,col=4,lty=2) # blue within blocks
	lines(med_position,morans_I3mean,col=2,lty=2) # red inter blocks
}

# difference
plot(med_position,morans_I2ref-morans_I3ref,type="l");
lines(med_position,mean.diff.mI,type="l",lty=2)

### plot of ref above boxplots
# only general autocorrelation
# here boxplot plot the quartiles: 25%, median, 75% + whiskers for the 95% CI

dev.new(width=3.5,height=3.5)
op <- par(mar = c(4,4,0.5,0.5))
boxplot.free(table.mI[,1:6],ylim=c(0,max(c(table.mI,simple_mIref$mI1))),xaxt="n",xlab="distance (m)",ylab=expression("Morans'I"))
axis(1,at=1:6,labels=legend_position)
lines(simple_mIref$mI1,col=4)
# dev.print(device=pdf,"MoransI_generated.pdf")

# only diff
dev.new(width=3.5,height=3.5)
op <- par(mar = c(4,4,0.5,0.5))
diff.table<-table.mI[,7:12]-table.mI[,13:18]
diff.ref<-morans_I2ref-morans_I3ref
boxplot.free(diff.table,ylim=c(min(c(diff.table,diff.ref)),max(c(diff.table,diff.ref))),xaxt="n",xlab="distance (m)",ylab=expression(I[S]-I[D]),range=0)
axis(1,at=1:6,labels=legend_position)
lines(diff.ref,type="l",col=4);
# dev.print(device=pdf,"Diff_IS_ID_generated.pdf")

# ### study the fits
# # get the fits
# fit_mI1=mat.or.vec(nbclasses,1)
# for(i in 1:nbsimul){
# 	simple_mI<-table.mI[i,]
# 	fit_mI1[i]<-sum((simple_mI[1:nbclasses]-simple_mIref$mI1)^2);
# }
# # display the trials from worse to best
# dev.new()
# par(mfrow=c(2,4))
# plot_reel(data$easting,data$northing,data$status,main="data");
# simple_mIref<-plot.structured.moransI(distances,mIref);
# for(best in 1:3){
# 	set.seed(saved.seeds[order(fit_mI1)[best]])
# 	cat("num best:",best,"fit:",fit_mI1[order(fit_mI1)[best]],"\n");
# 	u.r <- drop(rmvnorm.prec.pseudo(n=1, mu=rep(mu.r,dimension), Q=Q.r,Rstruct=cholQ.r));
# 	u.r<-u.r-mean(u.r)+mu.r
# 	if(use.v.gen){
# 		v.r <- drop(rmvnorm.prec.pseudo(n=1, mu=rep(mv.r,dimension), Q=Qv.r,Rstruct=cholQv.r));
# 		# v.r <- drop(rmvnorm.canonical(n=1, b=rep(mv.r,dimension), Q=Qv.r));
# 		w.r<-u.r+v.r
# 	}else{
# 		w.r<-u.r
# 	}
# 	y.r <- rnorm(n=dimension, mean=w.r, sd=1);
# 
# 	bivect<-rep(beta.r,dimension) # can be improved
# 	z.r<-generate_z(y.r,bivect,zNA); 
# 
# 	# generate NA
# 	zNA.r<-sample(1:dimension,nbNA)
# 	z.r[zNA.r]<-9
# 
# 	plot_reel(data$easting,data$northing,z.r,main=paste("best fit",best));
# 	simple_mI<-table.mI[order(fit_mI1)[best],]
# 	morans_I1<-simple_mI[1:nbclasses];
# 	morans_I2<-simple_mI[(1:nbclasses)+nbclasses];
# 	morans_I3<-simple_mI[(1:nbclasses)+2*nbclasses];
# 	plot(c(distances[1],distances[(length(distances))]),c(0.8*min(morans_I1,morans_I2,morans_I3),1.1*max(morans_I1,morans_I2,morans_I3)),type='n',xaxt='n',xlab="distances (m)",ylab="Morans's I")
# 	axis(1,at=med_position,labels=legend_position)
# 	lines(med_position,morans_I1,col=1) # black general
# 	if(include_streets_anal){
# 		lines(med_position,morans_I2,col=4) # blue within blocks
# 		lines(med_position,morans_I3,col=2) # red inter blocks
# 	}
# 
# }


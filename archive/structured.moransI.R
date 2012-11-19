source("spam_complement.r")
source("morans_functions.r")
dyn.load("useC2.so")
name<-"MIstruct_Infest_AllPau_0-135by15_100rep"
distances = seq(0,135,15)
use.map.gen<-FALSE
city="PAUCARPATA"
period=""
subsetcity=0
use.struct<-TRUE
nb_rep_sign<-100 # nb repetitions for permutation test
variable<-"infestation" # "infestation" or "opening"

if(variable=="opening"){
				value.NA<-9
				use.NA<-TRUE
}else if(variable=="infestation"){
				value.NA<-9
				use.NA<-FALSE
}else{
				stop("study for this variable not yet implemented")
}

source("import_data.r")

if(variable=="infestation"){
				variableObs<-data$status
}else if(variable=="opening"){
				variableObs<-as.numeric(data$status!=9)
}

# only moran's I
par(mfrow=c(2,3))
mats_neigh<-gen.mats.neigh(distances,data$easting,data$northing)
mIref<-structured.moransI(distances,mats_neigh,variableObs);
plot.structured.moransI(distances,mIref,neigh.struct=TRUE);

if(use.struct){
	mats_neigh<-gen.mats.neigh(distances,data$easting,data$northing,data$block_num)

	# mIref1<-structured.moransI(distances,mats_neigh,variableObs,nb_rep_sign=nb_rep_sign,rand.sym=TRUE);
	# plot.structured.moransI(distances,mIref1);

	mIref2<-structured.moransI(distances,mats_neigh,variableObs,nb_rep_sign=nb_rep_sign,rand.sym=FALSE);
	plot.structured.moransI(distances,mIref2);

	par(mfrow=c(2,2))
	# plot.structured.moransI(distances,mIref1);
	plot.structured.moransI(distances,mIref2);
}

# ### add neighbors structure on morans_I1 alone
# dev.new()
# par(mfrow=c(1,1),mar=c(4,4,4,4))
# plot(c(distances[1],distances[(length(distances))]),c(0.8*min(c(0,morans_I1)),1.1*max(morans_I1)),type='n',xaxt='n',xlab="distances (m)",ylab="Morans's I")
# axis(1,at=med_position,labels=legend_position)
# lines(med_position,morans_I1,col=1,type="b",pch=3) # black general
# rap.neigh<-mIref$nb.neigh[-1,2]/mIref$nb.neigh[-1,1]
# lines(med_position,rap.neigh*max(morans_I1),lty=2,col="blue")
# axis(4,at=seq(0,max(morans_I2),0.25*(max(morans_I2))),labels=seq(0,100,25),col="blue")
# mtext("% of neighbors on same city block", side=4, line=3, cex.lab=1,las=3, col="blue")


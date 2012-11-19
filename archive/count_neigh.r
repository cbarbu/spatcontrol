
if(!exists("GLOBALSETPARAMETERS")){
	source("parameters_sampler.r")
}else{
	if(GLOBALSETPARAMETERS==TRUE){
		source("parameters_sampler.r")
	}
}
subsetcity<-0
source("inspectors_detection.r")
## data
source("spam_complement.r")

if(use.map.gen){
	map<-make_data_points(nb_houses_per_sidex,nb_blocks_per_side,nb_blocks_per_side,inter_houses_space,nb_houses_per_sidey,ratio_street_dist)
	plot(map$easting,map$northing,col=map[,3],pch=15,cex=0.5)
	data<-map
}else{
	if(city=="PAUCARPATA"){
		data.base <- read.csv("DB_simple_Pau_cyclo1_19Jul2011.csv",header=TRUE);

		# data<-data.base[data.base$northing>8181750&data.base$northing<8182800&data.base$easting>232500&data.base$easting<234000,]
		# plot_reel(data$easting,data$northing,2*data$status-1)
			
		## small
		# height_sub<-500;
		# width_sub<-550;

		if(subsetcity==0){
			data<-data.base; 
		}else if(subsetcity==1){
			height_sub<-1000;
			width_sub<-1500;
			xmin<-232500
			ymin<-8181750
			# rect(xmin,ymin,xmin+width_sub,ymin+height_sub,border=4)
		}else if(subsetcity==2){
			height_sub<-1000;
			width_sub<-1500;
			xmin<-233500
			ymin<-8182800
			# rect(xmin,ymin,xmin+width_sub,ymin+height_sub,border=4)
		}else if(subsetcity==3){
			## small data set for testing
			height_sub<-500;
			width_sub<-300;
			xmin<-233000
			ymin<-8182000
		}
		data<-data.base[data.base$northing>ymin&data.base$northing<ymin+height_sub&data.base$easting>xmin&data.base$easting<xmin+width_sub,]
		data$oanimal<-as.integer((data$CO==1 | data$AV==1 | data$GA==1 | data$OV==1 | data$otros.animales!="-1"))

		## two windows
	}else{

		# data.base <- read.csv("DB_mm_blocks_05May2011.csv",header=TRUE);
		# data.base <- read.csv("data.csv",header=TRUE);
		# colnames(data.base) <- c("locality","status","collector","easting","northing");

		#q# subset the data arround transect
		# data<-data.base[data.base$northing>8185250&data.base$northing<8185500&data.base$easting>232000&data.base$easting<232250,]
		# data<-data.base[data.base$northing>8185000&data.base$northing<8185500&data.base$easting>233250&data.base$easting<233750,]
		# data<-data.base[data.base$northing>8185000&data.base$easting<233400,]
		# data<-data.base

		### cleaning of inspectors
		data$collector[data$collector==" Jorge Ampuero"] <- "Jorge Ampuero"
		data$collector[data$collector==" Jorge A"] <- "Jorge Ampuero"
		data$collector[data$collector=="Jose Velasquez "] <- "Jose Velasquez"
		data$collector[data$collector==" Hugo Vilcahuaman"] <- "Hugo Vilcahuaman"
		data$collector[data$collector=="Julio cesar Condori"] <- "Julio Cesar Condori"
		data$collector[data$collector=="Manuel Tamayo "] <- "Manuel Tamayo"
		data$collector<-factor(data$collector)
	}

	if(! (use.NA)){
		data<-data[data$status!=9,];
	}
	# ## remove too isolated houses
	# library("spdep");
	# dnb <- dnearneigh(as.matrix(data[,c("easting","northing")]), 0, threshold)
	# isolated<-which(card(dnb)<3)

	# if(length(isolated)>0){
	# 	data<-data[-isolated,];
	# 	cat("remove:",isolated,"\n");
	# }

}

## functions
spam.options(nearestdistnnz=c(9058076,400))
dist_mat <-nearest.dist(x=data[,c("easting","northing")], y=NULL, method="euclidian", delta=threshold, upper=NULL);          
dist_mat<-as.spam(dist_mat)

nbh<-length(data$easting)

count_inf_to<-function(vect,tr){
	return(length(which(vect<tr)))
}

dist_tr<-seq(5,60,5)
nbtr<-length(dist_tr)
isolated<-rep(0,nbtr)
under2<-rep(0,nbtr)
mean_num_neigh<-rep(0,nbtr)
med_num_neigh<-rep(0,nbtr)


for(num_tr in 1:nbtr){
	isolated[num_tr]<-length(which(apply_by_row_not_null.spam(dist_mat,min)>dist_tr[num_tr]));
	under2[num_tr]<-length(which(apply_by_row_not_null.spam(dist_mat,count_inf_to,dist_tr[num_tr])<2))
	mean_num_neigh[num_tr]<-mean(apply_by_row_not_null.spam(dist_mat,count_inf_to,dist_tr[num_tr]))
	med_num_neigh[num_tr]<-median(apply_by_row_not_null.spam(dist_mat,count_inf_to,dist_tr[num_tr]))
}
par(mfrow=c(2,3))
plot(dist_tr[-(1:3)],isolated[-(1:3)])
plot(dist_tr[-(1:3)],under2[-(1:3)])
plot(dist_tr,log10(isolated))
lines(dist_tr,log10(under2))
plot(dist_tr,isolated/nbh)
lines(dist_tr,under2/nbh)
abline(h=0.05)
abline(h=0.01)
plot(dist_tr,mean_num_neigh)
lines(dist_tr,med_num_neigh)

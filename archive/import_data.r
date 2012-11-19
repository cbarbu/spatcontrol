if(use.map.gen){
	map<-make_data_points(nb_houses_per_sidex,nb_blocks_per_side,nb_blocks_per_side,inter_houses_space,nb_houses_per_sidey,ratio_street_dist)
	plot(map$easting,map$northing,col=map[,3],pch=15,cex=0.5)
	db<-map
}else{
	if(city=="PAUCARPATA"){
		db.base <- read.csv("DB_simple_Pau_cyclo1_19Jul2011_blockSize.csv",header=TRUE);

		db.base$oanimal<-as.integer((db.base$CO==1 | db.base$AV==1 | db.base$GA==1 | db.base$OV==1 | db.base$otros.animales!="-1"))

		db.base$TrueStatus<-db.base$status
		db.base$status[db.base$status==9]<-value.NA

		# db<-db.base[db.base$northing>8181750&db.base$northing<8182800&db.base$easting>232500&db.base$easting<234000,]
		# plot_reel(db$easting,db$northing,2*db$status-1)
			
		## small
		# height_sub<-500;
		# width_sub<-550;

		if(subsetcity==0){
			height_sub<-max(db.base$northing)-min(db.base$northing);
			width_sub<-max(db.base$easting)-min(db.base$easting);
			xmin<-min(db.base$easting)
			ymin<-min(db.base$northing)
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
		# all
		# sel<- 1:dim(db.base)[1]
		# locality based
		loc<-as.integer(t(as.vector(as.data.frame(lapply(as.character(db.base$unicode),strsplit,".",fixed=TRUE))[3,])))
		# fall 2007
		if(period=="fall.2007"){
			sel<-(loc==46 | loc == 48 | loc == 50 | loc == 51 | loc == 52 | loc == 55)
		}else if(period=="fall.2008"){
			# fall 2008
			sel<-(loc==87 | loc==88 | loc == 90 | loc == 92 | loc == 94)
		}else if(period=="nofall.2007"){
			sel<-(loc!=46 & loc != 48 & loc != 50 & loc != 51 & loc != 52 & loc != 55)
		}else{
			sel<-rep(TRUE,length(loc))
		}

		db<-db.base[sel & db.base$northing>=ymin&db.base$northing<=ymin+height_sub&db.base$easting>=xmin&db.base$easting<=xmin+width_sub,]

		## two windows
	}else{

		db.base <- read.csv("DB_mm_blocks_05May2011.csv",header=TRUE);
		# db.base <- read.csv("data.csv",header=TRUE);
		# colnames(db.base) <- c("locality","status","collector","easting","northing");

		#q# subset the data arround transect
		# db<-db.base[db.base$northing>8185250&db.base$northing<8185500&db.base$easting>232000&db.base$easting<232250,]
		# db<-db.base[db.base$northing>8185000&db.base$northing<8185500&db.base$easting>233250&db.base$easting<233750,]
		# db<-db.base[db.base$northing>8185000&db.base$easting<233400,]
		db<-db.base

		### cleaning of inspectors
		db$collector[db$collector==" Jorge Ampuero"] <- "Jorge Ampuero"
		db$collector[db$collector==" Jorge A"] <- "Jorge Ampuero"
		db$collector[db$collector=="Jose Velasquez "] <- "Jose Velasquez"
		db$collector[db$collector==" Hugo Vilcahuaman"] <- "Hugo Vilcahuaman"
		db$collector[db$collector=="Julio cesar Condori"] <- "Julio Cesar Condori"
		db$collector[db$collector=="Manuel Tamayo "] <- "Manuel Tamayo"
		db$collector<-factor(db$collector)
	}

	if(subsetcity==0){
	}else if(subsetcity==1){
		# limits itself to the localities where there is at least 1 infested house
		sel<-which(db$status !=9)
		infestationByLoc<-aggregate(db$status[sel],by=list(db$locality[sel]),sum)
		keepLoc<-which(infestationByLoc$x>0)
		sel<-which(db$locality %in% infestationByLoc$Group.1[keepLoc])
		db<-db[sel,];
	}
	if(! (use.NA)){
		db<-db[db$status!=9,];
	}
	# ## remove too isolated houses
	# library("spdep");
	# dnb <- dnearneigh(as.matrix(db[,c("easting","northing")]), 0, threshold)
	# isolated<-which(card(dnb)<3)

	# if(length(isolated)>0){
	# 	db<-db[-isolated,];
	# 	cat("remove:",isolated,"\n");
	# }

}

## remove duplicates
dup<-which(duplicated(cbind(db$easting,db$northing)))
if(length(dup)>0){
db<-db[-dup,]
}
db$X<-db$easting
db$Y<-db$northing



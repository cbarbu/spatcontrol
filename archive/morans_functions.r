library("spam")
source("functions_intercept.r")
source("R/symmetric.random.entries.by.line.spam.c.R")
### better done by A*B
# logic_and_spam<-function(A,B){
# 	# A and B matrices spam of same dimension
# 	# with 0 or 1 values
# 	# return a matrix C of same dimension with 
# 	# C[i,j]<-1 if (A[i,j] == 1 && B[i,j] == 1)
# 	# C[i,j]<-0 else
# 	C<-A+B;
# 	and<-which(C@entries>1.1)
# 	no_and<-which(C@entries<1.1)
# 	C@entries[and]<-rep(1,length(and))
# 	C@entries[no_and]<-rep(0,length(no_and))
# 	C<-as.spam(C)
# 	return(C)
# }
moran_perso<-function(dmt,values){
	# dmt a spam matrix with the weights between houses
	# values the values for each location
	# easily 20 and probably up to 100 times faster than moran.test

	# correction of dmt for correct calculus of I
	dim_mp<-length(values)
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
gen.mats.neigh<-function(distances,x,y,group=NULL){
	dist_matA <-nearest.dist(x=cbind(x,y), y=NULL, method="euclidian", delta=max(distances), upper=NULL);          

	dmtA<-dist_matA
	dimension<-dim(dist_matA)[1]
	dmtA@entries<-rep(1,length(dmtA@entries))# [dmtA@entries!=0]<-1 # 1 only when dist_matA not 0
	if(!is.null(group)){
		SBA <- nearest.dist(x=cbind(group,rep(0,length(group))), method="euclidian", upper=NULL,delta=0.1)
		SBA@entries<-rep(1,length(SBA@entries))
		SBA@entries<-rep(1,length(SBA@entries))
		SBA<-SBA*dmtA;

		ASA<-dmtA-SBA; # get 1 whereever the distances matrix is defined(under threshold) and not same block
		ASA<-as.spam(ASA)
	}
	mats_neigh<-list()
	for (i in 2:(length(distances))){
		limiteinf=distances[i-1];
		limitesup=distances[i];
		cat("\ndist:",limiteinf,"->",limitesup);

		dmtr<-dmtA # will be the matrix of neighbours for a ring(1/0)
		unselect<-which(dist_matA@entries>=limitesup | dist_matA@entries<limiteinf)
		dmtr@entries[unselect]<-rep(0,length(unselect))# 1 when dist_matA defined (distance<max(distance))
		dmtr<-as.spam(dmtr)
		cat(" mean nb Neigh",round(sum(dmtr/2)/dimension,digits=2))

		if(!is.null(group)){
			# ASr<-logic_and_spam(ASA,dmtr);
			# SBr<-logic_and_spam(SBA,dmtr);
			ASr<-ASA*dmtr;
			SBr<-SBA*dmtr;
			cat(" (SG:",round(sum(SBr/2)/dimension,digits=2),";DG:",round(sum(ASr/2)/dimension,digits=2),")")
			mats_neigh[[i]]<-list(dmtr=dmtr,SBr=SBr,ASr=ASr);
		}else{
			mats_neigh[[i]]<-list(dmtr=dmtr);
		}
	}
	cat("\n")
	return(mats_neigh)
}

# compute structured morans I 
structured.moransI<-function(distances,mats_neigh=NULL,raw.values,nb_rep_sign=0,rm.NA=TRUE,rand.sym=FALSE){
	# distances: breaks for the classes of distance to be tested for Moran's I
	# mats_neigh: for each class of distances the matrices of neighborhood (general, same block, different blocks)
	# values: analysed values
	# include_streets_anal: if compute the structured by streets moran's I
	# nb_rep_sign: number of repetitions if computing the significance
	# rm.NA: should NA be removed/ignored?

	## test validity of values
	zNA<-which(is.na(raw.values))
	zNoNA<-which(!is.na(raw.values))
	if(length(zNA)>0){
		if(!rm.NA){
			stop("values contains NA, if should be ignored, please use rm.NA=TRUE")
		}else{
			values<-raw.values[zNoNA]
		}
	}else{
		values<-raw.values
	}
	if(length(mats_neigh[[2]])==3){
		include_streets_anal<-TRUE
	}else{
		include_streets_anal<-FALSE
	}
	
	## calculate moran's I
	mI1 <- list() ; # list of moran test results general
	mI2 <- list() ; # list of moran test results in blocks
	mI3 <- list() ; # list of moran test results across streets
	difft<-NULL;
	nb_neigh<-mat.or.vec(length(distances),3)

	if(nb_rep_sign>0){
		difft<-mat.or.vec(length(distances),nb_rep_sign)
	}
	for (i in 2:(length(distances))){
		limiteinf=distances[i-1];
		limitesup=distances[i];
		cat("\nlimiteinf=",limiteinf,"limitesup=",limitesup);
		# NB: any function can be put in the glist argument of nb2listw

		label<-paste(limiteinf,limitesup,sep="-")
		if(length(zNA)>0){
			dmtr<-mats_neigh[[i]][[1]][zNoNA,zNoNA]
		}else{
			dmtr<-mats_neigh[[i]][[1]]
		}
		mI1[[label]]<-moran_perso(dmtr,values);
		cat(" Moran's I:", mI1[[label]]);
		nb_neigh[i,1]<-sum(dmtr/2) # div by 2 as the matrix is sym
		if(include_streets_anal){
			if(length(zNA)>0){
				SBrtrue<-mats_neigh[[i]][[2]][zNoNA,zNoNA]

			}else{
				SBrtrue<-mats_neigh[[i]][[2]]
			}
			if(length(zNA)>0){
				ASr<-mats_neigh[[i]][[3]][zNoNA,zNoNA]
			}else{
				ASr<-mats_neigh[[i]][[3]]
			}
			mI2[[label]]<-moran_perso(SBrtrue,values);
			mI3[[label]]<-moran_perso(ASr,values);
			nb_neigh[i,2]<-sum(SBrtrue/2)
			nb_neigh[i,3]<-sum(ASr/2)
			if(nb_rep_sign>0){
				# get a 0 in entries of SBr any time dmtr is 1
				phantom<-dmtr
				phantom@entries<-rep(0,length(phantom@entries))
				SBrtrue<-phantom+as.spam(SBrtrue)
				for(repetition in 1:nb_rep_sign){
					#randomize SBrtrue
					if(rand.sym){
						SBr<-as.spam(symmetric.random.entries.by.line.spam.c(SBrtrue)) # randomize sym
					}else{
						SBr<-as.spam(random_spam_entries_by_row(SBrtrue)) # randomize not sym
					}
					#put complement in ASr
					ASr<-as.spam(dmtr-SBr);
					mI2t<-moran_perso(SBr,values);
					mI3t<-moran_perso(ASr,values);
					difft[i,repetition]<-mI2t-mI3t;
				}
			}
		}
	}
	cat("\n")
	return(list(mI1=mI1,mI2=mI2,mI3=mI3,difft=difft,nb.neigh=nb_neigh))
}
legend_position<-NULL
plot.structured.moransI<-function(distances,mI,add=FALSE,neigh.struct=FALSE,plot=TRUE){
	mI1<-mI[[1]]
	if(length(mI[[2]])>0){
		mI2<-mI[[2]]
		mI3<-mI[[3]]
		include_streets_anal<-TRUE
	}else{
		include_streets_anal<-FALSE
	}

	## plot moran's I
	morans_I1 = {}
	morans_I2 = {}
	morans_I3 = {}
	med_position = {}
	legend_position = {}
	for (i in 1:length(mI1)){
		med_position=c(med_position,distances[i]-eval(parse(text=(names(mI1[i]))))/2);
		legend_position=c(legend_position,names(mI1[i]));
		morans_I1 = c(morans_I1,mI1[[i]])
		if(include_streets_anal){
			morans_I2 = c(morans_I2,mI2[[i]])
			morans_I3 = c(morans_I3,mI3[[i]])
		}
	}

	# pdf(file="moran_fn_dist_intrablock.pdf",width=8.5,height=7);
	if(plot){
		if(add==FALSE){
			plot(c(distances[1],distances[(length(distances))]),c(0.8*min(morans_I1,morans_I2,morans_I3),1.1*max(morans_I1,morans_I2,morans_I3)),type='n',xaxt='n',xlab="distances (m)",ylab="Morans's I")
			axis(1,at=med_position,labels=legend_position)
		}
		lines(med_position,morans_I1,col=1) # black general
		if(include_streets_anal){
			lines(med_position,morans_I2,col=4) # blue within blocks
			lines(med_position,morans_I3,col=2) # red inter blocks
			if(!is.null(mI[[4]])){
				plot.signif.structured.moransI(med_position,mI[[4]],morans_I1,morans_I2,morans_I3,legend_position)
			}
		}
	}
		
	# dev.off();
	legend_position<<-legend_position;
	return(list(mI1=morans_I1,mI2=morans_I2,mI3=morans_I3,med_position=med_position))
}
plot.signif.structured.moransI<-function(med_position,difft,morans_I1,morans_I2,morans_I3,legend_position,main=NULL,pch=3,col=4,xlab="distance (m)",ylab=expression(paste(I[sg],"-",I[dg])),...){
	## plot diff AS/SB
	diffr<-morans_I2-morans_I3
	nb_rep_sign<-dim(difft)[2]
	if(is.null(main)){
		main=paste(nb_rep_sign,"permutations")
	}
	plot(c(min(med_position),max(med_position)),c(min(difft,na.rm=TRUE),max(diffr,na.rm=TRUE)),type="n",xlab=xlab,xaxt="n",ylab=ylab,main=main,pch=pch,...)
	axis(1,at=med_position, lab=legend_position)
	lines(med_position,diffr,col=col,type="b",pch=pch,...)
	# ## get the 95% confidence lines 
	# ordered_difft<-mat.or.vec(length(med_position),nb_rep_sign)
	# for(i in 2:length(med_position)){
	# 	ordered_difft[i,]<-difft[i,order(difft[i,])]
	# }
	boxplot.free(t(difft[-1,]),add=TRUE,xaxt="n",at=med_position,pars=list(boxwex=(med_position[2]-med_position[1])/3),pch=pch)
	
	return()
}
get.med.position.axis<-function(distances){
	left.lim<-distances[1:(length(distances)-1)]
	right.lim<-distances[2:length(distances)]
	med.position<-(right.lim+left.lim)/2
	labels.axis<-paste(as.character(left.lim),"-",as.character(right.lim),sep="")

	axis(1,at=seq(1:length(left.lim)), lab=labels.axis)
	return(med.position)
}

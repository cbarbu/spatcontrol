# make the mean of the fields to allow a generation of data according to it

getField<-function(filename){
	nblinesString<-system(paste("wc -l ",filename,sep=""),intern=TRUE)
	nblinesTable<-strsplit(nblinesString,split=" +")[[1]]
	nonVoid<-which(nzchar(nblinesTable))[1]
	nblines<-as.integer(nblinesTable[nonVoid[1]])
	cat("initially",nblines,"lines\n")

	skiped<-floor(nblines/2);

	nblinesread<-min(1000,nblines-skiped);
	# nblinesread<-nblines-skiped;
	cat("going to import",nblinesread,"lines beginning at",skiped,"...")
	fields<-scan(filename,skip=skiped,sep="\t",nlines=nblinesread); # can take some time (15s of for 750 lignes skipping 750 lines
	cat("Imported\n")
	fields<-matrix(fields,nrow=nblinesread,byrow=TRUE)

	return(fields)
}
if(exists("INTERMEDIARY")){
	if(INTERMEDIARY){
# importOk<-try(source("estimated.txt"),silent=TRUE)
# if(class(importOk)=="try-error"){
	us<-getField("usamples.txt")
	par(mfrow=c(1,2))
	est.u<-apply(us,2,mean)
	sd.u<-apply(us,2,sd)
	plot_reel(db$easting,db$northing,est.u,base=quantile(est.u,prob=c(0.05)),top=quantile(est.u,prob=c(0.95)),main="est.u")
	plot_reel(db$easting,db$northing,sd.u,base=0,top=quantile(sd.u,prob=c(0.95)),main="sd.u")
	rm(us)
	ws<-getField("wsamples.txt")
	est.w<-apply(ws,2,mean)
			w<-ws[dim(ws)[1],]
			rm(ws)
			est.comp<-c.map%*%est.c.val
			est.v<-est.w-est.u-est.comp
		}
# }
}
# source("pseudo_data_generation.r")

par(mfrow=c(2,2))
plot_reel(db$easting,db$northing,visudata*2-1,main="data")
CIest.u<-as.vector(quantile(est.u,probs=c(0.025,0.975)))
plot_reel(db$easting,db$northing,est.u,main="mean u",base=CIest.u[1],top=CIest.u[2]);
plot_reel(db$easting,db$northing,est.c.comp,main="mean c",base=CIest.u[1],top=CIest.u[2])
plot_reel(db$easting,db$northing,est.w,main="mean w",base=CIest.u[1],top=CIest.u[2])
# dump("est.u",file="meanu.r")

# ## for film
# dev.new()
# for(i in 1:nblinesread){
# 	plot_reel(db$easting,db$northing,mus[i,])
# 	# Sys.sleep(1)
# }

## get to minimal, publishable dataset
source("extrapol_field.R")
db<-read.csv("DB_simple_Pau_cyclo1_19Jul2011_blockSize.csv",header=TRUE);
db<-db[order(db$easting,db$northing,db$status,db$collector),]

db$oanimal<-as.integer((db$CO==1 | db$AV==1 | db$GA==1 | db$OV==1 | db$otros.animales!="-1"))

db.small<-db[,c("easting","northing","block_num","CU","PE","oanimal","I.NO","P.NO")]
db.small$open<-as.numeric(db$status != 9)
db.small$infested<-as.numeric(db$status == 1)
db.small$inspector<-as.numeric(db$collector)
loc<-as.integer(t(as.vector(as.data.frame(lapply(as.character(db$unicode),strsplit,".",fixed=TRUE))[3,])))
# better
sel<-(loc==46  | loc == 42| loc == 48 | loc == 50 | loc == 51 | loc == 52 | loc == 53  | loc == 54| loc == 55)

# # original
# sel<-(loc==46 | loc == 48 | loc == 50 | loc == 51 | loc == 52 | loc == 55)
db.small$fitSet<-as.numeric(sel)
db.small<-changeNameCol(db.small,"block_num","cityBlockNum")


## if adjustments needed
plot(db.small$easting,db.small$northing,col=db.small$fitSet+1,asp=1)
with(db.small[which(db.small$infested==1),],lines(easting,northing,col="yellow",type="p"))
locs<-aggregate(db.small[,c("easting","northing")],by=list(loc),mean)
text(locs$easting,locs$northing,locs$Group.1,col=4)


dup<-which(duplicated(cbind(db.small$easting,db.small$northing)))
db.small<-db.small[-dup,]

write.csv(db.small,"OriginalDataPaucarpata.csv",row.names=FALSE)


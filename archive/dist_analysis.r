use.map.gen=FALSE
city="PAUCARPATA"
subsetcity=0
period=""
use.NA=TRUE
library("spam")
source("spam_complement.r")
source("import_data.r")

dist_mat <-nearest.dist(x=data[,c("easting","northing")], y=NULL, method="euclidian", delta=300, upper=NULL);

SB <- nearest.dist(x=cbind(data$block_num,rep(0,length(data$block_num))), method="euclidian", upper=NULL,delta=0.1)
SB@entries<-rep(1,length(SB@entries))
diag(SB)<-0
SB<-as.spam(SB)
SBnoSelf<-SB
diag(SBnoSelf)<-0


maxdistSB<-apply_by_row_not_null.spam(dist_mat*SBnoSelf,max)
mindistSB<-apply_by_row_not_null.spam(dist_mat*SB,min)

par(mfrow=c(1,2))
hist(mindistSB)
hist(maxdistSB)

# should check that this varies 


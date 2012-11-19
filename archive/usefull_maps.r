source("pseudo_data_generation.r")
## years of spraying
dev.new()
plot(data.base$easting,data.base$northing,col=data.base$FR_A-1,cex=0.1,asp=1)
sel<-which(data.base$pos==1)
lines(data.base$easting[sel],data.base$northing[sel],col="red",type="p",pch=16,cex=0.3)
legend("topright",title="Year of spray",levels(as.factor(data.base$FR_A)),col=levels(as.factor(data.base$FR_A-1)),pch=16)
## month of spraying
dev.new()
plot(data.base$easting,data.base$northing,col=data.base$FR_M-1,cex=0.1,asp=1)
sel<-which(data.base$pos==1)
# lines(data.base$easting[sel],data.base$northing[sel],col="red",type="p",pch=16,cex=0.3)
legend("topright",title="Month of spray",levels(as.factor(data.base$FR_M)),col=levels(as.factor(data.base$FR_M-1)),pch=16)

## localities
# get the localities (take a few seconds)
loc<-as.integer(t(as.vector(as.data.frame(lapply(as.character(data.base$unicode),strsplit,".",fixed=TRUE))[3,])))

dev.new()
plot(data.base$easting,data.base$northing,col=loc,cex=0.1,asp=1)
sel<-which(data.base$pos==1)
# lines(data.base$easting[sel],data.base$northing[sel],col="red",type="p",pch=16,cex=0.3)
# legend("topright",title="Year of spray",levels(as.factor(data.base$FR_A)),col=levels(as.factor(data.base$FR_A)),pch=16)

## subset selection (fall 2007)
sel<-which(loc==46 | loc == 48 | loc == 50 | loc == 51 | loc == 52 | loc == 55)
plot(data.base$easting,data.base$northing,cex=0.1,asp=1)
lines(data.base$easting[sel],data.base$northing[sel],col="green",type="p",pch=16,cex=0.3)
## only one year:
# > levels(as.factor(data.base$FR_A[sel]))
# [1] "7"
## only 4 months of fall
# > levels(as.factor(data.base$FR_M[sel]))
# [1] "9"  "10" "11" "12"

# similar but for fall 2008
dev.new()
sel<-which(loc==87 | loc==88 | loc == 90 | loc == 92 | loc == 94)
plot(data.base$easting,data.base$northing,cex=0.1,asp=1)
lines(data.base$easting[sel],data.base$northing[sel],col="green",type="p",pch=16,cex=0.3)

## controlling for only one year
# > levels(as.factor(data.base$FR <- A[sel]))
# [1] "8"
## controlling for months
# > levels(as.factor(data.base$FR <- M[sel]))
# [1] "10" "11" "12"



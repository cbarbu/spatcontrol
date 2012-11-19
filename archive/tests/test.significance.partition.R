### parameters
graphics.off()
lim.dist.classes = seq(20,40,20)

### data
data.base<-read.csv("DB_simple_Pau_cyclo1_19Jul2011.csv")
### functions
library("testthat")
library("spam")
source("spam_complement.r")
source("functions_intercept.r")
source("R/make.partition.matrices.R")
source("R/agreement.with.neigh.R")

### treatment
notNA<-which(data.base$status!=9)
x<-data.base$easting[notNA]
y<-data.base$northing[notNA]
z<-data.base$status[notNA]
p<-data.base$block_num[notNA]

source("R/test.significance.partition.R")
# TSi<-test.significance.partition(x,y,z,p,lim.dist.classes,type="")
# TS1<-test.significance.partition(x,y,z,p,lim.dist.classes,type="1")
# TS0<-test.significance.partition(x,y,z,p,lim.dist.classes,type="0")

TSa<-test.significance.partition(x,y,z,p,lim.dist.classes)
print(TSa)
TS<-TSa

par(mfrow=c(2,3))
plot(TS$agreement.rate,ylim=c(min(TS$agreement.rate.DGr,TS$agreement.rate.SGr),max(TS$agreement.rate.SGr,TS$agreement.rate.DGr)),type="b")
lines(TS$agreement.rate.DGr)
lines(TS$agreement.rate.SGr)
plot(TS$mean.signif.p.values)
plot(TS$min.signif.p.values)
plot(log(TS$global.p.values))

dev.new(width=6,height=3)
par(mfrow=c(1,2),mar=c(4,4,0.5,0.5))
plot(TS$mean.signif.p.values,ylab="mean(p.value_i)",xlab="distance class")
plot(TS$min.signif.p.values,ylab="min(p.value_i)",xlab="distance class")
dev.print(device=pdf,"initialPvalue.pdf")

dev.new(width=3,height=3)
par(mfrow=c(1,1),mar=c(4,4,0.5,0.5))
plot(TS$agreement.rate,ylim=c(min(TS$agreement.rate.DGr),max(TS$agreement.rate.SGr)),type="b",xlab="distance class",ylab="agreement ratio")
lines(TS$agreement.rate.DGr,type="b",lty=2,pch=2)
lines(TS$agreement.rate.SGr,type="b",lty=3,pch=3)
legend("topright",c("rH","rD","rS"),lty=c(1,2,3),pch=c(1,2,3))
dev.print(device=pdf,"visuRatios.pdf")

dev.new(width=3,height=3)
par(mfrow=c(1,1),mar=c(4,4,0.5,0.5))
plot(log10(TS$global.p.values),ylab="Global p-value (log10)",xlab="Distance classes")
dev.print(device=pdf,"globalPvalues.pdf")



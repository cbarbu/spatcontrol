par(mfrow=c(2,2))
y<-sample_y_direct(w,zpos,zneg,zNA,bivect);
yprime <- (y>0);
yprimeA<-sampleyprime(dimension,w,bivect,zpos,zneg,zNA);
yA<-sampley(dimension,w,yprimeA);
hist(as.numeric(yprime))
hist(y,breaks=17)
hist(as.numeric(yprimeA))
hist(yA,breaks=17)
 
dev.new()
par(mfrow=c(2,2))
plot_reel(data$easting,data$northing,yA)
plot_reel(data$easting,data$northing,yprimeA)
plot_reel(data$easting,data$northing,y)
plot_reel(data$easting,data$northing,yprime)

# while the yprime are very similars the y are quite different
# as it may be a randomness irregularity let us check
# on n tirages
graphics.off()
ybreaks<-seq(-10,10,0.5)
dev.new()
nrep=250
par(mfrow=c(3,5))
sumHistyA<-rep(0,length(ybreaks)-1)
for(i in 1:nrep){
	yprimeA<-sampleyprime(dimension,w,bivect,zpos,zneg,zNA);
	yA<-sampley(dimension,w,yprimeA);
	histyA<-hist(yA,breaks=ybreaks)
	sumHistyA<-sumHistyA+histyA$counts
	cat(" ",sum(histyA$counts))
}
sumHisty<-rep(0,length(ybreaks)-1)
dev.new()
par(mfrow=c(3,5))
for(i in 1:nrep){
	y<-sample_y_direct(w,zpos,zneg,zNA,bivect);
	histy<-hist(y,breaks=ybreaks)
	sumHisty<-sumHisty+histy$counts
	cat(" ",sum(histyA$counts))
}

dev.new()
par(mfrow=c(1,2))
histy$counts<-round(sumHisty/nrep)
histy$counts<-round(sumHistyA/nrep)

plot(histy,ylim=c(0,200))
plot(histyA,ylim=c(0,200))



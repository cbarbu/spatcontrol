setwd("/Users/mzlevy/Documents/Epi_in_R/Corrected_prevalence_Melgar")
#install.packages("msm")
install.packages("LaplacesDemon")


INTERMEDIARY<-TRUE
library(maptools)
source("parameters_sampler.r")
source("functions_intercept.r")
source("pseudo_data_generation.r")
source("prep_sampling.r")
# system("tail -n 1 usamples.txt > usampleslast.txt")
# system("tail -n 1 wsamples.txt > wsampleslast.txt")
sampled<-read.table(file="sampled.txt")
nbsimul<-length(sampled[,1])

#=================================
#	-Old using the 'w' unstead of y''
#=================================
# u<-drop(t(scan(file="usampleslast.txt",sep="\t")))
w<-drop(t(scan(file="wsampleslast.txt",sep="\t")))
y<-drop(rnorm(length(w),w,1));
z<-generate_z(y,bivect);
zpos <- which(z.r==1);
zneg <- which(z.r==0);
zNA <- which(z.r==9);

#=================================
#New: get the results of the simulated y''s
#	-without streets
#	-with inspectors
#=================================
yprime <- drop(t(scan(file = "ypsamplesLast.txt", sep = "\t")))
yprime<- matrix(yprime,nrow=1000,byrow=TRUE) 

yp <- apply(yprime, 2, mean)

#=================================
#	Locality-level Prevalence estimates
#	1-model-based
#	2-raw data, Prevalence = infested/ total open or known
#=================================

#counting up the mean infestation from model
infest_status <- min(db$locality):max(db$locality)
numhouses <- infest_status

for(i in min(db$locality):max(db$locality))
{
	infest_status[i] <- mean(yp[which(db$locality == i)])
	numhouses[i] <- length(which(db$locality == i))
}

infest_status[which(is.nan(infest_status))] <- NA
# stat_vals <- (max(infest_status) - infest_status)/max(infest_status)
# stat_vals <- 1-stat_vals

#counting up the mean infestation from RAW DATA

infest_known <-open<-infested<- min(db$locality):max(db$locality)
for(i in min(db$locality):max(db$locality))
{
	infested[i] <- length(which(db[which(db$locality == i), "status"] == 1))
	notOpened <- length(which(db[which(db$locality == i), "status"] == 9))
	open[i]<-length(which(db$locality == i)) - notOpened
	percent <- infested[i]/open[i]
	
	infest_known[i] <- percent
}

infest_known[which(is.nan(infest_known))] <-  NA
# known_vals <- (max(infest_known) - infest_known)/max(infest_known)
# known_vals <- 1-known_vals

#=================================
#	Mapping Locality-level Prevalence estimates
#	1-model-based
#	2-raw data
#=================================

# 1. plots predicted infestation status versus measured infestation status by the ministry
dev.new()
par(mfrow = c(1, 2))
plot_reel(db$easting, db$northing, infest_status[db$locality], base=0,top=0.4)

#raw data
plot_reel(db$easting, db$northing, infest_known[db$locality], base=0,top=.4)


#=================================
#	Store prevalence estimates by locality
#	1-predict (model-based)
#	2-vigil (raw data)
#=================================
predict <- data.frame(locality = which(!is.na(infest_status)), p = infest_status[which(!is.na(infest_status))])
vigil <-data.frame(locality = which(!is.na(infest_known)), p = infest_known[which(!is.na(infest_known))])

#=================================
#	plot prevalence estimates by locality
#=================================
#  locality prevalence
dev.new()
par(mfrow = c(2, 1))

plot(predict$locality, predict$p, ylim = c(0, 1))
abline(h=mean(predict$p), col="red")
plot(vigil$locality, vigil$p, ylim = c(0, 1))
abline(h=mean(vigil$p), col="red")


#=================================
#	cross-plots comparing model-based and raw data prevalence estimates
#=================================

plot(vigil$p,predict$p, xlim=c(0,.7),ylim=c(0,.7))
lines(0:100/100,0:100/100)


#=================================
#	Table of model estimates and raw data
#=================================
sel=which(!is.na(infest_known))

TABLE<-data.frame(diff$locality, infested[sel],open[sel], vigil$p, predict$p, predict$p-vigil$p)
names(TABLE)<-c("Locality","Infested","Total Observed","Observed Prevalence", "Model Prevalence Estimate", "Estimated Difference")
write.csv(TABLE, "Estimated_Prevalence_Localities_Melgar.csv", row.names=False)





















#######################old below!

#=================================
#	plot difference in model - raw data prevalence estimates by locality
#=================================
# examine the differences between the two maps by locality
diff <- data.frame(locality = predict$locality, pdiff = predict$p - vigil$p)
diff <- cbind(diff, percentdiff = diff$pdiff/vigil$p)

dev.new()
par(mfrow = c(2, 1))
plot(diff$locality, diff$pdiff)
abline(h=mean(diff$pdiff), col = "red")
plot(diff$locality, diff$percentdiff)
abline(h=mean(diff$percentdiff), col="red")

dev.new()
par(mfrow = c(2, 1))
plot(predict$p[which(vigil$p == 0)]~numhouses[predict$locality[which(vigil$p == 0)]])
plot(numhouses)

dev.new()
par(mfrow = c(2, 1))
plot(numhouses[diff$locality], diff$percentdiff)
abline(h=mean(diff$percentdiff[which(vigil$p!=0)]), col="green")
abline(h=median(diff$percentdiff[which(vigil$p!=0)]), col="blue")
abline(v=mean(numhouses[vigil$locality[which(vigil$p!=0)]]), col = "green")
abline(v=median(numhouses[vigil$locality[which(vigil$p!=0)]]), col = "blue")

plot(numhouses[diff$locality], diff$pdiff)
abline(h=mean(diff$pdiff), col = "green")
abline(h=median(diff$pdiff), col = "blue")
abline(v=mean(numhouses), col = "green")
abline(v=median(numhouses), col = "blue")

dev.new()
par(mfrow = c(1, 2))
plot_reel(db$easting, db$northing, diff$pdiff[match(db$locality, diff$locality)], base=min(diff$pdiff), top=max(diff$pdiff))

want <- which(db$locality %in% vigil$locality[which(vigil$p != 0)])
plot_reel(db$easting[want], db$northing[want], diff$percentdiff[db$locality[want]], base = min(diff$percentdiff[which(vigil$p!=0)]), top = max(diff$percentdiff[which(vigil$p!=0)]))

dev.new()
noninf <- which(vigil$p!=0)
plot(diff$percentdiff[noninf]~diff$pdiff[noninf])
abline(h=mean(diff$percentdiff[which(vigil$p!=0)]), col="green")
abline(v=mean(diff$pdiff[which(vigil$p!=0)]), col = "green")
sto <- lm(diff$percentdiff[noninf]~diff$pdiff[noninf])
show(summary(sto))
abline(sto, col= "pink")


dev.new()
par(mfrow = c(1, 2))
plot(numhouses[diff$locality[order(diff$pdiff, decreasing = TRUE)]])
plot(numhouses[diff$locality[rev(order(diff$percentdiff)[1:length(which(vigil$p!=0))])]])

dev.new()
par(mfrow = c(1, 3))
plot(numhouses[which(numhouses!=0)]~diff$pdiff)
plot(numhouses[which(numhouses!=0)][noninf]~diff$pdiff[noninf])
plot(numhouses[which(numhouses!=0)]~diff$percentdiff)

##########
# ordering by locality

vl <- vigil$locality[order(vigil$p, decreasing = TRUE)]
pl <- predict$locality[order(predict$p, decreasing = TRUE)]

rankings <- match(vl, pl)
for(i in 1:length(rankings))
	if(rankings[i]!=i)
	{
		show(i)
		show(vl[i])
		show(vigil$p[i])
		show(pl[i])
		show(predict$p[i])
		print("\n")
	}

rankings <- data.frame(locality = vigil$locality, vplace = match(vigil$locality, vl), pplace = match(vigil$locality, pl))
rankings <- cbind(rankings, diff = rankings$vplace - rankings$pplace)
averagerank <- rankings$pplace + rankings$vplace
averagerank <- averagerank/2
weight <- length(averagerank) + 1 - averagerank
normweight <- 1 - ((max(weight) - weight)/max(weight))
rankings <- cbind(rankings, weighteddiff = abs(rankings$diff*normweight))
print(rankings)

implocs <- rankings$locality[order(rankings$weighteddiff, decreasing = TRUE)]
imppdiff <- diff$locality[order(diff$pdiff, decreasing = TRUE)]
impabsdiff <- rankings$locality[order(abs(rankings$diff), decreasing = TRUE)]
impperdiff <- diff$locality[order(diff$percentdiff, decreasing = TRUE)[(length(diff$percentdiff)+1-length(noninf)) : length(diff$percentdiff)]]
print(implocs)
print(impabsdiff)
print(imppdiff)
print(impperdiff)
print(vl)
print(pl)
vals1 <- sort(rankings$weighteddiff, decreasing = TRUE)
vals2 <- sort(abs(rankings$diff), decreasing = TRUE)
vals3 <- sort(diff$pdiff, decreasing = TRUE)
vals4 <- sort(diff$percentdiff, decreasing = TRUE)[(length(diff$percentdiff)+1-length(noninf)) : length(diff$percentdiff)]

dev.new()
par(mfrow = c(2, 4))

plot_reel(db$easting, db$northing, rankings$weighteddiff[match(db$locality, rankings$locality)], base=min(rankings$weighteddiff), top=max(rankings$weighteddiff))

plot_reel(db$easting, db$northing, abs(rankings$diff[match(db$locality, rankings$locality)]), base=min(abs(rankings$diff)), top=max(abs(rankings$diff)))

plot_reel(db$easting, db$northing, diff$pdiff[match(db$locality, diff$locality)], base=min(diff$pdiff), top=max(diff$pdiff))

want <- which(db$locality %in% vigil$locality[which(vigil$p != 0)])
plot_reel(db$easting[want], db$northing[want], diff$percentdiff[db$locality[want]], base = min(diff$percentdiff[which(vigil$p!=0)]), top = max(diff$percentdiff[which(vigil$p!=0)]))

plot(implocs, vals1)
plot(impabsdiff, vals2)
plot(imppdiff, vals3)
plot(impperdiff, vals4)

par(mfrow = c(1, 4))
plot_reel(db$easting, db$northing, infest_status[db$locality], base=0)
title("Initial Data")
plot_reel(db$easting, db$northing, infest_known[db$locality], base=0)
title("Model Predctions")
plot_reel(db$easting, db$northing, diff$pdiff[match(db$locality, diff$locality)], base=0)
title("Probability Difference")
plot(vigil$p, predict$p)
title("x = initial, y= model")

stop()

# # hist of nb houses per block
# db$count<-1
# count<-aggregate(db$count,by=list(db$block_num),sum)
# quantile(count$x,prob=c(0.05,0.25,0.5,0.75,0.95))
localities<-getKMLcoordinates("Poligonos__MM.kml") # don't use ignoreAltitude as it cause for some reason only the first two points of the polygon to be imported
get_list_poly<-function(listfromKML){
	truc <- {}
loc=list();
	for( i in 1:length(listfromKML)){
		loc[i]<-listfromKML[[i]][1]
		for( j in 2:(length(listfromKML[[i]])-1)){
			loc[[i]]<-rbind(loc[[i]],listfromKML[[i]][[j]])
		}
		# lines(loc[[i]])
		# locs[i]<-as(loc[[i]][,1:2],"gpc.poly")
		truc<-c(truc,Polygon(loc[[i]][,1:2]))
	}
	list_poly<-SpatialPolygons(list(Polygons(truc,c("1"))))
	# plot(list_poly)
	return(list(loc,list_poly))
}
get_spat_poly<-function(listfromKML){
	truc <- {}
	for( i in 1:length(listfromKML)){
		truc<-c(truc,Polygon(listfromKML[[i]][,1:2]))
	}
	list_poly<-SpatialPolygons(list(Polygons(truc,c("1"))))
	# plot(list_poly)
	return(list_poly)
}
loc_lists<-get_list_poly(localities)
loc_list<-loc_lists[[1]]
coord<-as.data.frame(matrix(c(-71.475,-16.415,-71.52,-16.385),ncol=2,byrow=TRUE))

# my way to plot
plot(coord,asp=1,type="n")
for(i in 1:length(loc_list)){
	polygon(loc_list[[i]],col=i)
}
	
# native way but doesn't color well
loc_for_plot<-loc_lists[[2]]
plot(loc_for_plot)
#


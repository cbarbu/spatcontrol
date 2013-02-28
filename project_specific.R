unicodeToPDLV<-function(unicodes){
	PDLV<-as.data.frame(t(simplify2array(strsplit(as.character(unicodes),"\\."))))
	names(PDLV)<-c("P","D","L","V")
	return(PDLV)
}
# remove the trailing a-Z in unicodes to be able to merge with arequipa_gps_google 
make_unicode_gps<-function(unicodes){
	return(gsub("[A-Z]$","",unicodes,ignore.case=TRUE))
}
# Clean P,D,L,V from spaces
cleanPDLV<-function(Table){
	Table$P<-gsub(" ","",Table$P)
	Table$D<-gsub(" ","",Table$D)
	Table$L<-gsub(" ","",Table$L)
	Table$V<-gsub(" ","",Table$V)
	Table$unicode<-paste(Table$P,Table$D,Table$L,Table$V,sep=".")
	return(Table)
}
# clean unicodes from spaces and lower case
cleanUnicodes<-function(unicodes){
	unicodes<-gsub(" ","",unicodes)
	unicodes<-toupper(unicodes)
	return(unicodes)
}

########################
# plot shapefiles of locality limits
########################
## regenerate the R object file:
# mapLim<-readShapeSpatial("LocalitiesBoundaries/Arequipa_limits")
# mapLim@data
# mapLim@polygons

# # remove the points at the end of the Code if needed
# names(mapLim@data)<-c("Name","Code")
# mapLim@data$Code<-gsub("\\.$","",mapLim@data$Code)
# 
# # split provincia, district and localily level codes
# mapLim@data$P<-gsub("^[0-9]*\\.","",mapLim@data$Code)
# # set the projection system
# proj4string(mapLim)<-CRS("+proj=longlat +datum=WGS84")

# save it as a separate file
# save(mapLim,file="ArequipaLim.img")

## can then be load using simply
load("ArequipaLim.img")

### function to plot a given set of localities
# locNames should be a vector of character of the form 1.9.4
# * and other characters for regex can be used
# additional standar arguments to plot can be passed 
# in particular add=TRUE to add a locality on the top of a plot
plot.loc.arequipa<-function(locNames,utm=TRUE,...){
	# get the number of the polygons
	listNum<-c()
	for(name in locNames){
		listNum<-c(listNum,grep(name,mapLim@data$Code))
	}
	if(utm){
		library("rgdal")
		mapLim<-spTransform(mapLim, CRS("+proj=utm +zone=19 +south"))
	}

	plot(mapLim[listNum,],...)
	return(invisible(listNum))
}


## then can be used directly
# # example
# par(mfrow=c(1,2))
# plot(mapLim) # to plot all the localities
# locNames<-c("1.8.12","1.8.30")
# plot.loc.arequipa(locNames,add=TRUE,col="blue")
# 
# plot.loc.arequipa("1.8.*")
# plot.loc.arequipa(locNames,add=TRUE,col="blue")
# # and to plot something (like text/data) at the centroid of the polygons you can use
# text(mapLim@polygons[[1]]@labpt[1],mapLim@polygons[[1]]@labpt[2],mapLim@data$Code[[1]])



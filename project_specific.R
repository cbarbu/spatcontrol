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
	unicodes<-toupper(bugs$UNICODE)
	return(unicodes)
}



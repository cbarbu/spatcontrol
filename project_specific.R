unicodeToPDLV<-function(unicodes){
	PDLV<-as.data.frame(t(simplify2array(strsplit(as.character(unicodes),"\\."))))
	names(PDLV)<-c("P","D","L","V")
	return(PDLV)
}
# remove the trailing a-Z in unicodes to be able to merge with arequipa_gps_google 
make_unicode_gps<-function(unicodes){
	return(gsub("[A-Z]$","",unicodes,ignore.case=TRUE))
}



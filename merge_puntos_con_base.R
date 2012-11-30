source("spatcontrol/spatcontrol.R",local=TRUE,chdir=TRUE)

# 1 - .csv.R generado por check_points_to_poly.sh sobre los kmz
archivoPuntos<-"Miraflores_point-blocks.csv.R"
# 2 - archivo con datos de encuestas, rociado etc...
archivoDatos<-"miraflores_eed.csv"
# 3 - nombre final
archivoFinal<-"encuesta_miraflores.csv"

# importacion
points<-read.csv(archivoPuntos,header=TRUE)
obs<-read.csv(archivoDatos,header=TRUE)

# limpieza / preparation
obs<-set_to(obs,init=c("NULL"),final=0)

obs$unicode<-paste(obs$P,obs$D,obs$L,obs$V,sep=".")

# Join
enc<-merge(points,obs,by.x="unicode",by.y="unicode",all.x=FALSE,all.y=FALSE)

# guardar
write.csv(enc,archivoFinal,row.names=FALSE)


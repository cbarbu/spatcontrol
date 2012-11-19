nameSimul<-"extrapol_rocI_exp_epsKc0.01NoStreetsTr50_noinsp_nocof_nofitspat_f8_Kok_fast_final"
set.seed(1)
# # generate the known
# unicode<-paste("1.13.88",c(481:482,484:485,553,555:558),sep=".")
# opened<-rbinom(length(unicode),1,prob=0.7)
# infestation<-rbinom(sum(opened),1,prob=0.5)
# infested<-rep(0,length(unicode))
# infested[which(opened==1)]<-infestation
# 
# FR_D<-rep(8,length(unicode))
# FR_M<-rep(7,length(unicode))
# FR_A<-rep(2011,length(unicode))
# known<-data.frame(cbind(unicode,opened,infested,FR_D,FR_M,FR_A))
# # edit(known)

# write.csv(known,file="knownSubset.csv",row.names=FALSE)

# execute the kernel
source("extrapol_field.r")
# library("realtimeadapt")
subsetAround <- function (priordatafullmap, reportsUnicodes, threshold, ...) {
	reportslines <- match(reportsUnicodes, priordatafullmap$unicode)
	if(sum(is.na(reportslines))>0){
		stop("new points not in initial map") # in fine better to send a warning
	}
	dist_mat <- nearest.dist(x = priordatafullmap[, c("X", "Y")], 
		y = priordatafullmap[reportslines, c("X", "Y")], method = "euclidian", 
		delta = threshold, upper = NULL)
	selectedlines <- which(!is.na(apply_by_row_not_null.spam(dist_mat, 
				min)))
	subprior <- priordatafullmap[selectedlines, ]
	return(subprior)
}

## data
data.base <- read.csv("DB_simple_Pau_cyclo1_19Jul2011.csv",header=TRUE);

data.base$oanimal<-as.integer((data.base$CO==1 | data.base$AV==1 | data.base$GA==1 | data.base$OV==1 | data.base$otros.animales!="-1"))

data.base$TrueStatus<-data.base$status
data.base$status[data.base$status==9]<-9
data.base<-data.base[-which(duplicated(cbind(data.base$easting,data.base$northing))),]

## new map
newdata<-read.csv("knownSubset.csv")
newknowns<-newdata[which(newdata$opened == 1),]
data.base$X<-data.base$easting
data.base$Y<-data.base$northing
torefit<-subsetAround(data.base,newknowns$unicode,100)
new.status<-rep(9,length(torefit$status))
newknowns$status<-rep(0,length(newknowns$unicode))
newknowns$status[newknowns$opened==0]<-9
newknowns$status[newknowns$infested==1]<-1
new.status[match(newknowns$unicode,torefit$unicode)]<-newknowns$status
torefit$new.status<-new.status
torefit$old.status<-torefit$status
torefit$status<-new.status
torefit$opened<-as.numeric(torefit$status!=9)
torefit$infested<-as.numeric(torefit$status==1)
torefit$p.i<-rep(0.07959277,length(torefit$status))

##### 
updated<-extrapol.spatautocorel(db=torefit[,c("X","Y","opened","infested","collector","block_num","CU","PE","oanimal","I.NO","P.NO","p.i")],pfile="parameters_extrapol.r")
p.i<-updated$p.i
#####
plot_reel(torefit$X,torefit$Y,p.i,base=0,top=1)

# ## post spray rules:
# p.i<-est.yp # base is the estimate by the model
# # post spray corrections
# pstayinfested<-0.00705 # probability of infested if sprayed infested
# p.i[torefit$status==1]<-pstayinfested
# # probability of newly infested 
# pnewinfestation6months<-0.001101 
# # probability of infested if not observed infested
# p.i[torefit$status==0]<-p.i[torefit$status==0]*pstayinfested+pnewinfestation6months
# 
# # probability of infested if not sprayed
# p.i[torefit$status!=1]<-p.i[torefit$status!=1]+pnewinfestation6months

X<-torefit$easting
Y<-torefit$northing
unicode<-as.character(torefit$unicode_gps)
opened<-as.numeric(torefit$status!=9)
infested<-as.numeric(torefit$status==1)
block_num<-torefit$block_num
FR_D<-torefit$FR_D
FR_A<-torefit$FR_A
FR_M<-torefit$FR_M

## get 
probamap<-data.frame(X,Y,unicode,block_num,FR_D,FR_M,FR_A,opened,infested,p.i)
dev.new()
plot_reel(probamap$X,probamap$Y,probamap$p.i,base=0,top=1)
sel<-which(probamap$opened==1)
lines(probamap$X[sel],probamap$Y[sel],type="p",pch=1)
sel<-which(probamap$infested==1)
lines(probamap$X[sel],probamap$Y[sel],type="p",pch=10)
dev.new()
plot(probamap$p.i)

dev.new()
out<-grid.from.kernel(probamap$X,probamap$Y,probamap$p.i,steps=100,f=8.5,tr=100)
dist.weight<-matrix(as.vector(as.matrix(out$z)),nrow=length(out$xs))
image(x=out$xs,y=out$ys,z=dist.weight,asp=1,col=heat.colors(100),xlab="x (m)",ylab="y (m)")
lines(probamap$X,probamap$Y,pch=".",type="p")

sel<-which(probamap$opened==1)
lines(probamap$X[sel],probamap$Y[sel],type="p",pch=1)
sel<-which(probamap$infested==1)
lines(probamap$X[sel],probamap$Y[sel],type="p",pch=10)



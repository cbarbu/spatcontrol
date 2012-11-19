library("spam")
library("testthat")
dyn.load("useC2.so")
source("R/symmetric.random.entries.by.line.spam.R")
source("R/symmetric.random.entries.by.line.spam.c.R")

graphics.off()
A<-matrix(c(0,0,0,1,0,0,
	    0,0,1,0,0,0,
	    0,1,0,0,0,0,
	    1,0,0,0,0,0,
	    0,0,0,0,0,1,
	    0,0,0,0,1,0),nrow=6)

N<-matrix(c(0,1,0,1,0,0,
	    1,0,1,0,1,0,
	    0,1,0,0,0,1,
	    1,0,0,0,1,0,
	    0,1,0,1,0,1,
	    0,0,1,0,1,0),nrow=6)
coord<-matrix(c(0,0,
		1,0,
		2,0,
		0,1,
		1,1,
		2,1),nrow=6,byrow=TRUE);
draw.sep<-function(xa,ya,xb,yb,...){
	xmid<-(xa+xb)/2
	ymid<-(ya+yb)/2
	x1<-xmid-abs(ymid-ya)/2
	x2<-xmid+abs(ymid-ya)/2
	y2<-ymid-abs(xmid-xa)/2
	y1<-ymid+abs(xmid-xa)/2
	# cat(x1,y1,x2,y2,"\n")
	lines(c(x1,x2),c(y1,y2),...)
}
plot.sep<-function(coord,A,add=FALSE,...){
	# inverse entries in A to ease the plot
	pos<-which(A@entries==1)
	neg<-which(A@entries==0)
	A@entries[pos]<-0
	A@entries[neg]<-1

	# actual plot
	if(!add){
		plot(coord[,1],coord[,2],asp=1)
	}
	for(i in 1:dim(A)[1]){
		for(j in i:dim(A)[1]){
			if(A[i,j]==1){
				draw.sep(coord[i,1],coord[i,2],coord[j,1],coord[j,2],...)
			}
		}
	}
}
plot.link<-function(coord,A,add=FALSE,...){
	# A matrix with 1 if linked
	# actual plot
	if(!add){
		plot(coord[,1],coord[,2],asp=1)
	}
	for(i in 1:(dim(A)[1]-1)){
		to.plot<-which(A[i,]==1)
		to.plot<-to.plot[to.plot>i]
		for(j in to.plot){
			lines(c(coord[i,1],coord[j,1]),c(coord[i,2],coord[j,2]),...)
		}
	}
}

N<-as.spam(N)
phantom<-N
phantom@entries<-rep(0,length(phantom@entries))
A<-phantom+as.spam(A)
Arand<-symmetric.random.entries.by.line.spam.c(A)
# Arand<-symmetric.random.entries.by.line.spam.R(A)

# # test draw.sep
# plot(coord[,1],coord[,2],asp=1)
# draw.sep(0,0,0,1)
# draw.sep(0,0,1,0)
# draw.sep(0,0,1,1)
# draw.sep(1,0,2,1)
# check total number of 0/1 is conserved
expect_equal(length(Arand@entries),length(A@entries))
expect_equal(length(which(Arand@entries==0)),length(which(A@entries==0)))
expect_equal(length(which(Arand@entries==1)),length(which(A@entries==1)))

# check symmetry of A
expect_true(isSymmetric(Arand))

par(mfrow=c(1,2))
plot.link(coord,A);
plot.link(coord,Arand);

## on a true map
use.map.gen<-FALSE
city="PAUCARPATA"
period=""
subsetcity=0
use.NA<-FALSE
threshold=50
# source("import_data.r")

spam.options(nearestdistnnz=c(9058076,400))
dist_mat <-nearest.dist(x=data[,c("easting","northing")], y=NULL, method="euclidian", delta=threshold, upper=NULL);          
dimension <- nrow(data);
spam.options(nearestdistnnz=c(13764100,400))
SB <- nearest.dist(x=cbind(data$block_num,rep(0,length(data$block_num))), method="euclidian", upper=NULL,delta=0.1)
SB@entries<-rep(1,length(SB@entries))
dmt<-dist_mat
dmt@entries<-rep(1,length(dmt@entries))# [dmt@entries!=0]<-1 # 1 only when dist_mat not 0
SB@entries<-rep(1,length(SB@entries))
SB<-SB*dmt;

AS<-dmt-SB; # get 1 whereever the distances matrix is defined(under threshold) and not same block
AS<-as.spam(AS)

limiteinf<-20
limitesup<-40
dmtr<-dmt # will be the matrix of neighbours for a ring(1/0)
unselect<-which(dist_mat@entries>=limitesup | dist_mat@entries<limiteinf)
dmtr@entries[unselect]<-rep(0,length(unselect))# 1 when dist_mat defined (distance<max(distance))
dmtr<-as.spam(dmtr)
cat(" mean nbNeigh",sum(dmtr/2)/dimension)

ASr<-AS*dmtr;
SBr<-SB*dmtr;
# cat("pre calculus loaded\n")

N<-as.spam(dmtr)
phantom<-N
phantom@entries<-rep(0,length(phantom@entries))
A<-phantom+as.spam(SBr)
Arand<-symmetric.random.entries.by.line.spam.c(A)
# Arand<-symmetric.random.entries.by.line.spam.R(A)

# check total number of 0/1 is conserved
expect_equal(length(Arand@entries),length(A@entries))
expect_true(length(which(Arand@entries!=A@entries))>0)
expect_equal(length(which(Arand@entries==0)),length(which(A@entries==0)))
expect_equal(length(which(Arand@entries==1)),length(which(A@entries==1)))

# check symmetry of A
expect_true(isSymmetric(Arand))

# par(mfrow=c(1,2))
# 
# # plot.sep(cbind(data$easting,data$northing),A)
# # plot.sep(cbind(data$easting,data$northing),Arand)
# 
# plot.link(cbind(data$easting,data$northing),A)
# plot.link(cbind(data$easting,data$northing),Arand)

# #plot link across streets
# plot.link(cbind(data$easting,data$northing),N-A)
# plot.link(cbind(data$easting,data$northing),N-Arand)
# 
# plot.sep(cbind(data$easting,data$northing),N-A)
# plot.sep(cbind(data$easting,data$northing),N-Arand)
# 
# 
# plot.link(cbind(data$easting,data$northing),N)
# plot.link(cbind(data$easting,data$northing),Arand,add=TRUE,col=4)
# plot.sep(cbind(data$easting,data$northing),N)
# plot.sep(cbind(data$easting,data$northing),Arand,add=TRUE,col=4)


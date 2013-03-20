# ===================
# gpclib approach: pb, does whatever triangulation
# and add points to the limit of the polygons
# ===================

library(gpclib)
gpctri2tri<-function(gpctri,coords){
	colnames(coords)<-c("x","y")
	nvert<-dim(gpctri)[1]
	gpctriRef<-rep(0,nvert)
	for(i in 1:nvert){
		gpctriRef[i]<-identifyFromPoints(gpctri[i,1],gpctri[i,2],coords)
	}
	tri<-matrix(gpctriRef,ncol=3,byrow=TRUE)
	# tri<-list()
	# for(i in 0:(nrow(gpctri)/3 - 1)){
	# 	tri[[i+1]]<-gpctriRef[3*i+1:3]
	# }
	return(tri)
}
# test
x <-structure(c(0.0934073560027759, 0.192713393476752, 0.410062456627342,0.470020818875781, 0.41380985426787, 0.271408743927828, 0.100902151283831,
0.0465648854961832, 0.63981588032221, 0.772382048331416,0.753739930955121, 0.637744533947066, 0.455466052934407,0.335327963176065, 0.399539700805524,0.600460299194476), .Dim = c(8, 2))
y <- structure(c(0.404441360166551, 0.338861901457321, 0.301387925052047, 
		    0.404441360166551, 0.531852879944483, 0.60117973629424, 0.625537820957668, 
		    0.179976985040276, 0.341542002301496, 0.445109321058688,
		    0.610817031070196, 0.596317606444189, 0.459608745684695,
		    0.215189873417722), .Dim = c(7, 2))

x1 <- as(x, "gpc.poly")
y1 <- as(y, "gpc.poly")
## Show the triangulation
par(mfrow=c(1,2))
plot(append.poly(x1, y1))
triangles <- triangulate(append.poly(x1,y1))
for (i in 0:(nrow(triangles)/3 - 1)) 
	polygon(triangles[3*i + 1:3,], col="lightblue")

# transform to our format
coordsxy<-rbind(x,y)
colnames(coordsxy)<-c("x","y")
tri<-gpctri2tri(triangles,coordsxy)
plotMesh(vertices=coordsxy,triangles=tri)

# =======================
# triangulation based on constrained delaunay, with houses
# use tripack
# =======================
library(tripack)
library(splancs)

# from a constrained delaunay using tripack
# remove what is out of boundary region
removeOutSegments<-function(tpk){
	permut<-matrix(c(1,2,2,3,3,1),ncol=2,byrow=TRUE)
	# the constraint correspond to the external polygon
	beginConst<-min(tpk$lc)
	constx<-tpk$x[beginConst:tpk$n]
	consty<-tpk$y[beginConst:tpk$n]
	
	triangles<-triangles(tpk)[,1:3]
	nbtri<-dim(triangles)[1]
	# should limit to triangles with an two external points
	okTri<-rep(TRUE,nbtri)
	for(numtri in 1:nbtri){
		triangle<-triangles[numtri,]
		# triangle<-triangles[[numtri]]
		ptsMid<-mat.or.vec(3,2)
		triCurrentOK<-rep(FALSE,3)
		for(side in 1:3){
			refs<-triangle[permut[side,]]
			adRef<-abs(diff(refs))

			refsAreNext<-(adRef!=1 && adRef!=(length(tpk$x)-beginConst))  
			if(refsAreNext){ # don't test consecutive 
				# as consecutive points on border are ok
				pts1<-c(tpk$x[refs[1]],tpk$y[refs[1]])
				pts2<-c(tpk$x[refs[2]],tpk$y[refs[2]])
				ptsMid[side,]<-midPoint(pts1,pts2)
			}else{
				triCurrentOK[side]<-TRUE
			}
		}
		okTri[numtri]<-all(point.in.polygon(ptsMid[,1],ptsMid[,2],constx,consty,mode.checked=TRUE)!=0 | triCurrentOK)
	}
	pts<-cbind(tpk$x,tpk$y)
	colnames(pts)<-c("x","y")
	return(list(pts=pts,tri=triangles[which(okTri),],bc=beginConst))
}
# # Test
# xs<-c(232069.7,232094.1,232120.8,232126.8,232155.9,232192.0,232161.2,232119.0,232027.8,232052.4,232067.8,232081.9)
# ys<-c(8184673,8184688,8184669,8184667,8184642,8184601,8184585,8184556,8184644,8184663,8184641,8184653)

# # use the houses as interior points
# # could use midPoints of all delimitating points that are in poly 
# load("byHouseVig.img")
# a<-(point.in.polygon(byHouse$X,byHouse$Y,xs,ys)==1)
# pts<-byHouse[which(a),c("X","Y")]
# pts<-pts[!duplicated(pts),]
# names(pts)<-c("x","y")
# 
# triang1<-tri.mesh(pts$x,pts$y)
# triang2<-add.constraint(triang1,xs,ys) # coord of the points are added in order
# 
# par(mfrow=c(1,3))
# plot(xs,ys,asp=1)
# polygon(xs,ys)
# plot(triang2,asp=1)
# polygon(xs,ys)
# triOut<-removeOutSegments(triang2)
# plotMesh(vertices=triOut$pts,triangles=triOut$tri)
# polygon(xs,ys)
# # the blue should see not black lines outside

# make triangulation for all polygons in a polygon group
trianglePtsPoly<-function(ptsPolys){
	out<-ptsPolys
	triG<-ptsPolys$tri
	for(i in 1:length(ptsPolys$polys)){
		cat("i:",i,"\n")
		# i<-37
		# get coord polygon
		BoundInG<-ptsPolys$polys[[i]]
		initBoundInG<-BoundInG[1]
		xs<-ptsPolys$pts$x[BoundInG]
		ys<-ptsPolys$pts$y[BoundInG]

		# get coord houses in polygon
		a<-(point.in.polygon(byHouse$X,byHouse$Y,xs,ys)==1)
		pts<-byHouse[which(a),c("X","Y")]
		pts<-pts[!duplicated(pts),]
		names(pts)<-c("x","y")
		nbH<-dim(pts)[1]

		# if not enough houses (less than 3) add fake points
		nadd<-3-nbH
		if(nadd>0){
			dev.new()
			mes<-paste("Please add",nadd,"points in the polygon")
			print(mes)
			plot(xs,ys,main=mes,asp=1)
			pts<-rbind(pts,simplify2array(locator(n=nadd)))
			dev.off()
			nbH<-dim(pts)[1]
		}

		# triangulate houses
		print(pts)
		triang1<-tri.mesh(pts$x,pts$y)

		# test
		plot(xs,ys)
		polygon(xs,ys)
		lines(triang1)

		# add constraint polygon
		triang2<-add.constraint(triang1,xs,ys) # coord of the points are added in order

		# eliminate triangles out of polygon
		triOut<-removeOutSegments(triang2)
		plotMesh(vertices=triOut$pts,triangles=triOut$tri)

		# convert triOut in general coordinates
		previousPtsNum<-max(0,dim(out$pts)[1])
		out$pts<-rbind(out$pts,triOut$pts[1:nbH,])

		# convert ref of boundaries
		adjustBound<-initBoundInG-1-nbH
		boundIntriOut<-triOut$tri>nbH
		numBoundInTriOut<-which(boundIntriOut)
		outtri<-triOut$tri
		outtri[numBoundInTriOut]<-outtri[numBoundInTriOut]+adjustBound
		# expect_equal(range(outtri[numBoundInTriOut]),range(BoundInG))

		# convert refs of points
		numNotBound<-which(!boundIntriOut)
		adjustBound<-previousPtsNum
		outtri[numNotBound]<-outtri[numNotBound]+adjustBound
		triG<-rbind(triG,outtri)
	}
	out$tri<-triG
	return(out)
}
# # initial test
# # get houses of messy polygon 
# a<-(point.in.polygon(byHouse$X,byHouse$Y,xs,ys)==1)
# pts<-byHouse[which(a),c("X","Y")]
# pts<-pts[!duplicated(pts),]
# names(pts)<-c("x","y")
# ptsPolysMessyTri$pts<-rbind(ptsPolysMessyTri$pts,pts)
# 
# triang1<-tri.mesh(pts$x,pts$y)
# triang2<-add.constraint(triang1,xs,ys) # coord of the points are added in order
# 
# par(mfrow=c(1,3))
# plot(xs,ys,asp=1)
# polygon(xs,ys)
# plot(triang2,asp=1)
# polygon(xs,ys)
# triOut<-removeOutSegments(triang2)
# plotMesh(vertices=triOut$pts,triangles=triOut$tri)
# polygon(xs,ys)






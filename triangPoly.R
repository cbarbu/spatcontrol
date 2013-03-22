#=======================
## plotting and format transformation
#=======================
# plot SpatialPolygonsDataFrame in an editable fashion
plotSpatialPolygonsDataFrame<-function(spdf,add=FALSE){
	if(!add){
	plot(NULL,type="n",xlim=spdf@bbox[1,],ylim=spdf@bbox[2,],asp=1)
	}
	for(i in 1:length(spdf)){
		polygon(spdf@polygons[[i]]@Polygons[[1]]@coords)
	}
}
# plotSpatialPolygonsDataFrame(cityBlocks)

# make a list of nodes from coords in poly
coordsPolyFromSpatialPolygonsDataFrame<-function(spdf){
	polys<-list()
	pts<-list()
	polyName<-c()
	for(i in 1:length(spdf)){
		coordsPoly<-spdf@polygons[[i]]@Polygons[[1]]@coords
		firstPt<-max(0,dim(pts)[1])+1
		lastPtInCoords<-dim(coordsPoly)[1]-1 # last duplicate, trash it
		# print(coordsPoly)
		pts<-rbind(pts,as.data.frame(coordsPoly[1:lastPtInCoords,]))
		lastPt<-firstPt+lastPtInCoords-1
		# cat("firstPt",firstPt,"lastPt",lastPt,"\n")
		polys[[i]]<-seq(firstPt,lastPt)
		polyName<-c(polyName,rep(as.character(spdf@data$Name[i]),lastPtInCoords))
	}
	names(pts)<-c("x","y")
	dataPts<-as.data.frame(rep("P",dim(pts)[1]))
	names(dataPts)<-"nature"
	dataPts$polyName<-polyName
	return(list(pts=pts,polys=polys,dataPolys=spdf@data,dataPts=dataPts))
}
#==============================
# Triangulation, general
#==============================
midPoint<-function(pts1,pts2){
	midPts<-c((pts1[1]+pts2[1])/2,(pts1[2]+pts2[2])/2)
	return(midPts)
}
expect_equal(midPoint(c(1,0),c(0,1)),c(0.5,0.5))
midPoints<-function(pts1,pts2,nmid=1){
	pts1<-simplify2array(pts1)
	pts2<-simplify2array(pts2)
	if(nmid<1){
		error("nmid<0")
	}
	# compute the vectors
	factors<-seq(1:(nmid))/(nmid+1)
	midPts<-factors%*%t(pts2-pts1)

	# compute adjust for init
	midPts<-midPts+matrix(rep(pts1,nmid),ncol=2,byrow=TRUE)

	return(midPts)
}
coords<-matrix(c(1,0,0,1),ncol=2)
expect_equal(midPoints(coords[1,],coords[2,]),t(c(0.5,0.5)))
expect_equal(midPoints(c(1,0),c(0,1)),t(c(0.5,0.5)))
exp3<-matrix(c(0.75,0.25,0.50,0.50,0.25,0.75),ncol=2,byrow=TRUE)
expect_equal(midPoints(c(1,0),c(0,1),3),exp3,3)


# insert in vector at spot i, moving i if existing to the right
insert<-function(vals,vect,spot){
	if(spot>1){
		if(spot>length(vect)){
			out<-c(vect,vals)
			if(spot>length(vect)+1){
				warning(paste("insert at spot",spot,"but length is",length(vect)))
			}
		}else{
			out<-c(vect[1:(spot-1)],vals,vect[spot:length(vect)])
		}
	}else{
		out<-c(vals,vect)
	}
	return(out)
}
expect_equal(insert(c(1,2),c(5,4,3,2,1),2),c(5,1,2,4,3,2,1))
expect_equal(insert(c(1,2),c(5,4,3,2,1),1),c(1,2,5,4,3,2,1))
expect_equal(insert(c(1,2),c(5,4,3,2,1),6),c(5,4,3,2,1,1,2))
suppressWarnings(expect_equal(insert(c(1,2),c(5,4,3,2,1),7),c(5,4,3,2,1,1,2)))
expect_warning(insert(c(1,2),c(5,4,3,2,1),7),paste("insert at spot 7 but length is 5"))

# add intermediary points if necessary 
# between two given nodes limNodes
# in polygon pol 
# of structure inMap
# and return corresponding structure outMap
autoAddInterNodes<-function(inMap,pol,side,tr=40){
	outMap<-inMap
	initNodes<-outMap$polys[[pol]]
	initNodes<-c(initNodes,initNodes[1])
	coords<-outMap$pts[initNodes[side:(side+1)],]
	# print(coords)
	d<-dist(coords)
	nmid<-0
	if(d>tr){
		nmid<-as.vector(d%/%tr)
		# coord midpoints
		mid<-midPoints(coords[1,],coords[2,],nmid)
		# add the coordinates
		dmid<-outMap$dataPts[rep(initNodes[1],nmid),]
		outMap<-addPoints(outMap,mid,dataPts=dmid)
		lastmid<-dim(outMap$pts)[1]
		refmid<-seq((lastmid-nmid+1),lastmid)

		# add the node(s) in polygon
		outMap$polys[[pol]]<-insert(refmid,outMap$polys[[pol]],side+1)
	}
	attributes(outMap)$nmid<-nmid
	# cat("side:",side,"refmid:",refmid,"d:",d,"\n")
	return(outMap)
}
blurb<-list()
blurb$pts<-matrix(c(1,0,0,1,0,2),ncol=2,byrow=TRUE)
blurb$dataPts<-as.data.frame(matrix(seq(1,9),ncol=3))
blurb$polys<-list()
blurb$polys[[1]]<-c(1,2,3)
burp<-autoAddInterNodes(blurb,1,1,tr=0.4)
expect_equal(burp$pts[1,],blurb$pts[1,])
expect_equal(burp$polys[[1]],c(1,4,5,6,2,3))
expect_equal(burp$pts[1,],blurb$pts[1,])
exp3<-matrix(c(0.75,0.25,0.50,0.50,0.25,0.75),ncol=2,byrow=TRUE)
expect_equal(burp$pts[4:6,],exp3)
expect_equal(simplify2array(burp$dataPts[4,]),simplify2array(blurb$dataPts[1,]))
expect_equal(simplify2array(burp$dataPts[5,]),simplify2array(blurb$dataPts[1,]))
expect_equal(simplify2array(burp$dataPts[6,]),simplify2array(blurb$dataPts[1,]))


# add intermediary points if necessary at border of polygons
autoAddInterPoly<-function(ptsPolys,threshold=40){
	out<-ptsPolys
	for(pol in 1:length(ptsPolys$polys)){
		initNodes<-ptsPolys$polys[[pol]]
		prevnode<-1
		for(side in 1:length(initNodes)){
			# cat("side:",side,"prevnode:",prevnode,"\n")
			out<-autoAddInterNodes(out,pol,prevnode,tr=threshold)
			prevnode<-prevnode+attributes(out)$nmid+1
		}
	}
	return(out)
}
coords<-matrix(c(2,2,2,1),ncol=2,byrow=TRUE)
dmid<-blurb$dataPts[rep(1,2),]
blurb<-addPoints(blurb,coords,dataPts=dmid)
blurb$polys[[1]]<-(1:5)
burp<-autoAddInterPoly(blurb,0.4)
expect_equal(burp$polys[[1]],c(1,6,7,8,2,9,10,3,11,12,13,14,4,15,16,5,17,18,19))
# visual check
# par(mfrow=c(1,2))
# plotMesh(vertices=blurb$pts,polys=blurb$polys)
# plotMesh(vertices=burp$pts,polys=burp$polys)

#==============================
# Triangulation of polygons
#==============================

#--------------------
# gpclib approach: pb, does whatever triangulation
# and add points to the limit of the polygons
#--------------------

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
# # test
# x <-structure(c(0.0934073560027759, 0.192713393476752, 0.410062456627342,0.470020818875781, 0.41380985426787, 0.271408743927828, 0.100902151283831,
# 0.0465648854961832, 0.63981588032221, 0.772382048331416,0.753739930955121, 0.637744533947066, 0.455466052934407,0.335327963176065, 0.399539700805524,0.600460299194476), .Dim = c(8, 2))
# y <- structure(c(0.404441360166551, 0.338861901457321, 0.301387925052047, 
# 		    0.404441360166551, 0.531852879944483, 0.60117973629424, 0.625537820957668, 
# 		    0.179976985040276, 0.341542002301496, 0.445109321058688,
# 		    0.610817031070196, 0.596317606444189, 0.459608745684695,
# 		    0.215189873417722), .Dim = c(7, 2))
# 
# x1 <- as(x, "gpc.poly")
# y1 <- as(y, "gpc.poly")
# ## Show the triangulation
# par(mfrow=c(1,2))
# plot(append.poly(x1, y1))
# triangles <- triangulate(append.poly(x1,y1))
# for (i in 0:(nrow(triangles)/3 - 1)) 
# 	polygon(triangles[3*i + 1:3,], col="lightblue")
# 
# # transform to our format
# coordsxy<-rbind(x,y)
# colnames(coordsxy)<-c("x","y")
# tri<-gpctri2tri(triangles,coordsxy)
# plotMesh(vertices=coordsxy,triangles=tri)

#------------------------
# triangulation based on constrained delaunay, with houses
# use tripack
#------------------------
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
	if(tpk$nc>0){ # if no constraint, no check
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
trianglePtsPoly<-function(ptsPolys,byHouse){
	out<-ptsPolys
	triG<-ptsPolys$triIntra
	cat("triangulate poly:")
	cityBlocks$id_manz<-""
	for(i in 1:length(ptsPolys$polys)){
		cat("i:",i," ")
		# i<-37
		# get coord polygon
		BoundInG<-ptsPolys$polys[[i]]
		initBoundInG<-BoundInG[1]
		xs<-ptsPolys$pts$x[BoundInG]
		ys<-ptsPolys$pts$y[BoundInG]

		# get coord houses in polygon
		a<-(point.in.polygon(byHouse$X,byHouse$Y,xs,ys)==1)
		byHouse[which(a),"id_manz2"]<-ptsPolys$dataPolys[i,"Name"]
		cityBlocks$id_manz[i]<-byHouse[which(a)[1],"id_manz"]

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
			polygon(xs,ys)
			lines(pts[,"x"],pts[,"y"],type="p")
			pts<-rbind(pts,simplify2array(locator(n=nadd)))
			dev.off()
			nbH<-dim(pts)[1]
		}

		# triangulate houses
		# print(pts)
		triang1<-tri.mesh(pts$x,pts$y)

		# # test
		# plot(xs,ys)
		# polygon(xs,ys)
		# lines(triang1)

		# add constraint polygon
		triang2<-add.constraint(triang1,xs,ys) # coord of the points are added in order

		# eliminate triangles out of polygon
		triOut<-removeOutSegments(triang2)
		# plotMesh(vertices=triOut$pts,triangles=triOut$tri)

		# convert triOut in general coordinates
		# pts
		previousPtsNum<-max(0,dim(out$pts)[1])
		out$pts<-rbind(out$pts,triOut$pts[1:nbH,])

		# dataPts
		temp<-as.data.frame(cbind(rep("H",nbH),as.character(rep(out$dataPolys$Name[i]))))
		names(temp)<-names(out$dataPts)
		out$dataPts<-rbind(out$dataPts,temp)

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
	byHouseLim<<-byHouse
	cityBlocks<<-cityBlocks
	out$triIntra<-triG
	cat("\n")
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

#============================
# Triangulation of streets
#============================
#----------------------------
# Adding points/triangles
#----------------------------
addPoints<-function(ptsPolys,coords,dataPts=NULL){
	# ptsPolys: list with all points, polys etc...
	# coords: coordinates of the point
	# dataPts: data of the point to put in dataPts

	# coordinates
	ptsPolys$pts<-rbind(ptsPolys$pts,coords)

	# dataPts
	nbColDataPts<-dim(ptsPolys$dataPts)[2]
	if(is.null(dataPts)){
		toAdd<-as.data.frame(t(c("S",rep(NA,nbColDataPts-1))))
		names(toAdd)<-names(ptsPolys$dataPts)
	}else{
		if(dim(dataPts)[2]!=nbColDataPts){
			error("dataPts has not good number of columns")
		}
		toAdd<-dataPts
	}

	ptsPolys$dataPts<-rbind.general(ptsPolys$dataPts,toAdd)

	return(ptsPolys)
}

#----------------------------
# Calculus
#----------------------------
# are node and connected nodes in vect of refs
# nodes: nodes ref
# nodenum: item in nodes
# vect: ref of nodes on border
arcInVect<-function(nodenum,nodes,vect){
	# is border?
	# previous
	nodePrev<-nodenum-1
	if(nodePrev<1){ nodePrev<-length(nodes)}
	# next
	nodeNext<-nodenum+1
	if(nodeNext>length(nodes)){ nodeNext<-1}
	# test
	onBorder<-nodes[c(nodePrev,nodenum,nodeNext)] %in% vect
	return(onBorder)
}
# to evaluate which arc is the closest to a set of points
# make the meanOverRef(meanOverRef1(dist(pts,pts1)))
# and same for 2 and the smallest win
dists<-function(pts1,pts2){
	pivot<-dim(pts1)[1]
	totpts<-rbind(pts1,pts2)
	npts<-dim(totpts)[1]

	ds<-as.matrix(dist(totpts))[(1:pivot),((pivot+1):npts)]

	return(ds)
}
# Test
coordsRef<-matrix(c(0,0,0,1),ncol=2)
coords1<-matrix(c(0.5,0,0.5,0.5,0.5,1),ncol=2,byrow=TRUE)
res<-matrix(c(0.500000,0.7071068,1.118034,1.118034,0.7071068,0.500000),ncol=3,byrow=TRUE)
expect_true(mean(abs(res-dists(coordsRef,coords1)))<10e-6)

closerSet<-function(coordsRef,coords1,coords2){
	if(mean(dists(coordsRef,coords1))<=mean(dists(coordsRef,coords2))){
		return(1)
	}else{
		return(2)
	}
}
# Test
coordsRef<-matrix(c(0,0,0,1),ncol=2)
coords1<-matrix(c(0.5,0,0.5,0.5,0.5,1),ncol=2,byrow=TRUE)
coords2<-coords1
coords2[2,1]<-1
# plot(coordsRef,type="l")
# lines(coords1,col="orange")
# lines(coords2,col="blue")
expect_equal(closerSet(coordsRef,coords1,coords1),1)
expect_equal(closerSet(coordsRef,coords1,coords2),1)
expect_equal(closerSet(coordsRef,coords2,coords1),2)

# return the two arcs of nodes
arcsOfNodes<-function(nodes,node1,node2){
	arc1<-nodes[seq(node1,node2)]
	arc2<-nodes[c(seq(node2,length(nodes)),seq(1,node1))]
	return(list(arc1,arc2))
}
nodes<-c(3,1,2,44,4)
arc1<-c(1,2,44)
arc2<-c(44,4,3,1)
arcs<-arcsOfNodes(nodes,2,4)
expect_equal(arcs[[1]],arc1)
expect_equal(arcs[[2]],arc2)

# apply the attribute to 
initEndCross<-function(ptsPoly,cross){
	if(attributes(cross)$reverse){
		return(ptsPoly$crossings[cross,3:2])
	}else{
		return(ptsPoly$crossings[cross,2:3])
	}
}
blurb<-list()
blurb$crossings<-matrix(c(1231,1,2),ncol=3)
crrr<-1
attributes(crrr)$reverse<-FALSE
expect_equal(initEndCross(blurb,crrr),c(1,2))
attributes(crrr)$reverse<-TRUE
expect_equal(initEndCross(blurb,crrr),c(2,1))

# check if intermediary nodes on closest arc
# rely on the facts that 
# - attributes of crossing are corrects
# - nodes are actually on opposite polygons
# and then inits are consecutive
intermNodes<-function(ptsPolys,formerCross,currentCross){
	# get the nodes to test
	initEndForm<-initEndCross(ptsPolys,formerCross)
	initEndCur<-initEndCross(ptsPolys,currentCross)
	sysNodes<-c(initEndForm[1],initEndCur[1])
	transNodes<-c(initEndForm[2],initEndCur[2])

	# are they contiguous on arc closest to initNodes
	
	return(nodesPoly)
}

# check if exist intermediary points between two nodes
# on the 
# triangle a street portion between formerTriplet and current Triplet
triangle_street_portion<-function(ptsPoly,formerCrossing,currentCrossing){
	# recursively triangle between former and current
	# according to opposite intermediate
	while(countInterm(ptsPoly,formerCrossing,currentCrossing))

	# triangle former and current same block
	# to middle

	# triangle former and current middle
	# to previous same

	# create as many middle betwen current and former
	# as intermediary points on the opposite side
	# for op in intermediary points

	## for the following may want to test first what is the best 
	## triangulation: shortest line deviding the quadrangle
	# triangle previous op and op and previous middle  

	# triangle previous and current middle and op
	# to middle (if previous/current opposite different)

	# switch currents/previous

	return(ptsPoly)
}
# determine where to make the triangulation and if additional
# streets points are needed
# assume previous filter on the two nodes not being on the limits
additional_point_corner<-function(){
	# 

	# find best corner: between former and current which
	# one is closest to the block it is not already 
	# connected with

	# add middle corner to second cb

	# triangle missing part of street if corner is former
	if(cornerFormer){
		triangle_street_portion(cornerTriplet,currentTriplet)
	}else{
		triangle_street_portion(cornerTriplet,previousTriplet)
	}
	# add intersection
	return(list(triangles,points))
}

close_intersection<-function(nodes,triangles,middle_point,intersectionRefNode){
	# intersectionRefNode is the node that allowed to first flag the intersection

	currentNode<-intersectionRefNode
	while(! backToRef){
		# find next node clock wise 
		# triangle middle_point, nextnode, currentnode

		# 
		# find next node clock wise 

		backToRef<-nextnode==intersectionRefNode
	}

}
get_middle_point<-function(){
	intersectionOK<-FALSE
	while(!intersectionOK){
	# display the intersection

	# offer to select a point or enter a new one
	identify()

	# close the gaps around the choose middle_point 
	close_intersection(intersection,middle_point)

	# is it ok ?

	# want to retry?
	}
}



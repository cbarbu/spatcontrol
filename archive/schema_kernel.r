graphics.off()
source("functions_intercept.r")

extent<-105
T<-0.25
f<-15
Ku<-0.4 ### CAREFULL hereafter, doesn't take into account Ku and the normalizing factor. If it was, the correlations could be directly correlated. Has to be done before publication of this images
threshold<-sqrt(2*extent^2)*1.5
Kernel<-expKernel
steps<-300 # the more the more precise careful, 100 is already a lot
ref<-c(-5,-5)
maxz<-Kernel(1,12,f) # ceiling at the value of weight for mean distance of first neighbor inside block
border=NA # set cell borders in 3D graphs: NULL->border ; NA-> no limits
dark<-0.6
with.3D<-FALSE

# to plot houses and limits of a city block
# base polygon
base.x<-c(-35,0,40,5)
base.y<-c(-20,15,-25,-60)
# base houses in the block
ref.h.x<-c(0,12.5,25)
ref.h.y<-c(0,-12.5,-25)
plot.block<-function(xshift,yshift,housescol){
	# base houses coordinates
	ref.m.x<-ref.h.x[c(1,3)]
	ref.m.y<-ref.h.y[c(1,3)]
	# base polygon coordinates 

	polygon(base.x+xshift,base.y+yshift)
	lines(ref.h.x+xshift,ref.h.y+yshift,type="p",pch=5,col=housescol)
	lines(ref.h.x-20+xshift,ref.h.y-20+yshift,type="p",pch=5,col=housescol)
	lines(ref.m.x-10+xshift,ref.m.y-10+yshift,type="p",pch=5,col=housescol)
}
fullMapPlot<-function(initPlot=FALSE){
	if(initPlot==TRUE){
		plot(c(ref[1]),c(ref[2]),type="p",pch=3,asp=1,xlab="x (m)",ylab="y (m)",xlim=c(-extent,extent),ylim=c(-extent,extent)) 
	}
	for(diagshift in seq(-125,125,40)){
		for(antidiagshift in seq(-125,165,45)){
			xshift<-15+diagshift-antidiagshift
			yshift<--5+diagshift+antidiagshift
			if(xshift != 0 && yshift !=0){
				plot.block(ref[1]+xshift,ref[2]+yshift,1)
			}
		}
	}
	plot.block(ref[1],ref[1],4)
	points(c(ref[1]),c(ref[2]),type="p",pch=3)
}
# fullMapPlot()


# get the grid (as if all out of block)
out<-grid.from.kernel(ref[1],ref[2],10,Kernel,T=T,f,steps=steps,xlim=c(-extent,extent),ylim=c(-extent,extent),tr=threshold)

# modify inside block points
dist.weight<-matrix(as.vector(as.matrix((out$raw.weights))),nrow=length(out$xs))
pointsInPoly<-which(point.in.polygon(out$x,out$y,base.x+ref[1],base.y+ref[2])==1)
dist.weight[pointsInPoly]<-dist.weight[pointsInPoly]/T
dist.weight[dist.weight>maxz]<-maxz

####  plotting
# # 3D visualization with vrmlgen:
# cloud3d(out$x,out$y,dist.weight)
# # then open with ouvre out.wrl

# color plotting with plot_reel
dev.new(width=4,height=4)
par(mar = c(4,4,0.5,0.5))
zcol<-plot_reel(out$x,out$y,dist.weight,base=0,top=1)
# 3D ploting with persp (base)
dev.new(width=4,height=4)
par(mar = c(0.5,0.5,0.5,0.5))
pmat<-persp(out$xs,out$ys,dist.weight,asp=1,zlim=c(0,maxz),
	col=make.col.persp(dist.weight,color.function=heat.colors),
	phi=20,theta=-30,
	border=border, # set cell borders: NULL->border ; NA-> no limits
	xlab="x",ylab="y",zlab="correlation")
# lines(trans3d(c(ref[1],50),c(ref[2],50),c(1,0),pmat),col=4,lty=2)

# polygon(trans3d(base.x,base.y,rep(0,length(base.x)),pmat))

dev.new(width=4,height=4)
par(mar = c(4,4,0.5,0.5))
image(x=out$xs,y=out$ys,z=dist.weight,asp=1,col=heat.colors(100),xlab="x (m)",ylab="y (m)")
#plot of the transects
# lines(c(ref[1],ref[1]+60),c(ref[2],ref[2]+60),lty=1,col=4)
lines(c(ref[1]-80,ref[1]+80),c(ref[2]+80,ref[2]-80),lty=1,col=4)
fullMapPlot()

# 2D plot blue transect
dev.new(width=4,height=4)
par(mar = c(4,4,0.5,0.5))
xabs<-seq(0,extent,0.5)
distHSB<-sqrt(2*10^2) #dist to first house on transect
distHAS<-sqrt(2*20^2) #dist to second house on transect
limit<-(distHAS+distHSB)/2 # 12.5

yabs<-Kernel(1,xabs,f)
plot(xabs,yabs,type="l",xlab="distance (m)",ylab="correlation",lty=2,col=4,ylim=c(0,maxz))
yabs<-Kernel(T,xabs,f)
lines(xabs,yabs,type="l",lty=2,col=4)

yabs[xabs>limit]<-Kernel(T,xabs[xabs>limit],f)
yabs[xabs<=limit]<-Kernel(1,xabs[xabs<=limit],f)
lines(xabs,yabs,type="l",lty=1,col=4,ylim=c(0,maxz))

xH<-c(distHSB,distHAS,distHAS+distHSB)
yH<-c(Kernel(1,distHSB,f),Kernel(T,distHAS,f),Kernel(T,distHSB+distHAS,f))
points(xH,yH,pch=22)

# 2D plot green transect
dev.new(width=4,height=4)
par(mar = c(4,4,0.5,0.5))
xabs<-seq(-extent,extent,0.5)
limit1<- -(distHAS)/2 # 12.5
limit2<- (4*distHSB+distHAS)/2 # 12.5

yabs<-Kernel(1,abs(xabs),f)
plot(xabs,yabs,type="l",xlab="distance (m)",ylab="correlation",lty=2,col=4,ylim=c(0,maxz),lwd=2,cex=2,cex.lab=1.5)
yabs<-Kernel(T,abs(xabs),f)
lines(xabs,yabs,type="l",lty=2,col=4,lwd=2)

yabs[xabs<=limit1]<-Kernel(T,abs(xabs[xabs<=limit1]),f)
yabs[xabs>limit1&xabs<limit2]<-Kernel(1,abs(xabs[xabs>limit1&xabs<limit2]),f)
lines(xabs,yabs,type="l",lty=1,col=4,ylim=c(0,maxz),lwd=2)

xH<-c(-(distHAS+2*distHSB),-distHAS,0,distHSB,2*distHSB,2*distHSB+distHAS,4*distHSB+distHAS)
yH<-0*xH
yH[xH<=limit1|xH>=limit2]<-Kernel(T,abs(xH[xH<=limit1|xH>=limit2]),f)
yH[xH>limit1&xH<limit2]<-Kernel(1,abs(xH[xH>limit1&xH<limit2]),f)
points(xH,yH,pch=22,cex=2) 

### schema Moran
dev.new(width=4,height=4)
par(mar = c(4,4,0.5,0.5))
plot(c(ref[1]),c(ref[2]),type="p",pch=3,asp=1,xlab="x (m)",ylab="y (m)",xlim=c(-extent,extent),ylim=c(-extent,extent)) 
# plot a ring of neighborhood
library(plotrix)
draw.circle(ref[1],ref[2],45,col="lightgrey")
draw.circle(ref[1],ref[2],25,col="white")
fullMapPlot()

## schema indirect
# according to fit in 
# ~/Documents/recherche/Penn/anal_bug_spread_cross_sectionnal/corentin/R_analysis/streets_impact/paucarpata/start5/~/Documents/recherche/Penn/anal_bug_spread_cross_sectionnal/corentin/R_analysis/streets_impact/paucarpata/start5/spat_autocorrelogram_struct.R
# initially think in ~/Documents/recherche/Penn/anal_bug_spread_cross_sectionnal/corentin/R_analysis/general/encuesta_mariano_melg/spat_autcorrel_par_bande.R
distances<- seq(10,1010,10)
intra_prop<-function(distances){
	Rate<-1.9*exp(-(distances+5)/28)
	Rate[Rate>1]<-1
	return(Rate);
}

# sdCorrel<-200
# intraBlockDeparture<-185
# interBlockDeparture<-37
fSB<-80
fAS<-120
propSB<-0.42
propAS<-0.22
indirect.correl<-function(T,distances,f){
	# f and T are ignored, only used for consistency with Kernel
	if(class(distances)=="spam"){
		correl<-distances
		# correl@entries<-dnorm(distances@entries,0,sdCorrel)*(intra_prop(distances@entries)*intraBlockDeparture+(1-intra_prop(distances@entries))*interBlockDeparture)
		correl@entries<-exp(-distances@entries/fSB)*propSB*intra_prop(distances@entries)+exp(-distances@entries/fAS)*propAS*(1-intra_prop(distances@entries))
	}else{
		# correl<-dnorm(distances,0,sdCorrel)*(intra_prop(distances)*intraBlockDeparture+(1-intra_prop(distances))*interBlockDeparture)
		correl<-exp(-distances/fSB)*propSB*intra_prop(distances)+exp(-distances/fAS)*propAS*(1-intra_prop(distances))
	}
	return(correl)
}
dist.classes<-c(0,seq(5,1.5*sqrt(2*extent^2),20))
make.dist.bands<-function(distances){
	if(class(distances)=="spam"){
		sel<-which(distances@entries<dist.classes[2])
		distances@entries[sel]<-0
		for(i in 2:length(dist.classes-1)){
			medDist<-(dist.classes[i+1]+dist.classes[i])/2
			sel<-which(distances@entries<dist.classes[i+1]&distances@entries>=dist.classes[i])
			distances@entries[sel]<-medDist
		}
	}else{
		distances[distances<dist.classes[2]]<-0
		for(i in 2:length(dist.classes-1)){
			medDist<-(dist.classes[i+1]+dist.classes[i])/2
			sel<-which(distances<dist.classes[i+1]&distances>=dist.classes[i])
			distances[sel]<-medDist
		}
	}
	return(distances)
}
indirect.correl.band<-function(T,distances,f){
	# f and T are ignored, only used for consistency with Kernel
	# makes bands of distance
	
	distances<-make.dist.bands(distances)

	# apply normal indirect.correl
	correl<-indirect.correl(T,distances,f)

	return(correl)
}

# 1D plot
dev.new(width=4,height=4)
plot(distances,indirect.correl(f,distances,T),col=8,type="l")

# Kernel<-indirect.correl
Kernel<-indirect.correl.band

# get the grid (as if all out of block)
out<-grid.from.kernel(ref[1],ref[2],10,Kernel,T=T,f,steps=steps,xlim=c(-extent,extent),ylim=c(-extent,extent),tr=threshold)

# get a matrix from the output
dist.weight<-matrix(as.vector(as.matrix((out$raw.weights))),nrow=length(out$xs))

# 3D plot
dev.new(width=4,height=4)
par(mar = c(0.5,0.5,0.5,0.5))
pmat<-persp(out$xs,out$ys,dist.weight,asp=1,zlim=c(0,maxz),
	col=make.col.persp(dist.weight,color.function=heat.colors),
	phi=20,theta=-30,
	border=border, # set cell borders: NULL->border ; NA-> no limits
	xlab="x",ylab="y",zlab="correlation")
# lines(trans3d(c(ref[1],50),c(ref[2],50),c(1,0),pmat),col=4,lty=2)

# nice 2D plot
library(fields)
dev.new(width=4,height=4)
par(mar = c(4,4,0.5,0.5))
image.plot(x=out$xs,y=out$ys,z=dist.weight,asp=1,col=heat.colors(100),xlab="x (m)",ylab="y (m)",zlim=c(0,max(dist.weight)))

#### final comparison plot general gaussian/field streets
dev.new(width=9.2,height=4)

# Here is quick but quirky way to add a common legend to several plots. 
# The idea is leave some room in the margin and then over plot in this margin

par(oma=c( 0,0,0,4)) # margin of 4 spaces width at right hand side
par(mar = c(4,5,0.5,0.5))
set.panel( 1,2) # 2X2 matrix of plots

# now draw all your plots using usual image command
# Kernel<-indirect.correl
Kernel<-indirect.correl.band
out<-grid.from.kernel(ref[1],ref[2],10,Kernel,T=T,f,steps=steps,xlim=c(-extent,extent),ylim=c(-extent,extent),tr=threshold)
dist.weight1<-matrix(as.vector(as.matrix((out$raw.weights))),nrow=length(out$xs))
dist.weight1[dist.weight1>maxz]<-maxz

image(x=out$xs,y=out$ys,z=dist.weight1,asp=1,col=heat.colors(100),xlab="x (m)",ylab="y (m)",zlim=c(0,max(dist.weight1)))
fullMapPlot()

Kernel<-expKernel
out<-grid.from.kernel(ref[1],ref[2],10,Kernel,T=T,f,steps=steps,xlim=c(-extent,extent),ylim=c(-extent,extent),tr=threshold)

dist.weight2<-matrix(as.vector(as.matrix((out$raw.weights))),nrow=length(out$xs))
dist.weight2[pointsInPoly]<-dist.weight2[pointsInPoly]/T
dist.weight2[dist.weight2>maxz]<-maxz

image(x=out$xs,y=out$ys,z=dist.weight2,asp=1,col=heat.colors(100),xlab="x (m)",ylab="y (m)",zlim=c(0,max(dist.weight1,dist.weight2)))
fullMapPlot()

# plot the legend
par(oma=c( 0,0,0,1))# reset margin to be much smaller.
image.plot( legend.only=TRUE,col=heat.colors(100),zlim=c(0,max(dist.weight1,dist.weight2))) 

# image.plot tricked into  plotting in margin of old setting 
set.panel() # reset plotting device


#### structured morans'I
# using ~/Documents/recherche/Penn/anal_bug_spread_cross_sectionnal/corentin/R_analysis/streets_impact/paucarpata/start5/~/Documents/recherche/Penn/anal_bug_spread_cross_sectionnal/corentin/R_analysis/streets_impact/paucarpata/start5/spat_autocorrelogram_struct.R

# Kernel<-indirect.correl
Kernel<-indirect.correl.band
out<-grid.from.kernel(ref[1],ref[2],10,Kernel,T=T,f,steps=steps,xlim=c(-extent,extent),ylim=c(-extent,extent),tr=threshold)

dists<-matrix(as.vector(as.matrix((out$dists))),nrow=length(out$xs))
dists<-make.dist.bands(dists)
dist.weight3<-exp(-dists/fAS)*propAS
dist.weight3[pointsInPoly]<-exp(-dists[pointsInPoly]/fSB)*propSB
dist.weight3[dist.weight3>maxz]<-maxz

image.plot(x=out$xs,y=out$ys,z=dist.weight3,asp=1,col=heat.colors(100),xlab="x (m)",ylab="y (m)",zlim=c(0,max(dist.weight3)))

#### final comparison plot general moran/structured moran/field streets
if(with.3D){
dev.new(width=12,height=8)
layout(matrix(c(1,3,5,7,2,4,6,7),2,4,byrow=TRUE),widths=c(1,1,1,0.2))
}else{
	# # vert
	# dev.new(width=4.85,height=12)
	# layout(matrix(c(1,2,3,5,4,6),3,2,byrow=FALSE),widths=c(1,0.2))
	# horiz
	dev.new(width=11,height=4)
	layout(matrix(c(1,2,3,4),1,4,byrow=FALSE),widths=c(1,1,1,0.15))
}

# Here is quick but quirky way to add a common legend to several plots. 
# The idea is leave some room in the margin and then over plot in this margin

par(oma=c( 0,0,0,4),mar=c(5,5,0.5,0.5)) # margin of 4 spaces width at right hand side
# par(mar = c(4,5,0.5,0.5),mfcol=c(2,3))
# set.panel( 2,3) # 2X2 matrix of plots

## now draw all your plots using usual image command
# general Moran's I
image(x=out$xs,y=out$ys,z=dist.weight1,asp=1,col=heat.colors(100),xlab="x (m)",ylab="y (m)",zlim=c(0,max(dist.weight1,dist.weight2,dist.weight3)))
fullMapPlot()

if(with.3D){
pmat<-persp(out$xs,out$ys,dist.weight1,asp=1,zlim=c(0,maxz),
	col=make.col.persp(dist.weight1,color.function=heat.colors),
	shade=dark,
	phi=20,theta=-30,
	border=border, # set cell borders: NULL->border ; NA-> no limits
	xlab="x",ylab="y",zlab="correlation")
}

# Moran's I with streets
image(x=out$xs,y=out$ys,z=dist.weight3,asp=1,col=heat.colors(100),xlab="x (m)",ylab="y (m)",zlim=c(0,max(dist.weight1,dist.weight2,dist.weight3)))
fullMapPlot()

if(with.3D){
pmat<-persp(out$xs,out$ys,dist.weight3,asp=1,zlim=c(0,maxz),
	col=make.col.persp(dist.weight3,color.function=heat.colors),#"white" ,
	shade=dark,
	phi=20,theta=-30,
	border=border, # set cell borders: NULL->border ; NA-> no limits
	xlab="x",ylab="y",zlab="correlation")
}

# Gaussian field with streets
image(x=out$xs,y=out$ys,z=dist.weight2,asp=1,col=heat.colors(100),xlab="x (m)",ylab="y (m)",zlim=c(0,max(dist.weight1,dist.weight2,dist.weight3)))
Kernel<-expKernel
fullMapPlot()

if(with.3D){
pmat<-persp(out$xs,out$ys,dist.weight2,asp=1,zlim=c(0,maxz),
	col=make.col.persp(dist.weight2,color.function=heat.colors),
	shade=dark,
	phi=20,theta=-30,
	border=border, # set cell borders: NULL->border ; NA-> no limits
	xlab="x",ylab="y",zlab="correlation")
}

# plot the legend
par(oma=c( 0,0,0,0))# reset margin to be much smaller.
set.panel() # reset plotting device
image.plot( legend.only=TRUE,col=heat.colors(100),zlim=c(0,max(dist.weight1,dist.weight2,dist.weight3))) 

set.panel() # reset plotting device

## Nota: to get this images to display well in a pdf is a challenge do get it write:

# dev.print(device=pdf,"compared_2D_3D.pdf")
## in bash
# plosization.sh compared_2D_3D.pdf
## then convert it back to pdf
# convert compared_2D_3D.tiff ../../presentation_ASTMH/Images/compared_2D_3D.pdf


# only Moran
dev.new(width=4.85,height=7)
layout(matrix(c(1,2,3,3),2,2,byrow=FALSE),widths=c(1,0.2))
par(oma=c( 0,0,0,4),mar=c(4,4,0.5,0.5)) # margin of 4 spaces width at right hand side
image(x=out$xs,y=out$ys,z=dist.weight1,asp=1,col=heat.colors(100),xlab="",ylab="y (m)",zlim=c(0,max(dist.weight1,dist.weight2,dist.weight3)))
fullMapPlot()
image(x=out$xs,y=out$ys,z=dist.weight3,asp=1,col=heat.colors(100),xlab="x (m)",ylab="y (m)",zlim=c(0,max(dist.weight1,dist.weight2,dist.weight3)))
fullMapPlot()
par(oma=c( 0,0,0,1))# reset margin to be much smaller.
set.panel() # reset plotting device
image.plot( legend.only=TRUE,col=heat.colors(100),zlim=c(0,max(dist.weight1,dist.weight2,dist.weight3))) 

# dev.print(device=pdf,"2DMoran.pdf")


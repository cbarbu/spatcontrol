# get a map of the "subpopulations" by smoothing and countour plot

steps<-300
# inverse linear kernel
Kernel<-expKernel
out<-grid.from.kernel(block_data$easting,block_data$northing,est.prob.b.pos,Kernel,T=T,f,steps=steps,tr=threshold)

dist.weight<-matrix(as.vector(as.matrix(out$z)),nrow=length(out$xs))
image(x=out$xs,y=out$ys,z=dist.weight,asp=1,col=heat.colors(100),xlab="x (m)",ylab="y (m)")

contour(x=out$xs,y=out$ys,z=dist.weight,add=T)



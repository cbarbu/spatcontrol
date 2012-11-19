
library("spam")
library("fields")
source("R/grid.from.kernel.R")
# test
known.x<-c(2.2,3.3,7.5,7)
known.y<-c(2.3,7.4,8.5,4)
known.z<-c(1,-3,10,-5)

out<-grid.from.kernel(known.x,known.y,known.z,steps=100)
par(mfrow=c(2,3))
plot_reel(out$x,out$y,out$z)
lines(known.x,known.y,col=4,pch=1,type="p")
plot_reel(known.x,known.y,known.z)

contour(out$xs,out$ys,matrix(out$z,nrow=length(out$xs)),asp=1)
lines(known.x,known.y,col=4,pch=1,type="p")
# filled.contour(out$xs,out$ys,matrix(out$z,nrow=length(out$xs)),asp=1,nlevels=50)
# lines(known.x,known.y,col=4,pch=1,type="p")

image(x=out$xs,y=out$ys,z=matrix(out$z,nrow=length(out$xs)),asp=1,col=gray((0:100)/100))
lines(known.x,known.y,col=4,pch=1,type="p")
image(x=out$xs,y=out$ys,z=matrix(out$z,nrow=length(out$xs)),asp=1,col=heat.colors(100))
contour(out$xs,out$ys,matrix(out$z,nrow=length(out$xs)),asp=1,add=T)
# using "fields" to plot the corresponding legend
image.plot(x=out$xs,y=out$ys,z=matrix(out$z,nrow=length(out$xs)),asp=1,col=heat.colors(100))

# other good colors are:
# col= topo.colors(100) # to have bleue for negative, green->red for others
# col=terrain.colors(100) # have green to red/white



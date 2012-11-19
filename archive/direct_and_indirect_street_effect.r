# pull apart and combine the direct and indirect effect of streets
# after the final results of a simulation
### especially relies on:
# meanT<-mean(smallsampled$T)
# meanf<-mean(smallsampled$f)
ASDmat<-as.spam(Dmat*AS)
SBDmat<-as.spam(SB*Dmat)
par(mfrow=c(1,2))

# mean distance of the first neighbourg across street and inside the same block
mean_min_SB<-mean(apply_by_row_not_null.spam(SBDmat,min))
mean_min_AS<-mean(apply_by_row_not_null.spam(ASDmat,min))
cat(" mean distance of the first neighbourg, same block:",mean_min_SB," surrounding blocks",mean_min_AS,"\n");

# ratio influence first neighbourg
mean_ratio_first_neigh<-mean(meanT*exp(-apply_by_row_not_null.spam(ASDmat,min)/meanf)/exp(-apply_by_row_not_null.spam(SBDmat,min)/meanf))

# representation on the kernel
xabs<-seq(0,threshold,0.5)
titre<-paste("A-mean distance of the first neighbourg\n in relation to the kernel shape (f=",round(meanf,2),"; T=",round(meanT,2),")")
plot(xabs,exp(-xabs/f),type="l",col=4,xlab="Distance (m)",ylab=expression(omega[ij]),main=titre)
lines(xabs,T*exp(-xabs/f),type="l",col=2)
lines(c(mean_min_SB,mean_min_SB),c(-0.5,exp(-mean_min_SB/f)),col=4,lty=2)
lines(c(mean_min_AS,mean_min_AS),c(-0.5,meanT*exp(-mean_min_AS/meanf)),col=2,lty=2)
abline(h=0)
abline(v=0)

library(graphics)
arrows(mean_min_SB,exp(-mean_min_SB/f),mean_min_AS,exp(-mean_min_AS/f),code=3,cex=0.5,length=0.1)
arrows(mean_min_AS,exp(-mean_min_AS/f),mean_min_AS,T*exp(-mean_min_AS/f),code=3,cex=0.5,length=0.05)

text((mean_min_SB+mean_min_AS)/2,(exp(-mean_min_SB/f)+exp(-mean_min_AS/f))/2,"distance effect",pos=4)
text(mean_min_AS,(exp(-mean_min_AS/f)+T*exp(-mean_min_AS/f))/2,"direct effect",pos=4)

# composition of the spatial component
AS_spat_comp<-(meanT*exp(-ASDmat/meanf))%*%rep(1,ncol(ASDmat))
SB_spat_comp<-(exp(-SBDmat/meanf))%*%rep(1,ncol(SBDmat))
# plot(SB_spat_comp,col=4,ylab="sum of the weight",main="cumulated weight on same block(blue) or surrounding blocks(red)")
# lines(AS_spat_comp,col=2,type="p")
# plot(log(SB_spat_comp/AS_spat_comp),xlab="house number",main="log of the ratio of the participations of neighbourgs in the same block\n and on surrounding blocks in the spatial component")

## area version
plot(SB_spat_comp/(AS_spat_comp+SB_spat_comp),type="n",main="B-same block (blue) and surrounding blocks (red)\n participation to the spatial component",ylab="rate")
polygon(c(1,length(AS_spat_comp),length(AS_spat_comp),1),c(0,0,1,1),col=2,border=NA)
polygon(c(1,1:length(SB_spat_comp),length(SB_spat_comp)),c(0,SB_spat_comp/(SB_spat_comp+AS_spat_comp),0),col=4,border=NA)

# mean_ratio_AS_SB<-mean(SB_spat_comp/AS_spat_comp) # Inf as one or the other can have infinite values 

med_of_ratio_SB_AS<-median(SB_spat_comp/AS_spat_comp)
cat("median of the ratio of the participation of the neighbourg on the same block and neighbourgs of surrounding blocks in the spatial component:",med_of_ratio_SB_AS,"\n");

mean_rate_SB_in_spat_component<-mean(SB_spat_comp/(AS_spat_comp+SB_spat_comp))
med_rate_SB_in_spat_component<-median(SB_spat_comp/(AS_spat_comp+SB_spat_comp))
cat("rate of spatial component explained by neighbourgs on the same block:mean",mean_rate_SB_in_spat_component,"median:",med_rate_SB_in_spat_component,"\n")

dev.print(device=pdf,"total_effect_of_streets.pdf")

#### Comments:
# ## short version
# The sharp decrease of the associativity on a same block with the
# distance described by f=11.93  (blue line, fig XA) associated to the lower
# associativity on different blocks described by T=0.19, induce the spatial
# component to be largely dominated by the weights of the neighbors on
# the same block (fig.XB) with a mean of 91.9% of the spatial component
# explained by the neighbors on the same block.
# 
# ## long version
# The influence of streets if first revealed by the strong departure of
# T from its neutral prior (T=0.19). This means that for a same
# distance, the participation to the spatial component is 5 times lower
# for neighbors of surrounding square blocks than for neighbors of the
# same block (blue line versus red line on fig XA). This important
# effect describe only partially the influence of streets as streets
# also impose an increased distance to the neighbors. As we found the
# slope of the kernel to be quite sharp, this increased distance
# diminish the influence of houses on surrounding blocks.
# 
# For example, the mean distance of the first neighbor on the same block
# is 11 meters (dashed blue line fig. XA) when the distance of the first
# neighbor on a different block is 31 meters (dashed red line).
# Following the "on same block kernel" (blue line) this difference in
# distance induce an almost 5 fold decrease of the spatial weight
# (distance effect on fig.XA), as this first neighbor is also on a
# different block, the direct effect described by T induce an other 5
# fold decrease (direct effect). Together this effects lead the first
# neighbor on the same block to have a weight in the spatial component
# 17 times higher than the weight of the first neighbor on a different
# block (mean of the ratio).
# 
# This combined effects are so strong that, despite a much higher number
# of neighbors in surrounding blocks than in the same block, the spatial
# component is very largely dominated by the weights of neighbors on the
# same block (fig.XB-mean: 91.9%, median: 94.1%).


# in case reloading
library("LaplacesDemon")
source("morans_functions.r")
source("spam_complement.r")
source("prep_sampling.r")
source("DeltaSampling.r")

# block<-2008 
# sel<-which(data$block_num==block)
block_data<-data[sel,]
par(mfrow=c(2,2))
block_visudata<-block_data$pos
block_visudata[block_data$insp==0]<-0.5
plot_reel(block_data$easting,block_data$northing,block_visudata,base=0,top=1)

block_u_pred <- u_pred[sel]
plot_reel(block_data$easting,block_data$northing,block_u_pred,base=0,top=max(block_u_pred))
text(block_data$easting,block_data$northing,labels=round(block_u_pred,2),pos=4)

# for all houses known, remove it, get u predicted and by classes see if good
# doing in on the all field is a bit two much with full evaluation of the parameters but if we only guess the spatial field it may be possible
# it is at least possible only on a block

#### init the kernel_fit_space.r
## general
Ku<-mean(smallsampled$Ku)
Kv<-mean(smallsampled$Kv)
K<-c(Ku,Kv,Kc);
c.val<-est.c.val
f<-meanf
T<-meanT
Q.est<-QfromfT(dist_mat,AS,SB,meanf,meanT);

## before each block
starter<-1
dimension<-length(sel)
w<-est.w[sel]
u<-est.u[sel]
bivect<-est.detection[sel]
zposb<-which(data$status[sel]==1)
znegb<-which(data$status[sel]==0)
zNAb<-which(data$status[sel]==9)
y<-as.integer(rnorm(length(w),mean=w,sd=1)>0)
y[zposb]<-1
Q<-Q.est[sel,sel]
R <- makeR(length(u),Q,K);
cholR <- chol.spam(R);
cholQ<-chol.spam(Q)
c.comp<-drop(c.map[sel,]%*%c.val)
grid.stab<-seq(1,length(w),ceiling(length(w)/5))# values of the field tested for stability, keep 5 values

ItTestNum<- gibbsit(NULL,NminOnly=TRUE);
beginEstimate<-1
AdaptOK<-TRUE

nbsimul<-ItTestNum+beginEstimate
nbtraced<-2*(2+length(grid.stab))+4
spacer<-(2+length(grid.stab))

sampledb<-as.matrix(mat.or.vec(nbsimul+1,nbtraced));
sampledb[1,1]<-mean(u)
sampledb[1,2]<-sd(u)
sampledb[1,3:spacer]<-u[grid.stab]
sampledb[1,spacer+1]<-mean(w)
sampledb[1,spacer+2]<-sd(u)
sampledb[1,(spacer+3):(2*spacer)]<-w[grid.stab]
LLHu<-llh.ugivQ(dimension,u,Q,K[1])
sampledb[1,(2*spacer)+1]<-llh.ugivQ(dimension,u,Q,K[1])
sampledb[1,(2*spacer)+2]<-llh.ygivw(y,w);
sampledb[1,(2*spacer)+3]<-llh.zgivy(y,zposb,znegb,bivect);
sampledb[1,(2*spacer)+4]<-llh.zgivw(w,zposb,znegb,bivect);


sum.u.b<-rep(0,length(u))
sum.v.b<-rep(0,length(u))
sum.w.b<-rep(0,length(u))
sum.y.b<-rep(0,length(u))

Rprof()
source("kernel_fit_space.r")
Rprof(NULL)
# par(mfrow=c(2,1))
# plot_reel(block_data$easting,block_data$northing,block_visudata,base=0,top=1)
u_pred.b<-pnorm(est.u.b,0,1)
plot_reel(data$easting[sel],data$northing[sel],pnorm(c.comp,0,1))
text(block_data$easting,block_data$northing,labels=round(pnorm(c.comp,0,1),2),pos=4)
plot_reel(data$easting[sel],data$northing[sel],u_pred.b)
text(block_data$easting,block_data$northing,labels=round(u_pred.b,2),pos=4)

## general analysis
par(mfrow=c(2,2))
hist(est.w)
hist(est.w[data$status==9],add=T,col=4)

hist(pnorm(est.w))
hist(pnorm(est.w[data$status==9]),add=T,col=4)

# making the same for the omitted data

# making the same for the omitted data, taking into account only the same block


# run the model on the rociado data using previously fitted parameters 
nameSimul<-"extrapol_rocIIHunter_exp_epsKc0.01StreetsTr50_insp0.7_nocof_nofitspat_f9.0T0.3"

source("RanalysisFunctions.R")

## full rociado map
roc<-read.csv("maps_hunter_blocks.csv")
roc<-set_to(roc,init=c("NULL"),final=0)
roc$oanimal<-as.numeric(roc$CO!=0 | roc$OV!=0 | roc$AV!=0 | roc$GA!=0 | roc$otros.animales!=0) # careful, don't use cofactors for now: may not be present for all inspected

## apply kernel on cycloI, without streets, only urban and peri-urban
roc$opened<-roc$IIspray
roc$infested<-roc$IIinfested

sel<-which(!is.na(roc$X))
notsel<-which(! 1:length(roc$X) %in% sel)
if(length(notsel)>0){
	warning("Carefull, need to remove:",roc$unicode[notsel],"no X/Y\n")
}
roc<-roc[sel,]
sel<-which(!is.na(roc$block_num))
notsel<-which(! 1:length(roc$X) %in% sel)
if(length(notsel)>0){
cat("Carefull, need to remove:",roc$unicode[notsel],"no block_num.\n")
}
roc<-roc[sel,]

source("extrapol_field.r")
full.screen<-dev.new
updated<-extrapol.spatautocorel(db=roc[,c("X","Y","opened","infested","block_num")])
p.i<-updated$p.i
save.image("make_prior_map.img");

#####
full.screen()
par(mfrow=c(1,2))
plot_reel(roc$X,roc$Y,roc$infested,base=0,top=1,main="p.i")
plot_reel(roc$X,roc$Y,p.i,base=0,top=1,main="p.i")

## checking
db=roc[,c("X","Y","opened","infested")]
full.screen()
par(mfrow=c(1,3))
sel<-which(db$opened==1 & db$infested==0)
hist(p.i[sel],main="dist p.i non-infested",xlim=c(0,1))
sel<-which(db$opened==1 & db$infested==1)
hist(p.i[sel],main="dist p.i infested",xlim=c(0,1))
sel<-which(db$opened==0)
hist(p.i[sel],main="dist p.i non-open",xlim=c(0,1))

full.screen()
par(mfrow=c(2,2))
plot_reel(db$X,db$Y,p.i,base=0,top=1,main="p.i ")
sel<-which(db$opened==1 & db$infested==0)
plot_reel(db$X[sel],db$Y[sel],p.i[sel],base=0,top=1,main="p.i non-infested")
sel<-which(db$opened==1 & db$infested==1)
plot_reel(db$X[sel],db$Y[sel],p.i[sel],base=0,top=1,main="p.i infested")
sel<-which(db$opened==0)
plot_reel(db$X[sel],db$Y[sel],p.i[sel],base=0,top=1,main="p.i non-opened")

roc$est_p.i_db<-p.i # estimate of the day before
write.csv(roc,file="roc_upto2009_p.i.csv",row.names=FALSE)

## estimate of the day after
# post spray corrections

### probability of being observed positive if was observed positive
sel<-which(roc$Iinfested==1 & roc$IIspray==1)
# basic rate of observed infested at II when inf at I
pioIIiI<-length(which(roc$IIinfested[sel]==1))/length(sel)
# then this is the probability of remaining infested * proba of obs
# the probability of being remaining infested is:
source("parameters_extrapol.r")
piIIiI<-pioIIiI/priorinspquality

### probability of being observed positive if was observed negative
sel<-which(roc$Iinfested==0 & roc$Ispray==1 & roc$IIspray==1)
# basic rate of observed infested at II when not inf at I
poiIInoiI<-sum(roc$IIinfested[sel])/length(sel)
# the probability of being infested when not observed infested at I then is:
piIInoiI<-poiIInoiI/priorinspquality
# but the probabilty of being infested if not infested is poiIInoiI as
# the probability of non-detection is the same for I and II 

### probability of opening the door if not infested~ the overestimation of cycle I on cycle II
### that is 
# > byAll
#   Iinfested Ispray est_p.i_da IIinfested fact_overest
#   1         0      0  469.22743        129     3.637422
#   2         0      1   35.60341         30     1.186780
#   3         1      1   52.00000         52     1.000000
## so quickly done we need to devide our estimage by ??? correct by new value from
# 20120301-195007extrapol_rocI_gau_epsKc0.01NoStreetsTr50_insp0.7_nocof_nofitspat_f15.1637_Kok_fast_final
fact_overest_no<-1

estimatedp.i<-roc$est_p.i_db # estimate of the day before
sel<-which(roc$infested==1)
estimatedp.i[sel]<-pioIIiI
# probability of newly infested 
# probability of infested if not observed infested
sel<-which(roc$infested==0 & roc$opened==1)
# observation* (probability of being infested times remaining infested +prob not being infested and becoming infested)
estimatedp.i[sel]<-(estimatedp.i[sel]*piIIiI+poiIInoiI*(1-estimatedp.i[sel]))*priorinspquality

# probability of infested if not sprayed
sel<-which(roc$opened==0)
estimatedp.i[sel]<-priorinspquality*(estimatedp.i[sel]+poiIInoiI*(1-estimatedp.i[sel]))*fact_overest_no

## get 
roc$est_p.i_da<-estimatedp.i ## Nota this is an estimate of OBSERVED
write.csv(roc,file="roc_upto2009_p.i.csv",row.names=FALSE)

#### checking prediction of rociado II by day after I
## mapping
full.screen()
par(mfrow=c(1,2))
plot_reel(roc$X,roc$Y,roc$est_p.i_da,base=0,top=1,main="p.i day after I")
plot_reel(roc$X,roc$Y,roc$IIinfested,base=0,top=1,main="observed II")
# zm()

# overall
sel<-which(roc$IIspray==1) # work only on opened houses at II
byAll<-aggregate(roc[sel,c("est_p.i_da","IIinfested")],by=list(roc$Iinfested[sel],roc$Ispray[sel]),sum)
names(byAll)[1:2]<-c("Iinfested","Ispray")
byAll$fact_overest<-byAll$est_p.i_da/byAll$IIinfested
#=> overestimate by 3.68 something

### overall quality of prediction
library(pROC)
dev.new()
par(mfrow=c(1,2))
sel<-which(roc$IIspray==1)
rocNotCorrected<-roc(roc$IIinfested[sel],roc$est_p.i_da[sel],plot=TRUE)
text(0.2,0.2,label=paste("AUC:",rocNotCorrected$auc))

# correcting for not opened being overestimated
est_p.i_da_adjusted<-roc$est_p.i_da
sel<-which(roc$opened==0)
est_p.i_da_adjusted[sel]<-roc$est_p.i_da[sel]/byAll$fact_overest[1]
sel<-which(roc$opened==1& roc$infested==0)
est_p.i_da_adjusted[sel]<-roc$est_p.i_da[sel]/byAll$fact_overest[2]
sel<-which(roc$IIspray==1)
rocCorrected<-roc(roc$IIinfested[sel],est_p.i_da_adjusted[sel],plot=TRUE)
text(0.2,0.2,label=paste("AUC:",rocCorrected$auc))


full.screen()
par(mfrow=c(2,3))
#### quality base prediction
## quality of prediction by localities
sel<-which(roc$IIspray==1) # work only on opened houses at II
byLoc<-aggregate(roc[sel,c("est_p.i_da","IIinfested")],by=list(roc$Iinfested[sel],roc$Ispray[sel],roc$L[sel],roc$D[sel],roc$P[sel]),sum)
names(byLoc)[1:5]<-c("Iinfested","Ispray","L","D","P")

plot(byLoc$est_p.i_da,byLoc$IIinfested,col=byLoc$Iinfested+byLoc$Ispray+2,pch=3,asp=1)
abline(b=1,a=0,col=1,lty=2)

sel<-which(byLoc$Iinfested==1)
lminfI<-lm(byLoc$IIinfested[sel]~byLoc$est_p.i_da[sel])
abline(lminfI,col=4)
sel<-which(byLoc$Iinfested==0 & byLoc$Ispray ==1)
lmnegI<-lm(byLoc$IIinfested[sel]~byLoc$est_p.i_da[sel])
abline(lmnegI,col=3)
sel<-which(byLoc$Ispray ==0)
lmnoI<-lm(byLoc$IIinfested[sel]~byLoc$est_p.i_da[sel])
abline(lmnoI,col=2)

# ## quality of prediction by district
# sel<-which(roc$IIspray==1) # work only on opened houses at II
# byDis<-aggregate(roc[sel,c("est_p.i_da","IIinfested")],by=list(roc$Iinfested[sel],roc$Ispray[sel],roc$D[sel],roc$P[sel]),sum)
# names(byDis)[1:4]<-c("Iinfested","Ispray","D","P")
# plot(byDis$est_p.i_da,byDis$IIinfested,col=byDis$Iinfested+byDis$Ispray+2,pch=3,asp=1)
# abline(b=1,a=0,col=1,lty=2)
# 
# sel<-which(byDis$Iinfested==1)
# lminfI<-lm(byDis$IIinfested[sel]~byDis$est_p.i_da[sel])
# abline(lminfI,col=4)
# sel<-which(byDis$Iinfested==0 & byDis$Ispray ==1)
# lmnegI<-lm(byDis$IIinfested[sel]~byDis$est_p.i_da[sel])
# abline(lmnegI,col=3)
# sel<-which(byDis$Ispray ==0)
# lmnoI<-lm(byDis$IIinfested[sel]~byDis$est_p.i_da[sel])
# abline(lmnoI,col=2)
# 
# ## quality of prediction by region
# sel<-which(roc$IIspray==1) # work only on opened houses at II
# byRegion<-aggregate(roc[sel,c("est_p.i_da","IIinfested")],by=list(roc$Iinfested[sel],roc$Ispray[sel],roc$region[sel]),sum)
# names(byRegion)[1:3]<-c("Iinfested","Ispray","region")
# plot(byRegion$est_p.i_da,byRegion$IIinfested,col=byRegion$Iinfested+byRegion$Ispray+2,asp=1,pch=as.character(byRegion$region))
# abline(b=1,a=0,col=1,lty=2)
# 
# sel<-which(byRegion$Iinfested==1)
# lminfI<-lm(byRegion$IIinfested[sel]~byRegion$est_p.i_da[sel])
# abline(lminfI,col=4)
# sel<-which(byRegion$Iinfested==0 & byRegion$Ispray ==1)
# lmnegI<-lm(byRegion$IIinfested[sel]~byRegion$est_p.i_da[sel])
# abline(lmnegI,col=3)
# sel<-which(byRegion$Ispray ==0)
# lmnoI<-lm(byRegion$IIinfested[sel]~byRegion$est_p.i_da[sel])
# abline(lmnoI,col=2)

#### quality of prediction corrected
## by locality
roc$est_p.i_da_adjusted<-est_p.i_da_adjusted
sel<-which(roc$IIspray==1) # work only on opened houses at II
byLoc<-aggregate(roc[sel,c("est_p.i_da_adjusted","IIinfested")],by=list(roc$Iinfested[sel],roc$Ispray[sel],roc$L[sel],roc$D[sel],roc$P[sel]),sum)
names(byLoc)[1:5]<-c("Iinfested","Ispray","L","D","P")

plot(byLoc$est_p.i_da_adjusted,byLoc$IIinfested,col=byLoc$Iinfested+byLoc$Ispray+2,pch=3,asp=1)
abline(b=1,a=0,col=1,lty=2)

sel<-which(byLoc$Iinfested==1)
lminfI<-lm(byLoc$IIinfested[sel]~byLoc$est_p.i_da_adjusted[sel])
abline(lminfI,col=4)
sel<-which(byLoc$Iinfested==0 & byLoc$Ispray ==1)
lmnegI<-lm(byLoc$IIinfested[sel]~byLoc$est_p.i_da_adjusted[sel])
abline(lmnegI,col=3)
sel<-which(byLoc$Ispray ==0)
lmnoI<-lm(byLoc$IIinfested[sel]~byLoc$est_p.i_da_adjusted[sel])
abline(lmnoI,col=2)

# ## quality of prediction by district
# sel<-which(roc$IIspray==1) # work only on opened houses at II
# byDis<-aggregate(roc[sel,c("est_p.i_da_adjusted","IIinfested")],by=list(roc$Iinfested[sel],roc$Ispray[sel],roc$D[sel],roc$P[sel]),sum)
# names(byDis)[1:4]<-c("Iinfested","Ispray","D","P")
# plot(byDis$est_p.i_da_adjusted,byDis$IIinfested,col=byDis$Iinfested+byDis$Ispray+2,pch=3,asp=1)
# abline(b=1,a=0,col=1,lty=2)
# 
# sel<-which(byDis$Iinfested==1)
# lminfI<-lm(byDis$IIinfested[sel]~byDis$est_p.i_da_adjusted[sel])
# abline(lminfI,col=4)
# sel<-which(byDis$Iinfested==0 & byDis$Ispray ==1)
# lmnegI<-lm(byDis$IIinfested[sel]~byDis$est_p.i_da_adjusted[sel])
# abline(lmnegI,col=3)
# sel<-which(byDis$Ispray ==0)
# lmnoI<-lm(byDis$IIinfested[sel]~byDis$est_p.i_da_adjusted[sel])
# abline(lmnoI,col=2)
# 
# ## quality of prediction by region
# sel<-which(roc$IIspray==1) # work only on opened houses at II
# byRegion<-aggregate(roc[sel,c("est_p.i_da_adjusted","IIinfested")],by=list(roc$Iinfested[sel],roc$Ispray[sel],roc$region[sel]),sum)
# names(byRegion)[1:3]<-c("Iinfested","Ispray","region")
# plot(byRegion$est_p.i_da_adjusted,byRegion$IIinfested,col=byRegion$Iinfested+byRegion$Ispray+2,asp=1,pch=as.character(byRegion$region))
# abline(b=1,a=0,col=1,lty=2)
# 
# sel<-which(byRegion$Iinfested==1)
# lminfI<-lm(byRegion$IIinfested[sel]~byRegion$est_p.i_da_adjusted[sel])
# abline(lminfI,col=4)
# sel<-which(byRegion$Iinfested==0 & byRegion$Ispray ==1)
# lmnegI<-lm(byRegion$IIinfested[sel]~byRegion$est_p.i_da_adjusted[sel])
# abline(lmnegI,col=3)
# sel<-which(byRegion$Ispray ==0)
# lmnoI<-lm(byRegion$IIinfested[sel]~byRegion$est_p.i_da_adjusted[sel])
# abline(lmnoI,col=2)

#=> note that for non opened, 1 is above, 2 and 3 under the prediction: situation is worse behind city doors than peri-urban and rural doors: at equivalent ranking, should priviledge the city

#### Naive prediction post II, using directly I Nota: estimate of real, not observed infestation
est_p.i_daII<-est_p.i_da_adjusted
## if observed infested at II simply apply the clearance rate
sel<-which(roc$IIspray==1 & roc$IIinfested==1)
est_p.i_daII[sel]<- piIIiI

## what is probability of infested if observed non-infested: P(I|ONI)
# I: infested
# ONI: observed non infested
# OI: observed infested
# q = P(OI|I) = 1-P(ONI|I)
# Pi = P(I) probability given by the model that the house is infested
# P(I|ONI)=P(I inter ONI)/P(ONI) = P(I)P(ONI|I)/P(ONI)
# and according to model P(ONI)= Pi*(1-q)+(1-Pi)
# so P(I|ONI)= Pi *(1-q)/(Pi*(1-q)+(1-Pi))
# or P(I|ONI)= 1/(1+(1/Pi-1)/(1-q))
sel<-which(roc$IIspray==1 & roc$IIinfested==0)
est_p.i_daII[sel]<- (1/(1+(1/est_p.i_daII[sel]-1)/(1-priorinspquality))) # real infestation day before II
est_p.i_daII[sel]<- est_p.i_daII[sel]*piIIiI # real infestation day after

## if not opened, P(I) remains the same, no clearance, nothing

## gravity of the situation:
plot_reel(roc$X,roc$Y,est_p.i_daII,base=0,top=max(est_p.i_daII))

hist(est_p.i_daII)
text(0.15,5000,label=paste("Total:",round(sum(est_p.i_daII),1),"houses infested"))
#=> 265.3 houses infested, not bad compared to the naive application of infestation rate 
# at I to closed at I and II:
naiveRemain<-length(which(roc$Ispray==0 & roc$IIspray==0))*length(which(roc$Iinfested==1))/length(which(roc$Ispray==1))
#=> 1145.9

# of this houses almost all is in non sprayed houses at II:
sel<-which(roc$IIspray==0)
fracInNonOpenedII<-sum(est_p.i_daII[sel])/sum(est_p.i_daII)
#=> 98.9 %

# Can we focus on some of the non opened?
sel<-which(roc$IIspray==0)
cumOrdSum<-cumsum(sort(est_p.i_daII[sel],decreasing=TRUE))
plot(1:length(est_p.i_daII[sel]),cumOrdSum)
# it is a bit diffuse actually,

# To catch 95% of the infested we would need to go to
invCumOrdFn<-ecdf(cumOrdSum)
invCumOrdFn(0.95*max(cumOrdSum))*length(cumOrdSum)
#=> 5300 houses, this is a lot
# even to catch 75% of it we need to check:
invCumOrdFn(0.75*max(cumOrdSum))*length(cumOrdSum)
#=> 2609 houses

# so going proactively to these houses is difficult if not done at spraying:
# - we should insist on having houses returned to after spraying in risked area
# this will represent maybe 20% more in cost and much less problems then
# - we need to estimate the rate of new arrivals
# - we have 127 denuncia with observed bugs, this is only half the estimation of +
# so we are likely to catch only half of the foci
# - maybe the other prior map will be kind of different

roc$est_p.i_daII<-est_p.i_daII
write.csv(roc,"roc_p.i_fromIadjustedwithII.csv",row.names=FALSE)

## plot carte de risque da I vs II
plot_reel(roc$X,roc$Y,roc$est_p.i_da,base=0,top=max(roc$est_p.i_da))
sel<-which(roc$IIinfested==1)
lines(roc$X[sel],roc$Y[sel],pch="I",type="p",col=4)
save.image("make_prior_map.img");

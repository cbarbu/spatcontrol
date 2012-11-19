
nameSimul<-"Pau_exp_increasingThreshold_TimeAndEstimates"

thresholdDistances<-c(10,30,50,70,100)
times<-mat.or.vec(5,4)
Tvalues<-mat.or.vec(5,4)
fvalues<-mat.or.vec(5,4)
SBW<-mat.or.vec(5,4)

### setting values
# copy/pasted from R using load("EndVisu.img") and then:
# mainLoopStop-mainLoopStart
# dim(sampled)
# mean(smallsampled$T)
# quantile(smallsampled$T,prob=c(0.5,0.025,0.975))
# mean(smallsampled$f)
# quantile(smallsampled$f,prob=c(0.5,0.025,0.975))

timesNames<-c("user","system","elapsed","iterations")
valuesNames<-c("mean","median","q0.025","q0.975")

## for the exp kernel
times[1,]<- c(40495.22, 2829.21,48332.77,330024)
Tvalues[1,]<-c(34.96355,5.520372,0.137677,267.721728)
fvalues[1,]<-c(19.57052,11.970660,5.342736,87.451790)
SBW[1,]<-c(0.3119354,0.2115791,0.005663486,0.918239661)

times[2,] <- c(42860.377,8421.966,55800.714,211684)
Tvalues[2,]<-c(0.5231111,0.4705824,0.1780234,1.1957756)
fvalues[2,]<-c(10.72561,10.224765,7.651799,16.595552)
SBW[2,]<-c(0.9066174,0.9108886,0.8252585,0.9596655)

times[3,]<- c(101568.36,27234.59,136169.30,325553)
Tvalues[3,]<-c(0.3697076,0.3459865,0.1461466,0.7182645)
fvalues[3,]<-c(9.359771,9.057240,7.182404,13.063058)
SBW[3,]<-c(0.9296719,0.9321378,0.8804378,0.9653708)

times[4,] <- c(123049.64,29250.28,159926.20,217064)
Tvalues[4,]<-c(0.3552028,0.3184059,0.1347143,0.8134450)
fvalues[4,]<-c(8.845266,8.653502,6.890508,11.802760)
SBW[4,]<-c(0.93685,0.9390565,0.8946616,0.9698684)

times[5,] <- c(1297512.3,419406.2,1765792.9,1426562)
Tvalues[5,]<-c(0.3404479,0.3150142,0.1322863,0.7171205)
fvalues[5,]<-c(8.907073,8.760035,7.065517,11.659143)
SBW[5,]<-c(0.9336569,0.936383,0.8860981,0.9660987)
### plotting
dev.new(width=6.965182,height=2.750677)
par(mfrow=c(1,2),mar=c(4,4,1,2))
plot(thresholdDistances,times[,3]/3600,type="l",lty=1,ylab="Time to convergence (h)",xlab="Threshold distance (m)",ylim=c(0,160))
# lines(thresholdDistances,times[,2]/3600,type="l",lty=2)
# lines(thresholdDistances,times[,1]/3600,type="l",lty=1)

plot(thresholdDistances,times[,3]/times[,4],type="l",lty=1,ylab="Time per iteration (s)",xlab="Threshold distance (m)",ylim=c(0,1.8))
if(interactive())
  dev.print(device=pdf,"time.pdf")

dev.new(width=6.965182,height=2*2.750677)
par(mfrow=c(2,2),mar=c(4,4,1,2))
plot(thresholdDistances,log(Tvalues[,4]+1),type="l",lty=2,ylab=bquote(log(lambda)),xlab="Threshold distance (m)",ylim=c(0,5.5))
lines(thresholdDistances,log(Tvalues[,3]+1),type="l",lty=2)
lines(thresholdDistances,log(Tvalues[,2]+1),type="l",lty=1)
lines(thresholdDistances,log(Tvalues[,1]+1),type="l",lty=1,col=2)

plot(thresholdDistances,Tvalues[,4],type="l",lty=2,ylab=bquote(lambda),xlab="Threshold distance (m)",ylim=c(0,2))
lines(thresholdDistances,Tvalues[,3],type="l",lty=2)
lines(thresholdDistances,Tvalues[,2],type="l",lty=1)
lines(thresholdDistances,Tvalues[,1],type="l",lty=1,col=2)

plot(thresholdDistances,fvalues[,4],type="l",lty=2,ylab=bquote(delta),xlab="Threshold distance (m)",ylim=c(0,30))
lines(thresholdDistances,fvalues[,3],type="l",lty=2)
lines(thresholdDistances,fvalues[,2],type="l",lty=1)
lines(thresholdDistances,fvalues[,1],type="l",lty=1,col=2)

plot(thresholdDistances,log(fvalues[,4]+1),type="l",lty=2,ylab=bquote(log(delta)),xlab="Threshold distance (m)",ylim=c(0,4.5))
lines(thresholdDistances,log(fvalues[,3]+1),type="l",lty=2)
lines(thresholdDistances,log(fvalues[,2]+1),type="l",lty=1)
lines(thresholdDistances,log(fvalues[,1]+1),type="l",lty=1,col=2)
if(interactive())
  dev.print(device=pdf,"delta_lambda.pdf")

dev.new(width=6.965182/2,height=2.750677)
par(mfrow=c(1,1),mar=c(4,4,1,2))
plot(thresholdDistances,SBW[,4],type="l",lty=2,ylab="Same Block Index",xlab="Threshold distance (m)",ylim=c(0,1))
lines(thresholdDistances,SBW[,3],type="l",lty=2)
lines(thresholdDistances,SBW[,2],type="l",lty=1)
lines(thresholdDistances,SBW[,1],type="l",lty=1,col=2)
if(interactive())
  dev.print(device=pdf,"same_block_index.pdf")




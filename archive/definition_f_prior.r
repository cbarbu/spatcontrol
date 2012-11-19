# f is much more intuitive if we use its inverse, what we do from now
# the 2-200 range allows to describe a very large range of behaviour
# and limiting us to this avoid to run into unidentifiability
# of f if we loose ourselves further
xabs<-seq(0,200,1)
plot(xabs,exp(-xabs/100),ylim=c(0,1),type="n",xlab="Distance (m)",ylab=expression(omega[ij]))
lines(xabs,exp(-xabs/2),ylim=c(0,1))
text(20,0.1,"f=2")
lines(xabs,exp(-xabs/20),ylim=c(0,1))
text(45,0.4,"f=20")
lines(xabs,exp(-xabs/200),ylim=c(0,1))
text(90,0.8,"f=200")

# that can correspond to some biological thing for the insect
# it therefore make sense to have a prior that describe this
# the log normal can do this kind of flat tail shapes
# having a standard deviation of 1 on the log scale makes sens here
# then to mu parameter, we can utilize that the mode of the log normal 
# distribution follows exp(mu+sigma^2)
# if we want the log to be centered arround 20 (the strength divided approximately by 3 every 20 meters)
# we then use
par(mfrow=c(1,2))
mu_lognorm<- 1+log(20)
xabs<-seq(-1,7,0.1)
plot(xabs,dlnorm(exp(xabs),mu_lognorm,1),xlab="log(f)",ylab="density")
abline(v=log(2))
abline(v=log(20))
abline(v=log(200))

dev.print(device=pdf,"definition_f_prior.pdf")



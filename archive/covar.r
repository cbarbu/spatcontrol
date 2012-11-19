par(mfrow=c(2,5))
plot(sqrt(1/sampled[,5]),sampled[,11],type="l")
abline(a=0,b=-1,col=4)
plot(sqrt(1/sampled[,5]),sqrt(1/sampled[,10]),type="l")
plot(sqrt(1/sampled[1:nbsimul,5]),c.vals[1:nbsimul,1],type="l")
plot(sqrt(1/sampled[1:nbsimul,5]),c.vals[1:nbsimul,2],type="l")
plot(sqrt(1/sampled[1:nbsimul,5]),c.vals[1:nbsimul,3],type="l")
plot(sqrt(1/sampled[1:nbsimul,5]),c.vals[1:nbsimul,4],type="l")
plot(sqrt(1/sampled[1:nbsimul,5]),c.vals[1:nbsimul,5],type="l")

plot(sqrt(1/sampled[1:nbsimul,5]),sampled[1:nbsimul,1],type="l")
plot(sqrt(1/sampled[1:nbsimul,5]),sampled[1:nbsimul,3],type="l")


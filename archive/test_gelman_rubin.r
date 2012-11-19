#### import chains
chain1file<-"20110808-191305sub1PauTrue5Cof_NoNAInsp_eps0.01_exp_streets/sampled.txt"
chain2file<-"20110808-191150sub1PauTrue5Cof_NoNAInsp_eps0.01_exp_streets/sampled.txt"
burnIn<-50000
parameters<-c(1,3,5,10,15)

chain1<-scan(chain1file,sep="\t"); 
chain2<-scan(chain2file,sep="\t"); 
chain1<-matrix(chain1,ncol=15,byrow=TRUE)
chain2<-matrix(chain2,ncol=15,byrow=TRUE)

chain1Simple<-chain1[-1,parameters]
chain2Simple<-chain2[-1,parameters]

#### perform the gelman rubin
## coda too complex to get chains in matrices
# library(coda)
# chain1bis<-mcmc(chain1)
# chain2bis<-mcmc(chain2)
# gelman.diag(chain1

## boa perfect for a quick look based on the menu interface
library(boa)
# the gelman and rubin is amazingly good
boa.chain.gandr(list(chain1Simple,chain2Simple),alpha=0.05)

## for final figures, locfit plot nice density curves
library(locfit)
# NB removing the burnIn is important for locfit to do well it's job
#    2) small alpha may improve the precision of the final curve
Tfit<-locfit(~c(chain1Simple[-(1:burnIn),1],chain2Simple[-(1:burnIn),1]),xlim=c(0,15),alpha=1)
plot(Tfit,main="autocorrelation ratio for a same distance\naccross streets versus inside block",xlab="T posterior density")

ffit<-locfit(~c(chain1Simple[-(1:burnIn),2],chain2Simple[-(1:burnIn),2]),alpha=0.5)
plot(ffit,main="kernel's shape parameter",xlab="f posterior density")

# only to show the match between locfit and the histogram:
hist(c(chain1Simple[-(1:burnIn),2],chain2Simple[-(1:burnIn),2]),breaks=100,freq=FALSE)
lines(ffit,main="kernel's shape parameter",xlab="f posterior density")
# so locfit is clearly reliable if we want to compare the outputs of various simulation
# and much easier to read

## we now have to plot the comparisons between the kernels
## using this tools
# T/f
# origin from same block
# kernels shapes according to best T/f


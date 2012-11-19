# houses <-read.csv("rociado_ciclo1_xy.csv")
# houses$treated=houses$T_TRAT_RESID+houses$Rec_TRAT_RESID
# houses$tot_vect = houses$T_peridomicilio_Tcap+houses$T_INTRADOM_T_CAP
# houses$presence <- as.integer(houses$tot_vect>0)
# houses2<-houses[which(houses$treated == 1),]

source("pseudo_data_generation.r")
houses<-data
zpos <- which(data$status==1)
zneg <- which(data$status==0)
zNA <- which(data$status==9)
houses$treated<-0
houses$treated[c(zneg,zpos)]<-1
houses$presence<-0
houses$presence[zpos] <- 1;
houses2<-houses[which(houses$treated == 1),]
houses2$Nombre_brigadista<-houses2$collector

# aggregate by brigadista
brigadista_table<-aggregate(houses2$Nombre_brigadista,list(houses2$Nombre_brigadista),FUN=length)

# # add cuyes
# brig_add<-aggregate(houses2[["CU_ANIM_DOM"]],list(houses2$Nombre_brigadista),FUN=sum)
# if(min(brig_add[[1]] == brigadista_table[[1]])==TRUE){
# 	brigadista_table$CU<-brig_add[[2]]
# }else{
# 	cat("WARNING, CU not added, problem in aggregate\n");
# }

# add positive houses
brig_add<-aggregate(houses2$presence,list(houses2$Nombre_brigadista),FUN=sum)
if(min(brig_add[[1]] == brigadista_table[[1]])==TRUE){
	brigadista_table$Tr<-brig_add[[2]]
}else{
	cat("WARNING, presence not added, problem in aggregate\n");
}

# # let us study positive for cuyes for enough reviewed houses
# brigadista_table$rCU <- brigadista_table$CU/brigadista_table$x
# hist(brigadista_table$rCU[brigadista_table$x>20],breaks=10)

brigadista_table$rTr <- brigadista_table$Tr/brigadista_table$x
hist(brigadista_table$rTr[brigadista_table$x>20],breaks=10)

## before anything: look at the data
par(mfcol=c(2,3))
hist(brigadista_table$rTr,main="Histogram of inspectors detection rate")
plot(brigadista_table$rTr,brigadista_table$x,xlab="detection rate",ylab="number of interventions")
# doesn't looks as bimodal as Paucarpata let us map it
inspector <- matrix(0,nrow(houses2),nrow(brigadista_table))
for(i in 1:nrow(brigadista_table)){
	inspector[,i]<- (houses2$Nombre_brigadista==brigadista_table$Group.1[i])
}
houses2$qualInsp<-as.vector(inspector%*%brigadista_table$rTr);
plot_reel(houses2$easting,houses2$northing,houses2$presence,asmin=0,main="map of bugs presence")
plot_reel(houses2$easting,houses2$northing,houses2$qualInsp,asmin=0,asmax=0.5,main="map of inspectors detection rate")
limit<-0.11
sel<-houses2$qualInsp<limit
plot_reel(houses2$easting[sel],houses2$northing[sel],houses2$qualInsp[sel],asmin=0,asmax=0.5,main="only inspectors under 0.1")
sel<-houses2$qualInsp>limit
plot_reel(houses2$easting[sel],houses2$northing[sel],houses2$qualInsp[sel],asmin=0,asmax=0.5,main="only inspectors over 0.1")
## oups, the map shows kind of a clear cut...
dev.print(device=pdf,"inspectors_quality_mapping.pdf")


## EM norm norm modelling
# one can expect that over 20 a binomial is similar to a normal (but problem on variance so after with binomial would be better)

# for commodity, attach the CU rates column
brigadista_table2 <- brigadista_table[brigadista_table$x>20,]

# vector to be analysed
positive_hits <- brigadista_table2$Tr
# positive_hits <- brigadista_table2$CU
total_hits <- brigadista_table2$x

# init values of parameters
equal_var =1
freq_vect=positive_hits/total_hits
model="norm" # model used: "norm" or "binom"
curalpha <- 0.5
curmu0 <- 0.3 # only to be sure it's under curmu1
curmu1 <- 0.5
cursigsq0 <- sd(freq_vect)/5
cursigsq1 <- cursigsq0

#Expectation function
# compute the part of each observation that is in mode 1 normal version
Estep_norm <- function(y,alpha,mu0,mu1,sigsq0,sigsq1){

  # NB previous version was a loop on i (slower)
  prob0 <- (1-alpha)*dnorm(y,mean=mu0,sd=sqrt(sigsq0)) # probability that observation i pertain to the mode 0
  prob1 <- alpha*dnorm(y,mean=mu1,sd=sqrt(sigsq1))     # probability that observation i pertain to the mode 1
  ind <- prob1/(prob0+prob1)			    # part of i that is of mode 1

  ind 
}
# compute the part of each observation that is in mode 1 binomial version
Estep_binom <- function(y,totals,alpha,p0,p1){

  # NB previous version was a loop on i (slower)
  prob0 <- (1-alpha)*dbinom(y,totals,p0) # probability that observation i pertain to the mode 0
  prob1 <- alpha*dbinom(y,totals,p1)     # probability that observation i pertain to the mode 1
  ind <- prob1/(prob0+prob1)			    # part of i that is of mode 1

  ind 
}

#Maximization function
# use the fraction of each observation i that is of mode 1 to compute the best corresponding parameters
# normal version
Mstep_norm <- function(y,ind,equal_var){
  n <- length(y)
  alpha <- sum(ind)/n # alpha is simply the mean of fraction of each i pertaining to 1
  mu1 <- sum(ind*y)/sum(ind) # mu1 is the mean of the y value of the part of the observation of mode 1
  mu0 <- sum((1-ind)*y)/sum(1-ind) # idem mu0
  if(equal_var == 1){
	  sigsq <- sum(ind*((y-mu1)^2))+sum((1-ind)*((y-mu0)^2))
	  sigsq <- sigsq/n
	  sigsq0 <- sigsq
	  sigsq1 <- sigsq
  }else{
	  sigsq1 <- sum(ind*((y-mu1)^2))/sum(ind) # similarly mean of the sig^2 ponderated by pertenance to mode 1
	  sigsq0 <- sum((1-ind)*((y-mu0)^2))/sum(1-ind) # idem
  }
  c(alpha,mu0,mu1,sigsq0,sigsq1) # return the expected parameters
}
# binomial version
Mstep_binom <- function(y,totals,ind){
  n <- length(y)
  rates<-y/totals
  alpha <- sum(ind)/n # alpha is simply the mean of fraction of each i pertaining to 1
  p1 <- sum(ind*rates)/sum(ind) # p1 is the mean of the y value of the part of the observation of mode 1
  p0 <- sum((1-ind)*rates)/sum(1-ind) # idem p0
  
  c(alpha,p0,p1) # return the expected parameters
}

##observed data loglikelihood function
# only usufull to follow the evolution but not used by EM, EM treat "separetly" each parameter
# normal version
loglik.mix_norm <- function(y,ind,alpha,mu0,mu1,sigsq0,sigsq1){
	# for normal-normal model
  loglik <- sum(log(alpha*dnorm(y,mu1,sqrt(sigsq1))+(1-alpha)*dnorm(y,mu0,sqrt(sigsq0))))
  loglik
}
# binomial version
loglik.mix_binom <- function(y,totals,ind,alpha,p0,p1){
	# for binomial-binomial model
  loglik <- sum(log(alpha*dbinom(y,totals,p1)+(1-alpha)*dbinom(y,totals,p0)))
  loglik
}

# initialization
freq_vect <- positive_hits/total_hits
if(model == "binom"){
	curind <- Estep_binom(positive_hits,total_hits,curalpha,curmu0,curmu1)
	loglik <- loglik.mix_binom(positive_hits,total_hits,curind,curalpha,curmu0,curmu1)
	itermat <- c(curalpha,curmu0,curmu1,loglik)
	names_itermat <- c("alpha","mu0","mu1","loglik")
}else if (model == "norm"){
	curind <- Estep_norm(freq_vect,curalpha,curmu0,curmu1,cursigsq0,cursigsq1)
	loglik <- loglik.mix_norm(freq_vect,curind,curalpha,curmu0,curmu1,cursigsq0,cursigsq1)
	itermat <- c(curalpha,curmu0,curmu1,cursigsq0,cursigsq1,loglik)
	names_itermat <- c("alpha","mu0","mu1","sig0","sig1","loglik")
}else{
	warning("unknown model");
}
diff <- 1
numiters <- 1

#Running EM iterations
while (diff > 0.00001 || numiters <= 100){
	if(model == "binom"){
		curind <- Estep_binom(positive_hits,total_hits,curalpha,curmu0,curmu1)
		curparam <- Mstep_binom(positive_hits,total_hits,curind)
	}else if (model == "norm"){
		curind <- Estep_norm(freq_vect,curalpha,curmu0,curmu1,cursigsq0,cursigsq1)
		curparam <- Mstep_norm(freq_vect,curind,equal_var)
		cursigsq0 <- curparam[4]
		cursigsq1 <- curparam[5]
	}
	curalpha <- curparam[1]
	curmu0 <- curparam[2]
	curmu1 <- curparam[3]
	itermat <- rbind(itermat,c(curparam,loglik))

	if(model == "binom"){
		loglik <- loglik.mix_binom(positive_hits,total_hits,curind,curalpha,curmu0,curmu1)
	}else if (model == "norm"){
		loglik <- loglik.mix_norm(freq_vect,curind,curalpha,curmu0,curmu1,cursigsq0,cursigsq1)
	}
	numiters <- numiters + 1
	diff <- max(abs(itermat[numiters,]-itermat[numiters-1,])) 
	# print (c(numiters,loglik))
}

#Tracking iterations
ndisplay=length(itermat[1,]);
ncolumns=ceiling(ndisplay/2);
par(mfrow=c(2,ncolumns))
for (i in 1:ndisplay){
  plot(1:numiters,itermat[,i],type="l",xlab="iterations",ylab=names_itermat[i])
}
scan()
####### following to be revised and allow to check if mixture modelling is ok for binomial ##############################

# plotting fitted mixture density

finalparam<-itermat[numiters,]
alpha <- finalparam[1]
mu0 <- finalparam[2]
mu1 <- finalparam[3]
if(model == "norm"){
	par(mfrow=c(1,2))
	x <- ppoints(1000)*max(freq_vect)
	sigsq0 <- finalparam[4]
	sigsq1 <- finalparam[5]
	y1 <- (1-alpha)*dnorm(x,mu0,sqrt(sigsq0))
	y2 <- alpha*dnorm(x,mu1,sqrt(sigsq1))
}else{
	par(mfrow=c(1,3))
	# NB: for the binomial model this representation is quite bad as it doesn't reflect
	# the incertainty due to various total_hits
	ref_total_hits=ceiling(mean(total_hits));
	regular_hits <- 1:ceiling(max(freq_vect)*ref_total_hits)
	x<-regular_hits/ref_total_hits
	base_y=dbinom(regular_hits,ref_total_hits,mu0);
	scale_factor=max(histo_obs$density)/((1-alpha)*max(base_y))
	y1 <- (1-alpha)*base_y*scale_factor
	y2 <- alpha*dbinom(regular_hits,ref_total_hits,mu1)*scale_factor
}
histo_obs <- hist(freq_vect,prob=T)
lines(x,y1,col=2)
lines(x,y2,col=3)

#Getting Individual probabilities for each value
if(model == "norm"){
	finalindprops <- Estep_norm(freq_vect,alpha,mu0,mu1,sigsq0,sigsq1)
}else{
	# plotting in two dimension for the binomial model
	plot(positive_hits~total_hits)
	lines(c(0,max(total_hits)),c(0,mu0*max(total_hits)))
	lines(c(0,max(total_hits)),c(0,mu1*max(total_hits)))
	lines(1:max(total_hits),qbinom(0.025,1:max(total_hits),mu0))
	lines(1:max(total_hits),qbinom(0.975,1:max(total_hits),mu0))
	lines(1:max(total_hits),qbinom(0.025,1:max(total_hits),mu1))
	lines(1:max(total_hits),qbinom(0.975,1:max(total_hits),mu1))

	finalindprops <-  Estep_binom(positive_hits,total_hits,alpha,mu0,mu1)
}

hist(finalindprops)

lim_prob_good = 0.95 # threshold of "probability to consider an item to be in the "good part"
sum(finalindprops > lim_prob_good)
greatobs <- brigadista_table2[finalindprops>lim_prob_good,]

print(greatobs[order(freq_vect[finalindprops>lim_prob_good],decreasing=TRUE),])
allobs<-brigadista_table2[finalindprops>lim_prob_good,]

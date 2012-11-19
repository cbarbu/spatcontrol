
Delta_prop<-49.90858
Qprop<-QfromDelta(Delta_prop,dist_mat,AS);
Rprop <- makeR(dimension,Qprop,K);
center <- c(rep(0,dimension),y.r);
cholR<-chol(Rprop,memory=list(nnzcolindices=5e6))
x_prop <- drop(rmvnorm.canonical(n=1, b=center, Q=Rprop, memory=list(nnzcolindices=5e6),Rstruct=cholR)); 
LLHproposal <- LLHDeltaQx(Delta_prop,y.r,x_prop,Qprop,K,muDelta,sdDelta);

# to test the quality of the sampler we need to see if it samples with best LLH's an external test is to chech that a random x is not better than the xprop
# an internal test would be to plot the optimal distribution for each xi and see where we get xi's 
# complex, the ideal would be to sample direclty the best xi possible, with no Ku/Kv variance so that 

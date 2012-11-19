#### import chains
chain1file<-"../outputs/20110808-191305sub1PauTrue5Cof_NoNAInsp_eps0.01_exp_streets/sampled.txt"
chain2file<-"../outputs/20110808-191150sub1PauTrue5Cof_NoNAInsp_eps0.01_exp_streets/sampled.txt"

chain1<-scan(chain1file,sep="\t"); 
chain2<-scan(chain2file,sep="\t"); 
chain1<-matrix(chain1,ncol=15,byrow=TRUE)
chain2<-matrix(chain2,ncol=15,byrow=TRUE)

#### perform the gelman rubin
# library(coda)
# chain1bis<-mcmc(chain1)
# chain2bis<-mcmc(chain2)
# gelman.diag(chain1
libary(boa)
boa.chain.gandr(list(chain1,chain2))


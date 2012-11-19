source("intermediary.r")

starter<-nbsimul;
nbsimul<-200000;
f<-sampled[length(sampled),3]
T<-sampled[length(sampled),1]
Ku<-sampled[length(sampled),5]
Kv<-sampled[length(sampled),10]
Kc<-sampled[length(sampled),12]
c.val<-t(c.vals[length(c.vals),])
c.comp<-c.map%*%c.val
K<-c(Ku,Kv,Kc)
Q<-QfromfT(Dmat,AS,SB,f,T);
cholR<-NULL;
cholQ<-NULL;
source("kernel_sampler.r")

source("visualize_sampler.r")

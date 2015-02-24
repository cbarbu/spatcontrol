INTERMEDIARY<-TRUE
source("parameters_sampler.r")
source("functions_intercept.r")
source("pseudo_data_generation.r")
source("prep_sampling.r")
system("tail -n 1 usamples.txt > usampleslast.txt")
system("tail -n 1 wsamples.txt > wsampleslast.txt")
sampled<-read.table(file="sampled.txt")
nbsimul<-length(sampled[,1])

u<-drop(t(scan(file="usampleslast.txt",sep="\t")))
w<-drop(t(scan(file="wsampleslast.txt",sep="\t")))
y<-drop(rnorm(length(w),w,1));
z<-generate_z(y,bivect);
zpos <- which(z.r==1);
zneg <- which(z.r==0);
zNA <- which(z.r==9);

source("visualize_sampler.r")


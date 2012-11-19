library("testthat")
library("spam")
source("spam_complement.r")
source("R/make.partition.matrices.R")

distances = seq(0,145,10)
data.base<-read.csv("DB_simple_Pau_cyclo1_19Jul2011.csv")
notNA<-which(data.base$status!=9)
x<-data.base$easting[notNA]
y<-data.base$northing[notNA]
z<-data.base$status[notNA]
p<-data.base$block_num[notNA]

dist.mat <-nearest.dist(x=cbind(x,y),y=NULL, method="euclidian", delta=max(distances), upper=NULL); # the matrix of distances         
partition.matrices<-make.partition.matrices(dist.mat,p);
N<-partition.matrices[["N"]]
SG<-partition.matrices[["SG"]]
DG<-partition.matrices[["DG"]]
expect_false(length(which(DG@entries<0.9))>0)
expect_true(length(dist.mat@entries)==length(N@entries))
# expect_true(length(DG@entries)==length(N@entries))
# expect_true(length(N@entries)==length(SG@entries))


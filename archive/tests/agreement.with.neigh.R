source("R/agreement.with.neigh.R")
z<-c(0,1,1,1,0,0,0)
N<-matrix(c(1,1,1,0,0,0,0,
	  1,1,1,0,0,0,0,
	  1,1,1,0,0,0,0,
	  0,0,0,1,1,1,1,
	  0,0,0,1,1,1,1,
	  0,0,0,1,1,1,1,
	  0,0,0,1,1,1,1),
	  byrow=TRUE,ncol=7);
vect.agree<-agreement.with.neigh(z,N)

expect_equal(length(vect.agree),length(z))
expect_equal(vect.agree,c(0,1,1,0,2,2,2))

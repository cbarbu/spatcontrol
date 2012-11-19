# allows to fit only the fields, given the other parameters
# called particularly by predict_quality.r
sum.u.b<-rep(0,length(u))
sum.v.b<-rep(0,length(u))
sum.w.b<-rep(0,length(u))
sum.y.b<-rep(0,length(u))
sum.yp.b<-rep(0,length(u))
sum.w.b.sq<-rep(0,length(u))

num.simul<-starter;
while ((num.simul <= nbsimul || !adaptOK) && num.simul<max.nb.simul) {
	# cat("\n\n",num.simul,", out of:",nbsimul," ");
	
	if(use.yprime){
		yprime<-sampleyprime(dimension,w,bivect,zposb,znegb,zNAb);
		y<-sampley(dimension,w,yprime);
	}else{
		y<-sample_y_direct(w,zposb,znegb,zNAb,bivect);
		yprime <- (y>0);
		# cat("likugivQ after y sampling",llh.ugivQ(dimension,u,Q,K[1]),"\n");
	}

	if(use.cofactors && use.v){
		x <- fastsamplexuv(dimension,cholR,y-c.comp);
		u<-x[1:dimension];
		v<-x[dimension+(1:dimension)]-u;

		w<-x[dimension+(1:dimension)]+c.comp;

	}else if(use.cofactors && !use.v){
		x <- samplexcof(dimension,Q,K,y,Qc,cholR);
		u<-x[1:dimension];
		w<-x[dimension+(1:dimension)];
	}else if(!use.cofactors && use.v){
		x <- fastsamplexuv(dimension,cholR,y);
		u<-x[1:dimension];
		v<-x[dimension+(1:dimension)]-u;
		w<-x[dimension+(1:dimension)]
	}else{
		u<-sample_u(dimension,Q,K,y,cholQ);
		w<-u
		# cat("likugivQ after u sampling",llh.ugivQ(dimension,u,Q,K[1]),"\n");
		##

		# Ku<-sampleKu(dimension,Q,u,Kushape,Kuscale); 
		# K[1]<-Ku;
		# cat(" Ku:",Ku);
	}
	# cat("Ku",Ku,"Kv",Kv,"Kc",Kc,"\n");
	# LLHu<-fast.llh.ugivQ(dimension,u,Q,K[1],cholQ)

	if(beginEstimate<=num.simul){
		sum.u.b<-sum.u.b+u;
		if(use.v){
			sum.v.b<-sum.v.b+v;
		}
		sum.w.b<-sum.w.b+w;
		sum.w.b.sq<-sum.w.b.sq+(w)^2;
		sum.y.b<-sum.y.b+y;
		yp<-as.integer(y>0)
		sum.yp.b<-sum.yp.b+yp;
	}
	if(use.v){
		LLHv<-sum(dnorm(v,0,Kv,log=TRUE));
	}
	LLHyw<-llh.ygivw(y,w);
	LLHy<-llh.zgivy(y,zposb,znegb,bivect);
	LLH<-llh.zgivw(w,zposb,znegb,bivect);
	# cat("LLHyw:",LLHyw,"LLHy:",LLHy,"LLH:",LLH,"mu:",mean(u),"sdu:",sd(u));

	sampledb[num.simul+1,1]<-mean(u)
	sampledb[num.simul+1,2]<-sd(u)
	sampledb[num.simul+1,3:spacer]<-u[grid.stab]
	sampledb[num.simul+1,spacer+1]<-mean(w)
	sampledb[num.simul+1,spacer+2]<-sd(u)
	sampledb[num.simul+1,(spacer+3):(2*spacer)]<-w[grid.stab]
	sampledb[num.simul+1,(2*spacer)+1]<-LLHu
	sampledb[num.simul+1,(2*spacer)+2]<-LLHyw;
	sampledb[num.simul+1,(2*spacer)+3]<-LLHy;
	sampledb[num.simul+1,(2*spacer)+4]<-LLH;

	if(num.simul%%freqsave==0 || num.simul==(nbsimul)){
		## manual stop
		try(source("manual.stop.r"))
	}
	## auto stopping
	if((num.simul==ItTestNum+beginEstimate && num.simul>beginEstimate)){
		cat("\n\n",num.simul,", out of:",nbsimul," ");
		cb<-cb.diag(sampledb[(1+beginEstimate):num.simul,],logfile="convergence_tests.txt");
		cat("cb ok ")
		ItTestNum<-cb$newNbIt;
		if(use.autostop){
			nbsimul<-ItTestNum+beginEstimate;
			if(!cb$ok){
				sampledb<-resized(sampledb,nr=nbsimul+1)
				# c.valsamp<-resized(c.valsamp,nr=nbsimul+1)
			}
		}
		cat("new nbsimul:",nbsimul,"\n");

		cat("cb: num.simul \t begEst \t ok \t newNbIt \t new nbsimul \t nbItEff \t burnIn \t Kthin \t Kind \t min(resG)\n",file="convergence_tests.txt",append=TRUE)
		cat("\t",num.simul,"\t",beginEstimate,"\t",cb$ok,"\t",cb$newNbIt,"\t",nbsimul,"\t",cb$nbItEff,"\t",cb$burnIn,"\t",cb$Kthin,"\t",cb$Kind,"\t",cb$"min(resG)","\n",file="convergence_tests.txt",append=TRUE)
	}
	num.simul<-num.simul+1
}
# if(visu.progression){
# 	printdev(device=png,file="finalmaps.pdf",width=1000,height=500)
# }
nbsimul<-(num.simul-1)
est.u.b<-sum.u.b/(nbsimul-beginEstimate+1);
# dump("est.u.b",file=paste("estimated.txt",sep=""),append=TRUE)
est.y.b<-sum.y.b/(nbsimul-beginEstimate+1);
est.yp.b<-sum.yp.b/(nbsimul-beginEstimate+1);
# dump("est.y.b",file=paste("estimated.txt",sep=""),append=TRUE)
est.w.b<-sum.w.b/(nbsimul-beginEstimate+1);
# dump("est.w.b",file=paste("estimated.txt",sep=""),append=TRUE)
est.w.b.sq<-sum.w.b.sq/(nbsimul-beginEstimate+1)
est.sd.w<-sqrt(est.w.b.sq-(est.w.b)^2)
if(use.v){
est.v.b<-sum.v.b/(nbsimul-beginEstimate+1);
# dump("est.v",file=paste("estimated.txt",sep=""),append=TRUE)
}
# if(use.insp){
# 	est.beta<-sum.beta/(nbsimul-beginEstimate+1);
# 	dump("est.beta",file=paste("estimated.txt",sep=""),append=TRUE)
# }
# if(use.cofactors){
# 	est.c.val<-sum.c.val/(nbsimul-beginEstimate+1);
# 	dump("est.c.val",file=paste("estimated.txt",sep=""),append=TRUE)
# }
# # sampled<-sampled[1:nbsimul,];
# simulname<-paste(name,"_tr",threshold,"_T",T.r,"_f",f.r,"Ku",Ku.r,"_",nbsimul,sep="")
# # Rprof(NULL);
# cat("\n")

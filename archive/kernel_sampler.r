full.screen()

# initialization of record variables
nbploted<-2
sum.u<-rep(0,dimension);
if(use.cofactors){
	cat("\n Using cofactors:",cofs,"\n")
	nbploted<-nbploted+1
	sum.c.val<-0*c.val;
}
if(use.v){
	nbploted<-nbploted+1
	sum.v<-rep(0,dimension);
}
sum.w<-rep(0,dimension);
sum.y<-rep(0,dimension);
sum.yp<-rep(0,dimension)
if(use.insp){
	nbploted<-nbploted+1
	sum.beta<-beta*0;
}
if(visu.progression){
	full.screen()
	par(mfcol=c(nbploted,5))
}
# Rprof();
cat("b:",mean(bivect),"nb zneg",length(zneg))

mainLoopStart<-proc.time()
# main loop
num.simul<-starter;
while (num.simul <= nbsimul || (!adaptOK && final.run)) {
  cat("\n\n",num.simul,", out of:",nbsimul," ");
  # update y (latent infestation variable)
  y<-sample_y_direct(w,zpos,zneg,zNA,bivect);
  yprime <- (y>0);

  # update x (spatial component)
  if(use.cofactors && use.v){ #  cofactors and local error
    if(mh.cof){ # metropolis hasting sampling on cofactors
      if(use.MHx && num.simul>begin.MHx){
	x <-samplexMH(dimension,x,Q,K,y,zpos,zneg,bivect,c.comp,sdx.val,cholQ)
      }else{
	if(fit.spatstruct){
	  x <- samplexuv(dimension,Q,K,y-c.comp,cholR);
	}else{
	  x <- fastsamplexuv(dimension,cholR,y-c.comp);
	}
      }
      u<-x[1:dimension];
      v<-x[dimension+(1:dimension)]-u;
      c.all<-mhsamplec(c.val,c.comp,c.map,sdc.val,Kc,x[dimension+(1:dimension)],y);
      c.val<-c.all[[1]]
      c.comp<-c.all[[2]]
      # Kc <- sampleKg(dimension,Qc,c.comp,Kcshape,Kcscale);
      if(fit.spatstruct){
	if(use.MHK){
	  out<-sampleKuMHunifSigma(Ku,u,logsdKu,Q,LLHu=LLHu,cholQ=cholQ)
	  Ku<-out$Ku
	  LLHu<-out$LLHu
	  Kv<-sampleKvMHunifSigma(Kv,v,logsdKv)$Kv
	  K<-c(Ku,Kv)
	}else{
	  K <- sampleK(dimension,Q,x,K.hyper);
	  if(num.simul>num.simul.fixKu){
	    Ku<-K[[1]]
	  }else{
	    K[[1]]<-Ku
	  }
	  Kv<-K[[2]];
	}
      }
      w<-x[dimension+(1:dimension)]+c.comp;
    }else{
      print("not implemented yet")
      break()
    }

  }else if(use.cofactors && !use.v){
    if(mh.cof){
      u<-sample_u(dimension,Q,K,y-c.comp,cholQ);
      wnoc<-u
      if(fit.spatstruct){
	if(num.simul>num.simul.fixKu){
	  Ku <- sampleKu(dimension,Q,u,Kushape,Kuscale); 
	  K[1]<-Ku;
	}
      }
      c.all<-mhsamplec(c.val,c.comp,c.map,sdc.val,Kc,wnoc,y);
      c.val<-c.all[[1]]
      c.comp<-c.all[[2]]
      w<-wnoc+c.comp
    }else{
      x <- samplexuvcof(dimension,Q,K,y,Qc,cholR);
      u<-x[1:dimension];
      w<-x[dimension+(1:dimension)];
      # Kc <- sampleKg(dimension,Qc,w-u,Kcshape,Kcscale);
      if(num.simul.fixKu<num.simul){
	Ku <- sampleKu(dimension,Q,u,Kushape,Kuscale); 
      }
      Kv<-1
    }
  }else if(!use.cofactors && use.v){
    if(fit.spatstruct){
      x <- samplexuv(dimension,Q,K,y,cholR);
    }else{
      x <- fastsamplexuv(dimension,cholR,y);
    }
    u<-x[1:dimension];
    w<-x[dimension+(1:dimension)]
    if(fit.spatstruct){
      K <- sampleK(dimension,Q,x,K.hyper);
      if(num.simul<num.simul.fixKu){
	K[1]<-Ku;
      }else{
	Ku<-K[1];
      }
    }
  }else{
    u<-sample_u(dimension,Q,K,y,cholQ);
    w<-u
    # cat("likugivQ after u sampling",llh.ugivQ(dimension,u,Q,K[1]),"\n");
    ##

    if(fit.spatstruct){
      if(num.simul.fixKu<num.simul){
	Ku<-sampleKu(dimension,Q,u-mean(u),Kushape,Kuscale); 
	K[1]<-Ku;
      }
    }
  }
  cat("Ku",Ku,"Kv",Kv,"Kc",Kc,"\n");
  # LLHu<-fast.llh.ugivQ(dimension,u,Q,Ku,cholQ=cholQ) # for some mysterious reason, cannot be used here, some update of cholQ is not done correctly
  if(fit.spatstruct){
    LLHu<-llh.ugivQ(dimension,u,Q,Ku,cholQ=cholQ) # for some mysterious reason, cannot be used here, some update of cholQ is not done correctly
  }

  ## update the autocorrelation kernel
  if(fit.spatstruct){
    if(!use.cosamplingfT || !firstAdaptOK){
      if(use.f){
	out<-sample_f(u,Ku,T,logsdfprop,f,mf,sdlf,Q,LLHu,AS,SB,cholQ=cholQ,Dmat=Dmat); 
	f<-out$f;
	Q<-out$Q;
	LLHu<-out$LLHu;
	cholQ<-out$cholQ;
      }

      if(use.streets==TRUE){
	out<-sample_T(u,Ku,f,T,logsdTprop,mT,sdT,Q,LLHu,AS,SB,cholQ=cholQ,Dmat=Dmat);
	T<-out$T;
	Q<-out$Q;
	LLHu<-out$LLHu;
	cholQ<-out$cholQ;
      }
    }else{
      out<-sample_fT(u,K,T,logsdTprop,mT,sdlT,f,logsdfprop,sdCoupleFact,mf,sdlf,Q,LLHu,AS,SB,cholQ=cholQ,Dmat=Dmat);
      f<-out[[1]];
      T<-out[[2]];
      Q<-out[[3]];
      LLHu<-out[[4]];
      cholQ<-out[[5]];
    }
  }

  if(use.insp){
    beta <- samplebeta(zpos,zneg,inspector,yprime,abeta,bbeta);
    bivect <- as.vector(inspector %*% beta);
    cat("beta (",mean(beta),",",sd(beta),")");
  }
  if(beginEstimate<=num.simul){
    sum.u<-sum.u+u;
    if(use.v){
      sum.v<-sum.v+v;
    }
    if(use.cofactors){
      sum.c.val<-sum.c.val+c.val;
    }
    sum.w<-sum.w+w;
    sum.y<-sum.y+y;
    if(use.insp){
      sum.beta<-sum.beta+beta;
    }
    yp<-as.integer(y>0)
    sum.yp<-sum.yp+yp;
  }
  if(use.v){
    LLHv<-sum(dnorm(v,0,Kv,log=TRUE));
  }
  if(use.cofactors){
    LLHc<-sum(dnorm(c.val,0,Kc,log=TRUE));
  }
  if(use.insp){
    LLHb<-sum(dbeta(beta,abeta,bbeta,log=TRUE));
  }

  # adapt sampling
  if((num.simul)%%20==0 && !adaptOK){
    adaptOK<-TRUE
    if(fit.spatstruct){
      if(use.cosamplingfT && firstAdaptOK ){
	rateaccept<-mean(tail(acceptfT,20))
	if(rateaccept<lowAcceptRate){
	  sdCoupleFact<-0.9*sdCoupleFact;
	  cat("update of sdCoupleFact to:",sdCoupleFact);
	  adaptOK<-FALSE
	}else if(rateaccept>highAcceptRate){
	  sdCoupleFact<-1.1*sdCoupleFact;
	  cat("update of sdCoupleFact to:",sdCoupleFact);
	  adaptOK<-FALSE
	}
      }else{
	# adapt sampling of Ku
	if(use.MHK){
	  rateaccept<-mean(tail(acceptKu,20))
	  cat("accept rate Ku:",rateaccept);
	  if(rateaccept<lowAcceptRate){
	    logsdKu<-0.9*logsdKu;
	    cat("update of logsdKu to:",logsdKu);
	    adaptOK<-FALSE
	  }else if(rateaccept>highAcceptRate){
	    logsdKu<-1.1*logsdKu;
	    cat("update of logsdKu to:",logsdKu);
	    adaptOK<-FALSE
	  }

	  if(use.v){
	    rateaccept<-mean(tail(acceptKv,20))
	    cat("accept rate Kv:",rateaccept);
	    if(rateaccept<lowAcceptRate){
	      logsdKv<-0.9*logsdKv;
	      cat("update of logsdKv to:",logsdKv);
	      adaptOK<-FALSE
	    }else if(rateaccept>highAcceptRate){
	      logsdKv<-1.1*logsdKv;
	      cat("update of logsdKv to:",logsdKv);
	      adaptOK<-FALSE
	    }
	  }
	}
	# adapt sampling of T
	if(use.streets){
	  rateaccept<-mean(tail(acceptT,20))
	  cat("accept rate T:",rateaccept);
	  if(rateaccept<lowAcceptRate){
	    logsdTprop<-0.9*logsdTprop;
	    cat("update of logsdTprop to:",logsdTprop);
	    adaptOK<-FALSE
	  }else if(rateaccept>highAcceptRate){
	    logsdTprop<-1.1*logsdTprop;
	    cat("update of logsdTprop to:",logsdTprop);
	    adaptOK<-FALSE
	  }
	}
	# adapt sampling of f
	if(use.f){
	  rateaccept<-mean(tail(acceptf,20))
	  cat("accept rate f:",rateaccept);
	  if(rateaccept<lowAcceptRate){
	    logsdfprop<-0.9*logsdfprop;
	    cat("update of logsdfprop to:",logsdfprop);
	    adaptOK<-FALSE
	  }else if(rateaccept>highAcceptRate){
	    logsdfprop<-1.1*logsdfprop;
	    cat("update of logsdfprop to:",logsdfprop);
	    adaptOK<-FALSE
	  }
	}
	if(use.cosamplingfT){
	  if(adaptOK){
	    firstAdaptOK<-TRUE
	    adaptOK<-FALSE
	    cat("firstAdaptOK, num.simul:",num.simul,"\n",file="convergence_tests.txt")
	  }
	}
      }
    }

    if(use.MHx && num.simul<begin.MHx+20){
      adaptOK<-FALSE
    }
    if(use.MHx && num.simul>=begin.MHx+20){
      rateaccept<-mean(tail(acceptx.val,20))
      cat("accept rate x.val:",rateaccept);
      if(rateaccept<lowAcceptRate){
	sdx.val<-0.9*sdx.val;
	cat("update of sdx.val to:",sdx.val);
	adaptOK<-FALSE
      }else if(rateaccept>highAcceptRate){
	sdx.val<-1.1*sdx.val;
	cat("update of sdx.val to:",sdx.val);
	adaptOK<-FALSE
      }
    }
    # adapt sampling of c.val
    if(use.cofactors && mh.cof){
      rateaccept<-mean(tail(acceptc.val,20))
      cat("accept rate c.val:",rateaccept);
      if(rateaccept<lowAcceptRate){
	sdc.val<-0.9*sdc.val;
	cat("update of sdc.val to:",sdc.val);
	adaptOK<-FALSE
      }else if(rateaccept>highAcceptRate){
	sdc.val<-1.1*sdc.val;
	cat("update of sdc.val to:",sdc.val);
	adaptOK<-FALSE
      }
    }
    if(adaptOK){
      cat("\n adaptation of sampling variances OK\n");
      if(final.run){
	beginEstimate<-num.simul+1;
	nbsimul<-beginEstimate+ItTestNum;
	sampled<-resized(sampled,nr=nbsimul+1);
	if(use.cofactors){
	  c.valsamp<-resized(c.valsamp,nr=nbsimul+1)
	}
      }
    }else{
      # adapt sampled size if needed
      if(nbsimul<num.simul+1 && final.run){
	nbsimul<-num.simul+ItTestNum
	sampled<-resized(sampled,nr=nbsimul+1)
	if(use.cofactors){
	  c.valsamp<-resized(c.valsamp,nr=nbsimul+1)
	}
      }
    }
    cat(file="convergence_tests.txt","num.simul:",num.simul,"AdaptOK:",adaptOK,append=TRUE)
    if(use.cosamplingfT){
      cat(file="convergence_tests.txt","firstAdaptOK:",firstAdaptOK,"\n",append=TRUE)
    }
    cat(file="convergence_tests.txt","beginEstimate",beginEstimate,"nbsimul",nbsimul,"\n",append=TRUE)
  }
  LLHyw<-llh.ygivw(y,w);
  LLHy<-llh.zgivy(y,zpos,zneg,bivect);
  LLH<-llh.zgivw(w,zpos,zneg,bivect);
  LLHTotal<-LLH;
  llhf<-lik.f(f,mf,sdlf);
  llhT<-lik.T(T,mT,sdlT);
  llhKu<-dgamma(Ku,shape=Kushape,scale=Kuscale,log=TRUE)
  llhKv<-dgamma(Kv,shape=Kvshape,scale=Kvscale,log=TRUE)

  LLHTotal<-LLH+LLHu+LLHv+LLHc+llhf+llhT+llhKu+llhKv+LLHb


  cat("LLHtotal:",LLHTotal,"(",LLH,"+",LLHu,"+",LLHv,"+",LLHc,"+",llhf,"+",llhT,"+",llhKu,"+",llhKv,"+",LLHb,")\n");
  cat("mu:",mean(u),"sdu:",sd(u),"T:",T,"f:",f);

  if(fit.spatstruct){
    sampled[num.simul+1,1]<-T;
    sampled[num.simul+1,2]<-LLHu;
    sampled[num.simul+1,3]<-f;
    sampled[num.simul+1,4]<-LLHu;
    sampled[num.simul+1,5]<-K[1];
    sampled[num.simul+1,6]<-LLHy;
    sampled[num.simul+1,7]<-LLH;
    sampled[num.simul+1,8]<-LLHyw;
    sampled[num.simul+1,9]<-num.simul;
    sampled[num.simul+1,10]<-K[2];
    sampled[num.simul+1,11]<-mean(u);
    sampled[num.simul+1,12]<-Kc;
    sampled[num.simul+1,13]<-LLHv;
    sampled[num.simul+1,14]<-LLHc;
    sampled[num.simul+1,15]<-LLHb;
    sampled[num.simul+1,16]<-LLHTotal;
  }else{
    sampled[num.simul+1,1]<-mean(u)
    sampled[num.simul+1,2]<-sd(u)
    sampled[num.simul+1,3:spacer]<-u[grid.stab]
    sampled[num.simul+1,spacer+1]<-mean(w)
    sampled[num.simul+1,spacer+2]<-sd(u)
    sampled[num.simul+1,(spacer+3):(2*spacer)]<-w[grid.stab]
    sampled[num.simul+1,(2*spacer)+1]<-LLHu
    sampled[num.simul+1,(2*spacer)+2]<-LLHyw;
    sampled[num.simul+1,(2*spacer)+3]<-LLHy;
    sampled[num.simul+1,(2*spacer)+4]<-LLH;
  }
  if(use.cofactors){
    c.valsamp[num.simul+1,]<-c.val
  }

  # cat("freqsave:",freqsave,"num.simul:",num.simul,"enter stop possibility?", num.simul%%freqsave==0,"\n")
  if(num.simul%%freqsave==0 || num.simul==(nbsimul)){
    write.table(t(u), "usamples.txt", sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
    write.table(t(w), "wsamples.txt", sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
    write.table(t(as.numeric(yprime)), "ypsamples.txt", sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)

    write.table(sampled[(num.simul+1-(freqsave-1)):(num.simul+1),], "sampled.txt",append=TRUE, sep="\t",col.names=FALSE,row.names=FALSE)
    if(visu.progression){
      plot_reel(db$X,db$Y,y,main="y")
      plot_reel(db$X,db$Y,u,main="u")
    }
    if(use.cofactors){
      # plot_reel(db$X,db$Y,w-u,main="cofactors")
      if(visu.progression){
	if(use.generated){
	  plot(c.val~c.val.r)
	}else{
	  plot(c.val)
	}
      }
      write.table(c.valsamp[(num.simul+1-(freqsave-1)):(num.simul+1),], "cofactors.txt", sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
      cat("\n c.val:",c.val,"(",c.val.r,")")
    }
    if(use.v){
      v<-x[dimension+(1:dimension)]-u;
      if(visu.progression){
	plot_reel(db$X,db$Y,v,main="v")
      }
    }
    if(use.insp){
      write.table(t(beta), "betasamples.txt", sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
      if(visu.progression){
	plot_reel(db$X,db$Y,bivect*2-1,main="insp")
      }
    }
    ## manual stop
    try(source("manual.stop.r",local=TRUE))
  }
  ## auto stopping
  if((num.simul==ItTestNum+beginEstimate && num.simul>beginEstimate && final.run)){
    cb<-cb.diag(sampled[(1+beginEstimate):num.simul,-9],logfile="convergence_tests.txt");
    ItTestNum<-min(cb$newNbIt,num.simul*3);
    if(use.autostop){
      nbsimul<-ItTestNum+beginEstimate;
      if(!cb$ok){
	sampled<-resized(sampled,nr=nbsimul+1);
	if(use.cofactors){
	  c.valsamp<-resized(c.valsamp,nr=nbsimul+1)
	}
      }

    }

    cat("cb: num.simul \t begEst \t ok \t newNbIt \t new nbsimul \t nbItEff \t burnIn \t Kthin \t Kind \t min(resG)\n",file="convergence_tests.txt",append=TRUE)
    cat("\t",num.simul,"\t",beginEstimate,"\t",cb$ok,"\t",cb$newNbIt,"\t",nbsimul,"\t",cb$nbItEff,"\t",cb$burnIn,"\t",cb$Kthin,"\t",cb$Kind,"\t",cb$"min(resG)","\n",file="convergence_tests.txt",append=TRUE)
  }
  num.simul<-num.simul+1
}
mainLoopStop<-proc.time()
cat("Time main loop:",mainLoopStop-mainLoopStart,"\n")
save.image(file="EndSampleImage.img")
if(visu.progression){
	printdev(device=png,file="finalmaps.pdf",width=1000,height=500)
}
nbsimul<-(num.simul-1)
est.u<-sum.u/(nbsimul-beginEstimate+1);
dump("est.u",file=paste("estimated.txt",sep=""),append=TRUE)
est.y<-sum.y/(nbsimul-beginEstimate+1);
dump("est.y",file=paste("estimated.txt",sep=""),append=TRUE)
est.w<-sum.w/(nbsimul-beginEstimate+1);
dump("est.w",file=paste("estimated.txt",sep=""),append=TRUE)
est.yp<-sum.yp/(nbsimul-beginEstimate+1);
dump("est.yp",file=paste("estimated.txt",sep=""),append=TRUE)
if(use.v){
	est.v<-sum.v/(nbsimul-beginEstimate+1);
	dump("est.v",file=paste("estimated.txt",sep=""),append=TRUE)
}
if(use.insp){
	est.beta<-sum.beta/(nbsimul-beginEstimate+1);
	dump("est.beta",file=paste("estimated.txt",sep=""),append=TRUE)
}
if(use.cofactors){
	est.c.val<-sum.c.val/(nbsimul-beginEstimate+1);
	dump("est.c.val",file=paste("estimated.txt",sep=""),append=TRUE)
}
sampled<-sampled[1:nbsimul,];
simulname<-paste(name,"_tr",threshold,"_T",T.r,"_f",f.r,"Ku",Ku.r,"_",nbsimul,sep="")
# Rprof(NULL);
cat("\n")

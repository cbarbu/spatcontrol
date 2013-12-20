library(Hmisc); # for errbar
zoom_loaded<-TRUE

ord.zoom<-list()
num.zoom<-0
max.zoom<-0
plot.zoom<-function(zs=1,...){
	out<-try(plot(...))
	if(length(out)==0){ # success of the try
		ord.zoom<<-list()
		num.zoom<<-0
		max.zoom<<-0
		if(zs==1){
			original<-zm()
		}else{
			original<-recordPlot()
		}
		ord.zoom[[1]] <<- original
		num.zoom <<-1# do it after not to change it if failure
		max.zoom <<-1
	}
	cat("num.zoom",num.zoom,"max.zoom",max.zoom,"\n");
	return(original)
}
zplot<-plot.zoom

# # initial prototype
# lines.zoom<-function(...){
# 	out<-try(lines(...))
# 	if(length(out)==0){ # success of the try
# 		recplot<-session.zoom()
# 		ord.zoom[[num.zoom+1]] <<- recplot
# 		num.zoom<<-num.zoom+1 
# 		max.zoom <<-num.zoom
# 	}
# 	cat("num.zoom",num.zoom,"max.zoom",max.zoom,"\n");
# 	return(recplot)
# }
fn.zoom<-function(fn,zs=1,...){
	out<-try(fn(...))
	if(length(out)==0){ # success of the try
		if(zs==1){
			recplot<-session.zoom()
		}else{
			recplot<-recordPlot();
		}
		ord.zoom[[num.zoom+1]] <<- recplot
		num.zoom<<-num.zoom+1 
		max.zoom <<-num.zoom
	}
	cat("num.zoom",num.zoom,"max.zoom",max.zoom,"\n");
	return(recplot)
}
zabline<-function(...){
	fn.zoom(abline,...)
}
zlines<-function(...){
	fn.zoom(lines,...)
}
zplot<-function(...){
	fn.zoom(plot,...)
}
zerrbar<-function(...){
	fn.zoom(errbar,...)
}
# function that plot the recorded plot num_to_plot 
# with the same zoom than the plot num_w_zoom
# or not according to keepzoom (1 yes, other no)
plotsamezoom.zoom<-function(num_to_plot,num_w_zoom,keepzoom){
	rp<-ord.zoom[[num_to_plot]] # get recorded plot
	if(keepzoom==T){ # put zoom of now to recorded plot
		if(! is.null(num_w_zoom)){
			wz<-ord.zoom[[num_w_zoom]]  # get with zoom plot
		}else{
			wz<-recordPlot();
		}
		lim<-getalst(tmp=wz[[1]]);
		zoomplot.zoom(xlim=lim[[1]],ylim=lim[[2]],rp=rp);
	}else{
		replayPlot(rp)
	}
}
back.zoom<-function(keepzoom=F,...){
	if(num.zoom-1<1){
		cat("Sorry, already at first item\n");
	}else{
		plotsamezoom.zoom(num.zoom-1,NULL,keepzoom)
		num.zoom <<- num.zoom-1;
	}
	cat("num.zoom",num.zoom,"max.zoom",max.zoom,"\n");
}
forward.zoom<-function(keepzoom=F,...){
	if(num.zoom+1>max.zoom){
		cat("Sorry, already at final item\n");
	}else{
		plotsamezoom.zoom(num.zoom+1,NULL,keepzoom)
		num.zoom <<- num.zoom+1;
		
	}
	cat("num.zoom",num.zoom,"max.zoom",max.zoom,"\n");
}
bkeepzoom.zoom<-function(...){
	back.zoom(keepzoom=T)
}
fkeepzoom.zoom<-function(...){
	forward.zoom(keepzoom=T)
}

bf.zoom<-function(){
	orig <- recordPlot()
	go_on<-T
	cat("You can go back an forth in saved graphs:\n");
	cat("1: backward\n");
	cat("2: forward\n");
	cat("4: backward same zoom\n");
	cat("5: forward same zoom\n");
	cat("9: original\n");
	while(go_on){
		sel<-scan(n=1)
		if(length(sel)==0){
			go_on<-FALSE;
		}else{
			if(exists("exec.zoom")){
				rm(exec.zoom)
			}
			exec.zoom<-switch(sel,
				back.zoom,
				forward.zoom,
				NULL,
				bkeepzoom.zoom,
				fkeepzoom.zoom,
				NULL,
				NULL,
				NULL,
				orig.zoom
				)
			if(! is.null(exec.zoom)){
				exec.zoom(orig=orig);
			}else{
				cat("Say it again?\n")
			}
		}
	}
	return(recordPlot())
}
bf<-bf.zoom


# # example:
# zplot(runif(1000))
# zlines(c(100,300),c(0.2,0.4))
# zabline(v=600)
#
# # then to go back 
# back.zoom()
# # or to go forward 
# forward.zoom()


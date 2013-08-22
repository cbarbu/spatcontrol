# allow to to zoom on last existing plot using zm()
# optionally, graphs can be plotted with zplot and lines with zlines
# allowing to go back and forth with bf()

# to use it: copy it in your folder and source it:
# source("zoom.r")

# WARNING: when you zoom on a multiplot window, you can only zoom on the last plot and the other plots will be rescaled
# accordingly, can be very usefull or anoying but I don't plan to change this soon

library(Hmisc); # for errbar

zoom_loaded<-1
multipancPoint<-function(ancien,fact=1,lim,point=NULL){
	# cat("ancien",ancien,"fact",fact,"lim:\n")
	# print(lim)
  	if(is.null(point)) point<-mean(ancien)
	newRange<- (1-fact)*point+fact*ancien

	return(newRange);
}
multipanc<-function(ancien,fact,lim,point=NULL){
	# cat("ancien",ancien,"fact",fact,"lim:\n")
	# print(lim)
	meanRange<-mean(ancien)
	newRange<-ancien+(1/fact-1)*(ancien-meanRange)
	return(newRange);
}
keepanc<-function(ancien,fact,lim,point=NULL){
	# cat("ancien",ancien,"fact",fact,"lim:\n")
	# print(lim)
	return(ancien);
}
usenew<-function(ancien,fact,lim,point=NULL){
	# cat("ancien",ancien,"fact",fact,"lim:\n")
	# print(lim)
	return(lim);
}
# to avoid repeating oneself, specially on something quickly changing
# in case of need for change here, also check "locator" in zoomplot.zoom()
is.plot.window<-function(alst){
  ## for former versions of R
  # tmp2 <- all.equal(.Primitive("plot.window"), fn)
  # tmp2 <- (length(grep("plot.window",deparse(fn)))>0)
  ## beginning 3.0.1 
  tmp2 <- all.equal("C_plot_window",alst[[1]]$name)
  return(tmp2)
}

# get the limits of the plot in argument
getalst<-function(tmp=recordPlot()[[1]]){
	for (i in seq(along = tmp)) {
		fn <- tmp[[i]][[1]]
		alst <- as.list(tmp[[i]][[2]])
		tmp2<- is.plot.window(alst)
		if (is.logical(tmp2) && tmp2) {
			alstwin<-alst
		}
	}
	return(alstwin);
}

zoomplot.zoom <- function (xlim=NULL, ylim = NULL,fact=NULL,rp=NULL,x=NULL,y=NULL,...) 
{
  # cat("using zoomplot.zoom")
	# rp is a recorded plot
	# fact is a factor of magnification/outzoom
	# fact has priority on xlim/ylim
	if(! is.null(fact)&& is.numeric(fact) && length(fact)==1){
		xlimfn <-multipancPoint;
		ylimfn <-multipancPoint;
	}else{
		if(!is.null(xlim)&& is.numeric(xlim) && length(xlim)==2){
			xlimfn <- usenew;
		}else{
			xlimfn <-keepanc;
		}
		if(!is.null(ylim) && is.numeric(ylim) && length(ylim)==2){
			ylimfn <-usenew;
		}else{
			ylimfn <-keepanc;
		}
	}

	if(is.null(rp)){
		tmp <- recordPlot()[[1]]
	}else{
		tmp<-rp[[1]]
	}

	for (i in seq(along = tmp)) {
		# cat("i:",i,"\n")
		fn <- tmp[[i]][[1]]
		alst <- as.list(tmp[[i]][[2]])
		tmp1 <- all.equal("C_locator",alst[[1]]$name)
		if (is.logical(tmp1) && tmp1) {
			next # will not like do.call
		}
		tmp2<- is.plot.window(alst)
		if (is.logical(tmp2) && tmp2) {
			# print(alst)
			# cat("alst orig:",alst[[1]],alst[[2]],"\n")
			alst[[2]] <- xlimfn(alst[[2]],fact,xlim,x)
			alst[[3]] <- ylimfn(alst[[3]],fact,ylim,y)
		}
		plotOk<-try(do.call(fn, alst))
	}
	if(class(plotOk)=="try-error"){
	   box() # finish the graph even if errors in between
	}
}

inout.zoom<-function(...){
	cat("Enter magnification factor: \n")
	f<-scan(n=1)
	zoomplot.zoom(fact=f,...);
	return()
}
in.zoom<-function(...){
  # Ideally later should center arround the point selected
  cat("Left click to zoom in\n")
  cat("Right click for other options\n")

  center<-locator(1)
  if(length(center$x)==1){
    zoomplot.zoom(fact=2,...);
    in.zoom()
  }
  return()
}
out.zoom<-function(...){
  # Ideally later should center arround the point selected
  cat("Left click to zoom out\n")
  cat("Right click for other options\n")
  center<-locator(1)
  if(length(center$x)==1){
    zoomplot.zoom(fact=0.5,...);
    out.zoom()
  }
  return()
}
orig.zoom<-function(orig){
	replayPlot(orig)
	return()
}
labelButton<-function(buttons){
  # label the buttons with easy to remember names
  label<-""
# cat("buttons:")
# print(buttons)
  if(length(buttons)>=2){ # rightbutton or scrolling
    if(buttons[2]==2){ # scroll down
      label<-"scrollDown"
    }else if(buttons[2]==1){ # right button
      label<-"right"
    }
  }else if(buttons==1){ # middle button
      label<-"middle"
  }else if(buttons==0){
      label<-"left"
  }else if(buttons==2){# scroll up
      label<-"scrollUp"
  }
  # cat("mevent:",label,"\n")
  return(label)
}
devset <- function(){
  if (dev.cur() != eventEnv$which) dev.set(eventEnv$which)
}


mouseDownsqOut <- function(buttons, x, y) {
    startx <<- x
    starty <<- y
    devset()
    usr <<- par("usr")
    # cat("buttonPress:",buttons,"\n")
    mevent<-labelButton(buttons)
    if(mevent=="scrollDown"){
      zoomDyn(buttons,x,y)
    }else if(mevent=="scrollUp"){ 
      zoomDyn(buttons,x,y)
    }else if(mevent=="middle"){ 
      # cat("Turn on zoomDyn\n")
      eventEnv$onMouseMove <- zoomDyn
    }else if(mevent=="left"){ 
      # cat("Turn on dragmousemove\n")
      eventEnv$onMouseMove <- dragmousemove
    }else if(mevent=="right"){ 
      # cat("Closing...")
      # return(invisible(1))
    }
    NULL
  }

mouseup <- function(buttons, x, y) {
  # plotMesh(coords,triangles, xlim = xlim, ylim = ylim, xaxs = xaxs, yaxs = yaxs)
  eventEnv$onMouseMove <- NULL
}	
keydown <- function(key) {
    
    NULL
  }


sqvsCallBack<-function(...){
    setGraphicsEventHandlers(prompt = "Zoom in and out, hit Esc to quit",
			     onMouseMove = NULL,
			     onMouseDown = mouseDownsqOut,
			     onMouseUp = mouseup,
			     onKeybd = keydown)
    eventEnv <<- getGraphicsEventEnv()
}
sqvsOut.zoom<-function(...){
	cat("Click left over opposite corners of zoom area.\n");
	cat("Click right for zoom out.\n")
	cat("Hit Esc to exit\n")
	sqvsCallBack()

	g<-getGraphicsEvent()
}
sq.zoom<-function(...){
	# use locator to zoom with the mouse (two left clicks)
	# specially, ... can be used to pass a recorded plot rp
	cat("Click left over opposite corners of zoom area.\n");
	cat("Click right for other options.\n")
	square<-locator(2)
	if(length(square)==2){
	  xmin<-min(square$x)
	  xmax<-max(square$x)
	  ymin<-min(square$y)
	  ymax<-max(square$y)
	  zoomplot.zoom(xlim=c(xmin,xmax),ylim=c(ymin,ymax),...)
	  sq.zoom(...)
	}
}
session.zoom<-function(...){
	orig <- recordPlot()
	go_on<-TRUE
	sq.zoom(...)
	while(go_on){
		cat("Do you want to?\n")
		cat("     zoom in: 1\n")
		cat("    zoom out: 2\n")
		cat(" zoom square: 3\n")
		cat(" set magnif.: 4\n")
		cat("back to init: 9\n")
		cat("        Exit: Enter\n")
		sel<-scan(n=1)
		if(length(sel)==0){
			go_on<-FALSE;
		}else{
			if(exists("exec.zoom")){
				rm(exec.zoom)
			}
			exec.zoom<-switch(sel,
				in.zoom,out.zoom,
				sq.zoom,inout.zoom,
				NULL,NULL,NULL,NULL,
				orig.zoom)
			if(!is.null(exec.zoom)){
				exec.zoom(orig=orig,...);
			}else{
				cat("Say it again?\n")
			}
		}
	}
	return(recordPlot());
}
setCallBack<-function(..., xlim = NULL, ylim = NULL, xaxs = "r", yaxs = "r"){
  # X11(type = "Xlib")

  # plot(rnorm(100),rnorm(100))
  startx <- NULL
  starty <- NULL
  usr <- NULL

  #---------------------
  # general functions
  #---------------------
  devset <- function(){
    if (dev.cur() != eventEnv$which) dev.set(eventEnv$which)
  }
  setNav<-function(){
    setGraphicsEventHandlers(prompt = "Scroll Zoom, click and drag, q to quit",
			     onMouseMove = dragmousemove,
			     onMouseDown = mouseDownNavig,
			     onMouseUp = mouseup,
			     onKeybd = keydown)
    eventEnv <<- getGraphicsEventEnv()
  }

  #--------------------
  # Navigation functions
  #--------------------
  multipanc<-function(ancien,fact,lim){
    # return new range according to ancient range and factor of magnification
    cat("ancien",ancien,"fact",fact,"lim:\n")
    # print(lim)
    meanRange<-mean(ancien)
    newRange<-ancien+(1/fact-1)*(ancien-meanRange)
    return(newRange);
  }

  dragmousedown <- function(buttons, x, y) {
    startx <<- x
    starty <<- y
    devset()
    usr <<- par("usr")
    eventEnv$onMouseMove <- dragmousemove
    NULL
  }

  mouseDownNavig <- function(buttons, x, y) {
    startx <<- x
    starty <<- y
    devset()
    usr <<- par("usr")
    # cat("buttonPress:",buttons,"\n")
    mevent<-labelButton(buttons)
    if(mevent=="scrollDown"){
      zoomDyn(buttons,x,y)
    }else if(mevent=="scrollUp"){ 
      zoomDyn(buttons,x,y)
    }else if(mevent=="middle"){ 
      # cat("Turn on zoomDyn\n")
      eventEnv$onMouseMove <- zoomDyn
    }else if(mevent=="left"){ 
      # cat("Turn on dragmousemove\n")
      eventEnv$onMouseMove <- dragmousemove
    }else if(mevent=="right"){ 
      # cat("Closing...")
      # return(invisible(1))
    }
    NULL
  }

  dragmousemove <- function(buttons, x, y) {
    devset()
    # cat("In dragmousemove\n")
    deltax <- diff(grconvertX(c(startx, x), "ndc", "user"))
    deltay <- diff(grconvertY(c(starty, y), "ndc", "user"))
    xlim<<-usr[1:2]-deltax
    ylim <<-usr[3:4]-deltay
    zoomplot.zoom(xlim=xlim,ylim=ylim,...)
    # plot(..., xlim = xlim, xaxs = "i",
    # 	 ylim = ylim, yaxs = "i",pch=".")
    NULL
  }
  zoomDyn <- function(buttons, x, y) {
    devset()
    usr <<- par("usr")
    mevent<-labelButton(buttons)
    if(mevent=="scrollDown"){
      fact<-1.5
    }else if(mevent=="scrollUp"){
      fact<-0.7
    }else {
      deltay <- diff(grconvertY(c(starty, y), "ndc", "user"))
      fact<-max(min(1+deltay/(usr[2]-usr[1]),10),0.1)
    }
    xuser<-grconvertX(x, "ndc", "user")
    yuser<-grconvertY(y, "ndc", "user")
    xlim<<- (1-fact)*xuser+fact*usr[1:2]
    ylim<<- (1-fact)*yuser+fact*usr[3:4]
    zoomplot.zoom(fact=fact,x=xuser,y=yuser)
    # cat("fact:",fact,"\n")
    # cat(usr[1:2],"->",xlim,"\n")
    # cat(usr[3:4],"->",ylim,"\n")
    # plot(..., xlim = xlim, xaxs = "i",
# 	 ylim = ylim, yaxs = "i",pch=".")
    NULL
  }

  mouseup <- function(buttons, x, y) {
    # plotMesh(coords,triangles, xlim = xlim, ylim = ylim, xaxs = xaxs, yaxs = yaxs)
    eventEnv$onMouseMove <- NULL
    NULL
  }	

  #--------------------
  # Modify mesh functions
  #--------------------
  addTriangle<-function(){
    # get the triangle
    currentTriangle<<-identify(n=3,coords$x,coords$y)
    # check it is a triangle

    # check it is not an existing triangle
  i<-1
  filterPoint<-which(!is.na(triangles))
  while(length(filterPoint)>0 && i<4 ){
	  print(filterPoint)
	  filterPoint<-filterPoint[grep(currentTriangle[[i]],triangles[filterPoint])]
	  i<-i+1
  }

  triangles[[length(triangles)+1]]<<-currentTriangle # add it

  polygon(coords$x[currentTriangle],coords$y[currentTriangle]) # plot it
  }
  mouseDownEdit<-function(buttons,x,y){
    startx <<- x
    starty <<- y
    devset()
    usr <<- par("usr")
    # cat("buttonPress:",buttons,"\n")
    mevent<-labelButton(buttons)
    if(mevent=="scrollDown"){
      # zoomDyn(buttons,x,y)
    }else if(mevent=="scrollUp"){ 
      # zoomDyn(buttons,x,y)
    }else if(mevent=="middle"){ 
      # cat("Turn on zoomDyn\n")
      # eventEnv$onMouseMove <- zoomDyn
    }else if(mevent=="left"){ 
      beginMovePoint(buttons,x,y)
    }else if(mevent=="right"){ 
      # cat("Closing...")
      # return(invisible(1))
    }
    NULL
  }

  beginMovePoint<-function(buttons,x,y){
    devset()
    cat("just before identify currentPoint")
    x <- grconvertX(x, "ndc", "user")
    y <- grconvertY(y, "ndc", "user")
    currentPoint<<-identifyFromPoints(x,y,coords)
    # currentPoint<<-identify(n=1,coords)
    text(coords$x[[currentPoint]],coords$y[[currentPoint]],currentPoint)
    print(currentPoint)
    eventEnv$onMouseMove<<-NULL
  }
  movePoint<-function(buttons,x,y){
  }
  mouseUpEdit<-function(buttons,x,y){
    startx <<- x
    starty <<- y
    devset()
    usr <<- par("usr")
    # cat("buttonPress:",buttons,"\n")
    mevent<-labelButton(buttons)
    if(mevent=="scrollDown"){
      # zoomDyn(buttons,x,y)
    }else if(mevent=="scrollUp"){ 
      # zoomDyn(buttons,x,y)
    }else if(mevent=="middle"){ 
      # cat("Turn on zoomDyn\n")
      # eventEnv$onMouseMove <- zoomDyn
    }else if(mevent=="left"){ 
      endMovePoint(buttons,x,y)
    setNav()
    }else if(mevent=="right"){ 
      # cat("Closing...")
      # return(invisible(1))
    }
    NULL
  }

  endMovePoint <-function(buttons,x,y){
    devset()
    cat("just before changing points")
    cat("coords:",coords$x[[currentPoint]],",",coords$y[[currentPoint]],"\n")
    x <- grconvertX(x, "ndc", "user")
    y <- grconvertY(y, "ndc", "user")
    coords$x[[currentPoint]]<<-x
    coords$y[[currentPoint]]<<-y
    cat("coords:",coords$x[[currentPoint]],",",coords$y[[currentPoint]],"\n")
    # plotMesh(coords,triangles, xlim = xlim, ylim = ylim, xaxs = xaxs, yaxs = yaxs)
    eventEnv$onMouseMove <- NULL
  }

  #===============================
  # Keyboard mode commuter
  #===============================

    keydown <- function(key) {
    if (key == "q") return(invisible(1))
    if (key == "n") {
      cat("navigation mode")
      setNav()
      
    }else if (key == "e") {
      cat("editing mode")
      setGraphicsEventHandlers(prompt = "Editing mode",
			       onMouseMove = NULL,
			       onMouseDown = mouseDownEdit,
			       onMouseUp = mouseUpEdit,
			       onKeybd = keydown)
      eventEnv <<- getGraphicsEventEnv()
    }else if (key == "t") {
      setGraphicsEventHandlers(prompt = "Editing mode",
			       onMouseMove = NULL,
			       onMouseDown = mouseDownNavig,
			       onMouseUp = mouseup,
			       onKeybd = keydown)
      eventEnv <<- getGraphicsEventEnv()

      addTriangle()
    }
    NULL
  }

  setGraphicsEventHandlers(prompt = "Click and drag, hit q to quit",
			   onMouseDown = mouseDownNavig,
			   onMouseUp = mouseup,
			   onMouseMove = NULL,
			   onKeybd = keydown)

  eventEnv <<- getGraphicsEventEnv()
}
navigation.zoom<-function(...){
  cat("Scroll to zoom\nLeft click to move\nq to quit\n")
  g<-getGraphicsEvent()
}
# try to replot the graph into a Xlib device to allow events handling
XlibReplot<-function(rp=NULL){
	if(is.null(rp)){
		rp <- recordPlot()
	}

	X11(type = "Xlib")
	try(replayPlot(rp),silent=TRUE)
	test<-try(setCallBack(),silent=TRUE)
	return(test)
}

# Main function, choosing between navigation or old "session" interactions
zm<-function(type=NULL,rp=NULL){
  if(is.null(type)){
    test<-try(setCallBack(),silent=TRUE)
    if(class(test)=="try-error"){
      cat("Device doesn't support event handling. \nTrying to replot in new interface.\n")
      test<-XlibReplot(rp=rp)
      if(class(test)=="try-error"){
	      dev.off()
	      cat("Fall back to classical interface.\n")
	      cat("On Linux/Mac use X11(type = \"Xlib\") to enable navigation.\n\n")

	      type<-"session"
      }else{
	      cat("Replot sucessful, use X11(type = \"Xlib\") to avoid this step.\n")
	      type<-"navigation"
      }
    }else{
      type<-"navigation"
    }
  }
    if(length(grep(type,"session"))>0){
      session.zoom()
    }else if(length(grep(type,"navigation"))>0){
      navigation.zoom()
    }
}
# need to check possibility to 
# orig <- recordPlot()
# X11(type = "Xlib")
# replayPlot(orig)


# # example:
# X11(type = "Xlib") # to allow google map like navigation
# plot(rnorm(1000),rnorm(1000))
# zm()
# # to force old behavior
# zm(type="s")


# # example old:
# plot(runif(1000))
# plot(runif(1000))
# abline(0,1)
# original<-session.zoom()
# 
# # and to go back to the original
# replayPlot(original);

# in addition, launch any plot with a zoom session
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
			original<-session.zoom()
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


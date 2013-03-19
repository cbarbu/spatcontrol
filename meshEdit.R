# example for getGraphicsEvent() from the help()
# adapted to add zoom capability

graphics.off()

#--------------------
# Generic functions
#--------------------
labelButton<-function(buttons){
  # label the buttons with easy to remember names
  label<-""
  if(length(buttons)==2){ # rightbutton or scrolling
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
  cat("mevent:",label,"\n")
  return(label)
}
plotMesh<-function(vertices=NULL,triangles=NULL, xlim = NULL, ylim = NULL, xaxs = "r", yaxs = "r"){
  plot(vertices$x,vertices$y, xlim = xlim, ylim = ylim, xaxs = xaxs, yaxs = yaxs)
  if(length(triangles)>0){
    for(numtri in 1:length(triangles)){
      polygon(vertices$x[triangles[[numtri]]],vertices$y[triangles[[numtri]]])
    }
  }
}
identifyFromPoints<-function(x,y,coords){
  pos<-which.min((coords$x-x)^2+(coords$y-y)^2)
  return(pos)
}

#======================================
# Call back loop
#======================================
meshCallBack <- function(..., xlim = NULL, ylim = NULL, xaxs = "r", yaxs = "r") {
  plotMesh(coords,triangles, xlim = xlim, ylim = ylim, xaxs = xaxs, yaxs = yaxs)
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
	setGraphicsEventHandlers(prompt = "Click and drag, hit q to quit",
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
    # cat("ancien",ancien,"fact",fact,"lim:\n")
    # print(lim)
    meanRange<-mean(ancien)
    newRange<-ancien+(1/fact-1)*(ancien-meanRange)
    return(newRange);
  }
  mouseDownNavig <- function(buttons, x, y) {
    startx <<- x
    starty <<- y
    devset()
    usr <<- par("usr")
    cat("buttonPress:",buttons,"\n")
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
    plot(..., xlim = xlim, xaxs = "i",
	 ylim = ylim, yaxs = "i")
    NULL
  }
  zoomDyn <- function(buttons, x, y) {
    devset()
    mevent<-labelButton(buttons)
    if(mevent=="scrollDown"){
      fact<-0.7
    }else if(mevent=="scrollUp"){
      fact<-1.5
    }else {
      deltay <- diff(grconvertY(c(starty, y), "ndc", "user"))
      fact<-max(min(1+deltay/(usr[2]-usr[1]),10),0.1)
    }
    cat("fact:",fact,"\n")
    xlim<<-multipanc(usr[1:2],fact,1) # 1 is dummy
    ylim<<-multipanc(usr[3:4],fact,1) # 1 is dummy
    plot(..., xlim = xlim, xaxs = "i",
	 ylim = ylim, yaxs = "i")
    NULL
  }

  mouseup <- function(buttons, x, y) {
    plotMesh(coords,triangles, xlim = xlim, ylim = ylim, xaxs = xaxs, yaxs = yaxs)
    eventEnv$onMouseMove <- NULL
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
    cat("buttonPress:",buttons,"\n")
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
    cat("buttonPress:",buttons,"\n")
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
    plotMesh(coords,triangles, xlim = xlim, ylim = ylim, xaxs = xaxs, yaxs = yaxs)
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
meshEdit<-function(){
  g<-0
  while(length(g)!=1 || g!=1){
    meshCallBack(coords,triangles)
    g<-getGraphicsEvent()
  }
}



# #===============================
# # Main loop
# #===============================
# X11(type = "Xlib")
# savepar <- par(ask = FALSE)
# 
# triangles<-list()
# coords<-as.data.frame(cbind(rnorm(1000), rnorm(1000)))
# names(coords)<-c("x","y")
# # This currently only works on the Windows
# # and X11(type = "Xlib") screen devices...
# meshEdit()
# par(savepar)

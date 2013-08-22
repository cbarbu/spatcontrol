\name{zoomplot.zoom}
\alias{zoomplot.zoom}
\title{
Central low level function of the zoom package.
}
\description{
This function allow to replot the current or a saved plot with specific boundaries, magnification factor and possibly arround a user defined x/y.
}
\usage{
zoomplot.zoom(xlim = NULL, ylim = NULL, fact = NULL, rp = NULL, x = NULL, y = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xlim}{
	A vector with min and max x
}
  \item{ylim}{
	A vector with min and max y
}
  \item{fact}{
	A scalar giving the magnification factor (>1 brings you closer)
}
  \item{rp}{
	A previously recorded plot with recordPlot(). With all the corresponding warnings in ?recordPlot. 
}
  \item{x}{
	x of a fix point when rescaling, by default the center.
}
  \item{y}{
	y of a fix point when rescaling, by default the center.
}
  \item{\dots}{
	Additional parameters not implemented, just in case.
}
}
\details{
This function is not necessarily easy to use by hand. It is designed to work well when called from higher level functions. End user should always use zm(). 
}
\value{
Not guaranted for now.
}
\author{
Corentin M. Barbu
}
\note{
This function is the heart of the zoom package and the one that can be affected by R version changes. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
zm, in.zoom
}
\examples{
plot(rnorm(1000),rnorm(1000))
zoomplot.zoom(fact=2,x=0,y=0)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ zoom }
\keyword{ plot }% __ONLY ONE__ keyword per line

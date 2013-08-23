\name{zoom-package}
\alias{zoom-package}
\alias{zoom}
\docType{package}
\title{
Zoom/navigate any plot directly in R.
}
\description{
The single zm() function allows to launch an interactive viewing session with any plot. You can zoom in and out and move the scope (just as in Google maps). Additionally, if a device contains several plots it will simultanuously navigate all the plots.

It fully works under Windows and Linux. 

On Mac it fully works if launching R from the Xquartz terminal.
It will default to the "session" mode if launched from Rgui. R launched from the system terminal may have a problem with X11 fonts that prevents even the simplest plots to be executed well, it will not get any better with this package.

See help(zm) for more details on how to use the package.
}
\details{
\tabular{ll}{
Package: \tab zoom\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2013-08-22\cr
License: \tab LGPL\cr
}
The only function end users should ever use is zm(). It orchestrate the access to other functions like:
\itemize{
\item{navigation.zoom(): }{launch a mouse interaction}
\item{session.zoom(): }{launch a console menu interaction}
}
Themselves orchestrating lower level functions:
\itemize{
\item{in.zoom(): }{zoom in}
\item{out.zoom(): }{zoom out}
\item{inout.zoom(): }{zoom to a magnification factor}
\item{sq.zoom(): }{zoom on a user defined square}
\item{zoomplot.zoom(): }{the heart function reploting everything as directed by higher level functions.}
}
}
\note{
This package is maintained as part of the wider spatcontrol project on github: https://github.com/cbarbu/spatcontrol

Bug reports/suggestions/patches can be directly submitted in this web interface.
}
\author{
Corentin M. Barbu

Maintainer: Corentin M. Barbu <corentin.barbu@gmail.com>
}

\keyword{ zoom }
\keyword{ plot }
\keyword{ navigate }
\keyword{ navigation }
\examples{
plot(rnorm(1000),rnorm(1000))
zm()
}

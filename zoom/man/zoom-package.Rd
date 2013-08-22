\name{zoom-package}
\alias{zoom-package}
\alias{zoom}
\docType{package}
\title{
Zoom/navigate any plot directly in R.
}
\description{
The single zm() function allows to launch an interactive session with any active plot.
See ?zm for details.
}
\details{
\tabular{ll}{
Package: \tab zoom\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2013-08-22\cr
License: \tab LGPL\cr
}
The only function end users should ever use is zm(). It orchestrate the access to other functions like :
- navigation.zoom(): launch a mouse interaction
- session.zoom(): launch a console menu interaction
Themselves orchestrating lower level functions:
- in.zoom(): zoom in
- out.zoom(): zoom out
- inout.zoom(): zoom to a magnification factor
- sq.zoom(): zoom on a user defined square
- zoomplot.zoom(): the heart function reploting everything as directed by higher level functions.
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

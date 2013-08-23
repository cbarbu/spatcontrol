\name{navigation.zoom}
\alias{navigation.zoom}
\alias{session.zoom}
\title{
Opening of an interactive zoom/navigate session.
}
\description{
To launch an interactive session you should use zm() but if you are sure of your device you can launch directly one of these functions.
}
\usage{
navigation.zoom(...)
session.zoom(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
	Everything that can be accepted by sq.zoom.
}
}
\details{
session.zoom launch an interactive console menu to navigate a plot.

navigation.zoom allows to interactively navigate a plot with the mouse.
}
\value{
Returns the final plot, as saved by recordPlot().
}
\author{
Corentin M. Barbu
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
zm().
}
\examples{
\dontrun{
navigation.zoom()
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ session }
\keyword{ navigation }% __ONLY ONE__ keyword per line

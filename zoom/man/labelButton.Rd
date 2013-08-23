\name{labelButton}
\alias{labelButton}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform the button event code passed by the event handler into a meaningful string.
}
\description{
The event handler code returns a hard to understand code. I made a lot of trials and errors and came up with the corresponding function to triage these codes. This may evolve as I better understand this system.
}
\usage{
labelButton(buttons)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{buttons}{
	button event code as passed to event functions by getGraphicsEvent()
}
}
\details{
Used for all analysis of the events in navigation.zoom
The returned values should be self-explaining.
}
\value{
The returned value is a string that can be:
\item{left}{ left button of the mouse}
\item{right}{ right button of the mouse}
\item{middle}{ middle button or scrolling will *pressed*}
\item{scrollDown}{ scrolling weel turned down}
\item{scrollUp}{ scrolling weel turned up}
}
\author{
Corentin M. Barbu
}
\seealso{
setCallBack
}
\examples{
labelButton(c(0,1)) # should return "right"
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ event }

\name{rsgns.param-class}
\Rdversion{1.1}
\docType{class}
\alias{rsgns.param-class}

\title{Class \code{"rsgns.param"}}
\description{
	defines start-time, stop time, readout interval, and other parmaeters for the reactions.
}
\section{Objects from the Class}{
	defines start-time, stop time, readout interval, and other parmaeters for the reactions.
Objects can be created by calls of the form \code{new("rsgns.param", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{time}:}{Object of class \code{"numeric"}}
    \item{\code{stop_time}:}{Object of class \code{"numeric"} }
    \item{\code{readout_interval}:}{Object of class \code{"numeric"} }
    \item{\code{save_interval}:}{Object of class \code{"numeric"} }
    \item{\code{coop}:}{Object of class \code{"numeric"} }
  }
}
\section{Methods}{
No methods defined with class "rsgns.param" in the signature.
}
\author{
	Shailesh Tripathi	
}
\examples{
rp<-new("rsgns.param",time=0,stop_time=1000,readout_interval=500)
showClass("rsgns.param")
}
\keyword{classes}

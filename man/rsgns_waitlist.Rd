\name{rsgns.waitlist-class}
\Rdversion{1.1}
\docType{class}
\alias{rsgns.waitlist-class}

\title{Class \code{"rsgns.waitlist"}}
\description{
	include waiting list parameters to the reactions.
}
\section{Objects from the Class}{
include waiting list parameters to the reactions.
Objects can be created by calls of the form \code{new("rsgns.waitlist", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{nodes}:}{Object of class \code{"numeric"}, defines which nodes to be released molecules at a certain point of time of reaction. }
    \item{\code{time}:}{Object of class \code{"numeric"} }
    \item{\code{mol}:}{Object of class \code{"numeric"}, number of molecules to be released. }
    \item{\code{type}:}{Object of class \code{"character"}}
  }
}
\section{Methods}{
No methods defined with class "rsgns.waitlist" in the signature.
}
\author{
	Shailesh Tripathi
}


\examples{
showClass("rsgns.waitlist")
}
\keyword{classes}

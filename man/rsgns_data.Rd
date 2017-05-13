\name{rsgns.data-class}
\Rdversion{1.1}
\docType{class}
\alias{rsgns.data-class}

\title{Class \code{"rsgns.data"}}
\description{
	Initialize a reaction object consists network information and other reaction related parameters.
}
\section{Objects from the Class}{
	Initialize a reaction object consists network information and other reaction related parameters.
	Objects can be created by calls of the form \code{new("rsgns.data", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{network}:}{Object of class \code{"igraph"}, describes a network}
    \item{\code{rconst}:}{Object of class \code{"numeric"}, reaction rate constant}
    \item{\code{submol}:}{Object of class \code{"numeric"}, substrate molecule}
    \item{\code{promol}:}{Object of class \code{"numeric"}, product molecules}
    \item{\code{dly.rn}:}{Object of class \code{"numeric"}, delay reactions}
    \item{\code{init_pop}:}{Object of class \code{"numeric"}, initial population }
    \item{\code{pop_op}:}{Object of class \code{"character"}}
    \item{\code{pop}:}{Object of class \code{"numeric"}}
    \item{\code{dly.samp}:}{Object of class \code{"list"}}
    \item{\code{rn.rate.function}:}{Object of class \code{"list"}, reaction rate function}
    \item{\code{inhib}:}{Object of class \code{"character"}, inhibitor}
  }
}
\author{
	Shailesh Tripathi
}
\examples{
g <- barabasi.game(10)
rsg <- new("rsgns.data",network=g, rconst=25, inhib=c(rep("*",vcount(g))), rn.rate.function=list("invhill", c(10,2)))
#showClass("rsgns.data")
}
\keyword{classes}

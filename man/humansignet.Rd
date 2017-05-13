\name{humansignet}
\alias{humansignet}
\docType{data}
\title{
	A curated human signaling network data.
}
\description{
	The subnetwork of curated human signaling network  is provided as an "igraph" class object.
	A user can use different igraph functions to access different features of the network.
}
\usage{data("humansignet")}
\format{	
	igraph object of human signalling network.
}
\details{
	This is a subnetwork of human signalling network which contains 325 vertices and 1194 edges. There respective activation and repression of genes can be accessed using E(humansignet)$op. Which provides a vector of +1 and -1. +1 represents the activation and -1 represents supression. 
Example:

	> E(humansignet)[1]
	
	EIF5B->RPSA # if E(humansignet)[1], 
	
	> E(humansignet)$op[1]
	
	1
	
	#In this example the first edge of humansignet network is between EIF5B gene and RPSA gene and,
	
	# +1 shows that AcrR activates the expression of acrA.

    }
\value{loads an igraph object}
\source{
	http://chenglongli.cn/people/EdwinWang/index.html
}
\references{
	Wang, E.: Cancer Systems Biology. Chapman \& Hall/CRC Mathematical and Computational Biology. CRC
Press,(2010). https://books.google.fi/books?id=sdWLbAPo0dwC
}
\examples{
	data(humansignet)
}
\keyword{datasets}

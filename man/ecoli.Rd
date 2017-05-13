\name{ecoli}
\alias{ecoli}
\docType{data}
\title{
	A transcription regulatory network of Ecoli.
}
\description{
	The subnetwork of ecoli transcription regulatory network is provided as an "igraph" class object.
	A user can use different igraph functions to access different features of ecoli network.
}
\usage{data("ecoli")}
\format{
	igraph object of ecoli network.
}
\details{
	This is a subnetwork of ecoli transcription regulatory network which contains 1692 vertices and 3591 edges. There respective activation and repression of genes can be accessed using E(ecoli)$op. Which provides a vector of +1 and -1. +1 represents the activation and -1 represents supression. 
Example:

	> E(ecoli)[1]
	
	AcrR->acrA # if E(ecoli)[1], 
	
	> E(ecoli)$op[1]
	
	-1
	
	#In this example the first edge of ecoli network is between AcrR gene and acrA gene and, -1 shows
	
	#that AcrR supress the expression of acrA.

    }
\value{loads an igraph object}
\source{
	http://regulondb.ccg.unam.mx/menu/download/datasets/files/network_tf_gene.txt
}
\references{
	Gama-Castro S. et al. (2016). "RegulonDB version 9.0: high-level integration of gene regulation, coexpression, motif clustering and beyond". 
   Nucleic Acids Res. 2016 Jan 4;44(D1):D133-43. doi: 10.1093/nar/gkv1156. Epub 2015 Nov 2. 
}
\examples{
	data(ecoli)
}
\keyword{datasets}

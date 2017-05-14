\name{plot.sgnesR}
\alias{plot.sgnesR}
\title{
	plots expression profiles of different molecules
}
\description{
	it is used for the visualization of  distribution and time series plots for different molecules.
}
\usage{
plot.sgnesR(x, mol = "R", method = c("boxplot", "exp"), gn = NULL,ln=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	is a list object returned by 'rsgns.rn' function
}
  \item{mol}{
	it is a characterter variable which profiles to be plotted as a boxplot.

}
  \item{method}{
	boxplot: it is for visualizing distributions for different molecules
	exp: for the visualization of line plot
}
  \item{gn}{
	it is a integer value to select the vertex and its neighbors to visualize the expression vs time. 
}
 \item{ln}{
	this options is used for visualizing logarithmic values of expression profiles
}
\item{...}{Other plot parameters.}
 
}
\value{plots expression profiles.}
\author{
	Shailesh Tripathi, Frank Emmert-Streib
}


\examples{
	 g <- barabasi.game(8)
     
     
     #Declaring initial concentration of RNA and Protein
     
     V(g)$Ppop <- sample(4,vcount(g), rep=TRUE)
     V(g)$Rpop <- sample(4, vcount(g), rep=TRUE)
     
     # Declaring edges as activation or supression for different nodes randomly
     sm <- sample(c(1,-1), ecount(g), rep=TRUE, p=c(.8,.2))
     E(g)$op <- sm
     
     #Declaring reaction rate vector randomly for following reactions
     # 1) Translation rate 2) RNA degradation rate
     # 3) Protein degradation rate 4) Protein binding rate
     # 5) unbinding rate 6) transcription rate
     rc <- sample(4,6, rep=TRUE)
     
     
     # Declaring inhibitor or catalyst for protein-promoters connection unbinding reactions
     inhb <- c("","","","*")
     
     # Declaring reaction rate function for protein-promoters connection ubinding reactions
     rn1 <- list("invhill", c(25,5), c(0,1))
     rn2 <- list("","")
     rn <- list(rn2,rn2,rn1)
     ##### declaring delay vector ##
     dl1 <- list("gaussian", c(10,5))
     dl2 <- list("gaussian", c(15,5))
     dl3 <- list("gaussian", c(20,5))
     
     dl <- list(dl1,dl2,dl3)
     
     # Declaring input data object
     rsg <- new("rsgns.data",network=g,inhib=inhb, rn.rate.function=rn,pop=0, rconst=rc, dly.samp=dl)
     
     # Declaring initial parameters for the reactions
     rp <- new("rsgns.param", time=0, stop_time=60,readout_interval=.2)
     
     #Call the R function for SGN simulator
     xx <- rsgns.rn(rsg, rp, delay=TRUE)
     plot.sgnesR(xx, mol="R", method="boxplot")	
}

\name{rsgns.rn}
\alias{rsgns.rn}
\alias{rsgns.rn}
\title{ Stochastic Gene Networks Simulator }
\description{
	Provides an R interface of SGNSim, 'Stochastic Gene Networks Simulator' developed by Andre S. Ribeiro1,  and Jason Lloyd-Price.
}
\usage{
rsgns.rn(data, param, sample=1, waitlist=NULL, timeseries = TRUE, 
delay = FALSE, induce = FALSE, indpop = NULL) 
}
\arguments{
  \item{data}{

	it is an "S4" object from 'rsgns.data' class which contains following slots
	'network', 'rconst', 'dly.rn', 'dly.samp', 'rn.rate.function', 'inhib'. 
	Each component of the class assigns attributes to the reaction parameters of a network as reactions.

       'network': is a graph object of 'igraph' class, each edge is described as one reaction.

       'rconst':  is a numeric value or a numreic vector of the size '6',  which assigns reaction rates to 
	1) Gene translation and degradation reactions, 2) Connection binding and unbind reactions, 3) transcription rate 	  

       'dly.rn':  is a numeric vector of size 3, which describes the delay in release in product of each reaction in 
	gene translation reaction as ribosomal delay, protein delay and in transcription reaction as a promoter delay.

       'dly.samp': defines a delay release distribution function for the product(s),
	is a list of list objects of size 3, each component is a list object, 
	which defines delays for ribosomal delay, protein delay and promoter delay
    	for first component of list contains function name and the second component is a parameter 
	vector of corresponding distribution of first component of list.

       'rn.rate.function': defines a reaction rate function for the substrates for different reactions. 
	This is a list list object of size 4, each component is a list object, first component of list contains reaction 
	rate function name and the second component is a parameter
        vector of corresponding function of first component of list.

	'inhib': This is a character or a character vector of "*" of size '6' that defines 
	a substrate is acting like a catalyst in a reaction. This is defined for per-gene translation degradation
	reactions, connection biinding and unbinding reactions and transcriptions.
	
	Three more attributes are assigned to network of 'igraph' class for initial populations of proteins and RNAs as follows:
	We can assign a numerical vector to V(g)$Rpop, and V(g)$Pop as initial population,
	and a vector of -1 and 1 to each edge is assigned as activation and repression to E(g)$op in order to 
	represent that a directed edge,  if a molecule activates or supresses to other molecule to which the direction of an
	edge is pointed to.
	
	
		
	Name of each component of the list is the identifier of the species. 
	This details of each input parameter is described here.
	(http://www.cs.tut.fi/~sanchesr/SGN/manual_sgns.pdf)   
	
	
}
  \item{param}{

	it is an "S4" object from 'rsgns.param' class which contains following slots
	'time', 'stop_time', 'readout_interval', 'save_interval', 'save_index', 'coop'. 
	Each component of the class define initial parameters of the reacton.

	
	'time': set the initial clock time

	'stop_time': tell the simulaton to stop when the time reaches stoptime.

	'readout_interval': set the interval between samples of the system in the output file.

  	'save_interval': set the interval between snapshots of the system. 


}
  \item{waitlist}{
	it can be used when the molecules are placed on a waiting list and the priority
queue of molecules to be released later.	
	it is a rsgns.waitlist class object which has following slots:
	nodes: is a numeric vector for the vertices of an input netwrok to be kept in queue.
	time: is a numeric vector that describes the waiting time of molecules mentioned in the 'nodes' slot. The time should be between start time and stop time of reactions. 
	mole: is a numeric vector that describes the number of molecules to be released at a particular time for the  molecules mentioned in the 'nodes' slot.
	type: is a characte type vector which describes that the waiting molecules is whether "RNA" or "Protein" molecules.	
}
	\item{timeseries}{boolean datatype for returning timeseries points or a sample at the end of reaction.}
	\item{delay}{boolean dataype for icluding delay parameter in the reaction.}
\item{sample}{
	is a numeric value to generate the number of samples from a given set of reactions.	
}

  \item{induce}{
	Boolean variable to include a inducer molecule to start translation reaction with the RNA molecule.
}
\item{indpop}{
	is a numerical value to include a molecule as an inducer in translation and degradation reaction}
\item{...}{\code{...} parameter for other inputs.}
}
\value{
  This function returns a list object which has following headers
  	 matrix where rows represents nodes and columns are samples.
}
\references{
	Andre S. Ribeiro; Jason Lloyd-Price, (2007) SGN Sim, a Stochastic Genetic Networks Simulator, Bioinformatics, 23(6):777-779. doi:10.1093/bioinformatics/btm004.
	\url{http://www.cs.tut.fi/~sanchesr/SGN/SGNSim.html}
}
\author{
	Shailesh Tripathi, Frank Emmert-Streib
}
\examples{
	
########################## Example 1 ##################################
#Initialize a network
g <- barabasi.game(10)


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
rc <- sample(3,6, rep=TRUE)


# Declaring inhibitor or catalyst for RNAs for  Translation rate reactions
inhb <- c("*")

# Declaring reaction rate function for RNAs for translation rate reactions
rn1 <- list("invhill", c(25,5), c(0,1))

rn <- list(rn1)
# Declaring input data object
rsg <- new("rsgns.data",network=g,inhib=inhb, rn.rate.function=rn,pop=0, rconst=rc)

# Declaring initial parameters for the reactions
rp <- new("rsgns.param", time=0, stop_time=1,readout_interval=.2)

#Call the R function for SGN simulator
xx <- rsgns.rn(rsg, rp)



########################## Example 2 ##################################

#Initialize a network
g <- barabasi.game(10)


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
rc <- sample(3,6, rep=TRUE)


# Declaring inhibitor or catalyst for RNAs in gene translation reaction, and connection binding reactions
inhb <- c("*","*","*")

# Declaring reaction rate function for RNAs for translation rate reactions
rn1 <- list("invhill", c(25,5), c(0,1))
rn2 <- list("invhill", c(30,5), c(0,1))
rn <- list(rn1, rn2)
# Declaring input data object
rsg <- new("rsgns.data",network=g,inhib=inhb, rn.rate.function=rn,pop=0, rconst=rc)

# Declaring initial parameters for the reactions
rp <- new("rsgns.param", time=0, stop_time=1,readout_interval=.2)

#Call the R function for SGN simulator
xx <- rsgns.rn(rsg, rp)


########################## Example 3 ##################################
##### Using delay option for RBS delay, Protein Delay, and promotors delay
#Initialize a network
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
dlyrn <- sample(50, 3)


# Declaring input data object
rsg <- new("rsgns.data",network=g,inhib=inhb, rn.rate.function=rn,pop=0, rconst=rc, dly.rn=dlyrn)

# Declaring initial parameters for the reactions
rp <- new("rsgns.param", time=0, stop_time=60,readout_interval=.2)

#Call the R function for SGN simulator
xx <- rsgns.rn(rsg, rp, delay=TRUE)


########################## Example 4 ##################################
##### Using delay option for RBS delay, Protein Delay, and promotors delay, 
#by passing a gaussian function with parameters
#Initialize a network
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


	
	#// Example input file for SGNSim's stochastic simulator
	#// Model of an Exclusive Toggle Switch subject to perturbations
	#time 0;
	#stop_time 1000;
	#readout_interval 1;
	#output_file results_toggleswitch.txt;
	#save_interval 250;
	#save_index 1;
	#save_file save_toggleswitch_%%.sim;
	#lua !{
	#	// Cooperativity of the repression
	#	coop = 2;
	#	-- Point where a species would reach
	#	-- equilibrium with no repression
	#	equil = 100;
	#}!
	#population {
	#	A = 0;
	#	B = 0;
	#}
	#reaction {
	#	// Production
	#	*B(invhill:equil/10,coop) --[equil]--> A;
	#	*A(invhill:equil/10,coop) --[equil]--> B;
	#	// Decay
	#	A --[1]--> ;
	#	B --[1]--> ;
	#}
	#// Perturb the system at t = 500 by releasing lots of As
	#queue [2*equil]A(250);
	#// Perturb again at t = 750 by releasing lots of Bs
	#// We should now see at least one toggling
	#queue [2*equil]B(500);
	#// Perturb again at t = 750 by releasing another lots of As
	#queue [2*equil]A(750);

}

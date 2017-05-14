\name{getreactions}
\alias{getreactions}
\title{
	The function allows to obtain a given network topology as a set of chemical reactions. 
}
\description{
	This functions uses the network topology and other parameters and returns a list of reactions defined as 
'initial population', 'activation', 'binding unbinding', and 'transcription and degradation' reactions of a given gene networks. 
}
\usage{
getreactions(data, delay = FALSE, induce = FALSE, indpop = NULL, formula = FALSE, waitlist = NULL)
}
%- maybe also 'usage' for other objects documented here.
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
	represent that a directed edge,  if a molecule activates or supresses to other molecule to which the direction of an edge is pointed to.
	
	
		
	Name of each component of the list is the identifier of the species. 
	This details of each input parameter is described here.
	(http://www.cs.tut.fi/~sanchesr/SGN/manual_sgns.pdf)   
	
	
}
  \item{delay}{
		Boolean variable to include delay parameter.
}
  \item{induce}{
	Boolean variable to include a inducer molecule to start translation reaction with the RNA molecule.
}
  \item{indpop}{
	A numeric variable to define the initial population of inducer molecule.
}
  \item{formula}{
	Boolean variable for the output reactions to be in a matrix form or in a character type. 
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
}
\value{
	returns a list object of rsgns.reactions class which has following sections:
	population: initial population of RNAs, proteins, promoter sites and promoters.
	trans_degradation: translation and degradation reactions in a matrix form or character form.
	binding_unbinding:Binding-unbinding reactions
	activation: This is a set of reactions of the transcription process
	 of the gene to which all possible combinations of the intermediary products
	of the activators of the genes contributes to the expression of genes
}
\author{
	Shailesh Tripathi, Frank Emmert-Streib.
}


\seealso{
	rsgns.rn
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
     
     # Declaring waiting list parameter
     rw <- new("rsgns.waitlist", nodes=c(1,5,7), time=c(20,30,40), mol=c(50,100,200), type=c("RNA", "Protein","RNA"))		
     #Call the R function for SGN simulator
     xx <- getreactions(rsg, delay=TRUE, waitlist=rw)

     ############### Example 2 ############
     ############## user can input its own set of reactions ###
     ### Here is the example of toggle switch reactions ######
     ### These reactions are shown as follows:
     ### *B(<reaction rate function>) --[<rate constant>]--> A;
     ### *A(<reaction rate function>) --[<rate constant>]--> B;
	 #// Decay
     ###  A --[1]--> ;
     ###  B --[1]--> ;
     ###  With sgnesR package it can be run in a following way ##
    
    # 1 ) initialize a dataframe object
    toggle <- getrndf()
    
    # 2 ) set different molecule property and add to object "toggle"
 	
     setmolprop(dfobj="toggle", 1, "B", 1, "s", "invhill:10,2", "*", 2, pop=0)
     setmolprop(dfobj="toggle", 1, "A", 1, "p")
     setmolprop(dfobj="toggle", 2, "A", 1, "s", "invhill:10,2", "*", 2)
     setmolprop(dfobj="toggle", 2, "B", 1, "p")
     setmolprop(dfobj="toggle", 3, "A", 1, "s", rc=1)
     setmolprop(dfobj="toggle", 4, "B", 1, "s", rc=1)
  
   # 3) Declare a waitlist parameter
	# Declaring waiting list parameter
     rw <- new("rsgns.waitlist", time=c(20,30,40), mol=c(50,100,200), type=c("A", "B","A"))
 			
   # 4) get a list of reactions for the input to rsgns.rn function

	xx <- getreactions(toggle, waitlist=rw)	
  
}

\name{setmolprop}
\alias{setmolprop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Declare reaction molecules and assign reaction parameters to each molecule
}
\description{
	
	If a user has its own set of reactions instead of a network then the user can prepare its own list of reactions and input to rsgns.rn function for simulation. 
This function declares reaction molecules and assign reaction parameters to each molecule
}
\usage{
setmolprop(dfobj = NULL, rnindex = NULL, name = NULL, molcount = NULL, type = NULL, rrfn = NULL, inhib = NULL, rc = NULL, delayfn = NULL, delaytime = NULL, pop = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dfobj}{
	is a character type and a reference to empty dataframe object
}
  \item{rnindex}{
	is a numeric type and describe the index for each reaction.
}
  \item{name}{
	is a character type and a describes molecule name.
}
  \item{molcount}{
	is a numeric value,describe number of molecules either participate or released during the reaction.
}
  \item{type}{
	is a character type, describe that a particular molecule is a substrate("s") or a product("p") for a 
	particular reaction index (rnindex)
}
  \item{rrfn}{
	is a character type assigns reaction rate function along with their parameter
	e.g. suppose for a molecule the reaction rate function is inverse hill and has two parameter
	a and b then the it can be described as follows:
	rrfn ="invhill:a,b"
}
  \item{inhib}{
	is a character type and it is used to describe if a molecule by "*" if  a molecule is a 
	catalyst or an inhibitor. 
}
  \item{rc}{
	is a rate constant for a reaction.
}
  \item{delayfn}{
	  is a character type assigns delay function along with their parameter
        e.g. suppose for a molecule the delay function is a gaussian distribution
	 and has two parameter
        a and b then the it can be described as follows:
        delayfn ="gaussian:a,b"
}
  \item{delaytime}{
	is a numeric type, if a user wants to assign a numeric value for the reaction delay. 
}
  \item{pop}{
	is a numeric type, assigns initial population of a particular molecule.
}
}
\value{
	returns a data frame object. 
}
\author{
	Shailesh Tripathi
}
\examples{
 ### Here is the example of toggle switch reactions ######
     ### These reactions are shown as follows:
     ### *B(<reaction rate function>) --[<rate constant>]--> A;
     ###    *A(<reaction rate function>) --[<rate constant>]--> B;
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

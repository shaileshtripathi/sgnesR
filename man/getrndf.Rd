\name{getrndf}
\alias{getrndf}
\title{
	The function creates empty dataframe object and adds different reaction parameters. To prepare a set of reactions. 
}
\description{
	The function creates empty dataframe object and adds different reaction parameters. To prepare a set of reactions. 
}
\usage{
    getrndf()
}
%- maybe also 'usage' for other objects documented here.
\value{
	returns an empty data frame object.
}
\author{
	Shailesh Tripathi, Frank Emmert-Streib.
}


\seealso{
	rsgns.rn
}
\examples{
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

.showlab <- function(verbose=FALSE){
if(verbose){
cat("rnindex: reaction index for substrate and product \n")
cat("name: Molecule name \n")
cat("molcount: number of molecules participate in the reaction or produced in reaction \n")
cat("type: write as substrate or product \n")
cat("rrfn: is a reaction rate function, should be described as e.g. invhill:2,5 \n")
cat("inhib: describe if molecule is an inhibitor or a catalyst \n")
cat("rc: is a numeric value for describing rateconstant \n")
cat("delayfn: is a delay function, should be described as e.g. gaussian:2,5 or gamma:4,2 \n")
cat("delaytime: is a numeric, if user want to give input a delay value instead of delay function \n")
cat("pop: initial population for input molecules. \n")
}

invisible(c("rnindex","name", "molcount", "type","rrfn","inhib","rc", "delayfn","delaytime", "pop"  )) 
}
#toggle switch example:
#coop = 2;
#equil = 100;
#*B(invhill:equil/10,coop) --[equil]--> A;
#*A(invhill:equil/10,coop) --[equil]--> B;
#// Decay
#A --[1]--> ;
#B --[1]--> ;

#toggle <- getrndf()
#setmolprop(dfobj="toggle", 1, "B", 1, "s", "invhill:10,2", "*", 2)
#setmolprop(dfobj="toggle", 1, "A", 1, "p")
#setmolprop(dfobj="toggle", 2, "A", 1, "s", "invhill:10,2", "*", 2)
#setmolprop(dfobj="toggle", 2, "B", 1, "p")
#setmolprop(dfobj="toggle", 3, "A", 1, "s", rc=1)
#setmolprop(dfobj="toggle", 4, "B", 1, "s", rc=1)

.setrn <- function(df=NULL){
options(scipen=999)
tmp1 <- unique(df[,"rnindex"])
rntmp <- c()
pop=NULL
pop <- apply(df,1, function(x)if(!is.na(x[10]))paste(x[2],"=",x[10],sep=""))
pop <- unique(unlist(pop))
df[which(df[,"molcount"]==1),"molcount"] <- ""
df[!is.na(df[,"rrfn"]),"rrfn"] <- paste("(",df[!is.na(df[,"rrfn"]),"rrfn"],")", sep="")
df[!is.na(df[,"delayfn"]),"delayfn"] <- paste("(",df[!is.na(df[,"delayfn"]),"delayfn"],")", sep="")
df[is.na(df)] <- ""
tmpdf <- apply(df, 1, function(x)paste(x[6],x[3],x[2],x[5],x[8],x[9],sep=""))
for(i in tmp1){
        k1 <- which(i==df[,"rnindex"])
        k2 <- which(df[k1,"type"]=="s")
        k3 <- which(df[k1,"type"]=="p")
        rct <- paste("[",df[k1,"rc"][k2[1]],"]",sep="")
        rntmp <- c(rntmp,paste(paste(tmpdf[k1][k2],collapse="+"), "--", as.character(rct),"-->", paste(tmpdf[k1][k3],collapse="+"),sep=""))
}
list(population=pop, reactions=paste(rntmp,";"))
}




setmolprop <- function(dfobj =NULL,rnindex=NULL,name=NULL,  molcount=NULL, type=NULL, rrfn=NULL, inhib=NULL, 
rc=NULL, delayfn=NULL, delaytime=NULL, pop=NULL){
	tmp <- .isvnull(rnindex,name,  molcount, type, rrfn, inhib,
rc, delayfn, delaytime,pop)
	labx <- .showlab()
	names(tmp) <- labx
	df <- getrndf()
	for(i in labx)
	df[1, i] <- tmp[[i]] 
	if(!is.null(dfobj)){
		 assign(dfobj, rbind(get(dfobj),df), envir = .GlobalEnv)
	}
	else{
	return(df)
	}
}
.isvnull <- function(rnindex=NULL,name=NULL,  molcount=NULL, type=NULL, rrfn=NULL, inhib=NULL,
rc=NULL, delayfn=NULL, delaytime=NULL, pop=NULL)
{
	if(is.null(rnindex)){rnindex <- NA}
	if(is.null(name)){name <- NA}
	if(is.null(molcount)){molcount <- NA}
	if(is.null(type)){type <- NA}
	if(is.null(rrfn)){rrfn <- NA}
	if(is.null(inhib)){inhib <- NA}
	if(is.null(rc)){rc <- NA}
	if(is.null(delayfn)){delayfn <- NA}
	if(is.null(delaytime)){delaytime <- NA}
	if(is.null(pop)){pop <- NA}
	list(rnindex,name,  molcount, type, rrfn, inhib,
rc, delayfn, delaytime, pop)
}

getrndf <- function(){data.frame(rnindex=integer(),
                 name=character(),
                 molcount=integer(),
                 type=character(),
                 rrfn=character(),
		 inhib=character(),
		 rc=numeric(),
		 delayfn=character(),
		delaytime=numeric(),
		pop=numeric(),
                 stringsAsFactors=FALSE)
}
#appendrn <- function(p, name=NULL)
getreactions<- function(data,   delay=FALSE, induce=FALSE, indpop=NULL, formula=FALSE, waitlist=NULL){
	xx= NULL
	if(class(data)=="rsgns.data"){
		xx <- .getreactionsnet(data, delay=delay, induce=induce, formula=formula, waitlist=waitlist)
	}
	if(class(data)=="data.frame"){
		xx <- .setrn(data)
		xx$formula<- TRUE
		xx$waitlist <- .wlcheck(waitlist)
		class(xx) <- "rsgns.reactions"
	}
	xx
}
  .wlcheck <- function(waitlist=NULL){
	  wl <- NULL
	  if(!is.null(waitlist)){
	
    	typ <- waitlist@type
	mol <- waitlist@mol
	time <- waitlist@time
	nd <- waitlist@nodes
	nd <- as.numeric(nd)
	if(length(typ)==1){
		typ <- rep(typ, length(mol))
	}
	if(length(mol)==1){
		mol <- rep(mol, length(mol))
	}
	typm <- c()
	for(i in 1:length(typ)){
	if(typ[i]=="RNA"){
		typm <- c(typm,"R")
	} 
	else if(typ[i]=="Protein"){
		typm <- c(typm,"P")
	}
	else{
		typm <- c(typm, typ[i])
	}
	}
	if(length(time)==1){
		time <- rep(time, length(mol))
	}
	nd <- paste(typm, nd,sep="")
	wl <- cbind(mol,nd,time) 
    }
	wl
}

.getreactionsnet <- function(data,   delay=FALSE, induce=FALSE, indpop=NULL, formula=FALSE, waitlist=NULL){
    wl <- NULL	
    if(class(data)!="rsgns.reactions")
    {
        if(is.null(V(data@network)$name)){
            V(data@network)$name <- paste("N",c(1:vcount(data@network)),sep="")
            
        }
    }
    g <- data@network
    rc <- c(data@rconst)
    
    if(length(rc)<6){
     rc <- rep(rc,6)
    }
   if(is.null(rc)||(length(rc)==0)){
     rc <- c(0.002,.005,.005,.005,.01,.2)
    }
  	
   getbub <- function(rnsx){
    bub <- names(rnsx)
    bind <- which(bub=="binding")
    tmp <-c()
    for(i in 1:length(bind)){
        tmp1 <- rnsx[[bind[i]]]
        tmp2 <- rnsx[[(bind[i]+1)]]
        for(j in 1:length(tmp1)){
            tmp <- c(tmp, tmp1[j], tmp2[j])
        }
    }
    tmp
   }
   
  rb_d <- p_d <- pro_d <-  c()

#qind <- NULL;
#if(class(induce)=="numeric"){
#	if(induce[1]>=rp@stop_time){
#
#		induce[1] <- rp@stop_time*c(.9)
#	}
# 	qind <- paste("queue [", induce[2], "]ind(", induce[1],");", sep="")
#	induce=TRUE
# } 
	  
dly.rn <- data@dly.rn  
dly.samp <- data@dly.samp  
if((class(dly.samp)=="list")&& (length(dly.samp)==0)){
dly.samp <- list(list("",""))
}
rn.rate.fn <- data@rn.rate.function
if((class(rn.rate.fn)=="list")&& (length(rn.rate.fn)==0)){
rn.rate.fn <- list(list("",""))
}
inhib <- data@inhib
if(!is.null(rn.rate.fn)){
	if((class(rn.rate.fn[[1]])=="character")&& (class(rn.rate.fn[[2]])=="numeric")){
		rn.rate.fn <- list(rn.rate.fn)
	}

}

if(delay){
	if(!is.null(dly.rn)){
		if(length(dly.rn)<3){
			for(i in (length(dly.rn)+1):3){
				dly.rn <- c(dly.rn,"")
			}
		}
	}
	if(length(dly.samp)==0){
		dly.samp <- rep(list(list("",""),3))
	}
	if(length(dly.samp)!=0){
		if(length(dly.samp)==2){
			if((class(dly.samp[[1]])=="character")&&(class(dly.samp[[2]])=="numeric")){
				dly.samp <- list(dly.samp) 
			}
		}
		if(length(dly.samp)<3){
			tmpl <- list(list())
			for(i in 1:3){
				if(!is.null(dly.samp[i][[1]])){
					tmpl[[i]] <- dly.samp[[i]]
				}
				else{
					tmpl[[i]] <- list("","")
				}
			}
			dly.samp <- tmpl
		}
	}	
   }
else{
	dly.rn <- NULL
	dly.samp <- NULL
}
    rns <- NULL
    if(class(data)!="rsgns.reactions"){
        rns <- .get_inipop(g,rc,delay=delay,dly.rn=dly.rn, dly.samp=dly.samp[1:3], rn.rate.fn=rn.rate.fn, inhib=inhib, induce=induce, indpop=indpop, formula=formula)
    	rns$formula <- formula
	}
    else{
        rns <- data
    }
    if(!is.null(waitlist)){
    	typ <- waitlist@type
	mol <- waitlist@mol
	time <- waitlist@time
	nd <- waitlist@nodes
	nd <- as.numeric(nd)
	nd <- V(g)$name[nd]
	if(length(typ)==1){
		typ <- rep(typ, length(mol))
	}
	if(length(mol)==1){
		mol <- rep(mol, length(mol))
	}
	
	typm <- c()
	for(i in 1:length(typ)){
	if(typ[i]=="RNA"){
		typm <- c(typm,"R")
	} 
	else if(typ[i]=="Protein"){
		typm <- c(typm,"P")
	}
	else{
		typm <- c(typm, typ[i])
	}
	}
	if(length(time)==1){
		time <- rep(time, length(mol))
	}
	nd <- paste(typm, nd,sep="")
	wl <- cbind(mol,nd,time) 
    }
	rns[["waitlist"]] <- wl	
	rns[["network"]] <- g	
	class(rns) <- "rsgns.reactions"
	rns
}


.pf <- function(str){
    
    if(class(str)!="matrix"){
        str <- t(as.matrix(str))
    }
    k1 <- grep("substrate", colnames(str))
    k2 <- grep("product", colnames(str))
    k3 <- grep("rconstant", colnames(str))
    k4 <- grep("delayrn", colnames(str))
    k5 <- grep("dlysmp", colnames(str))
    k6 <- grep("rratefunction", colnames(str))
    k7 <- grep("inhibitor", colnames(str))
    
    
    
    sub <- str[1,k1]
    prod <- str[1,k2]
    rc <- unlist(str[1,k3])
    drn <- rep("", length(prod))
    if((length(k4)==2)){
        drn[1:2] <- str[1,k4]
    }
    if((length(k4)<2) && (length(k4)!=0)){
        drn[1:length(prod)] <- rep(str[1,k4], length(drn))
    }
    dsmp <- rep("", length(prod))
    if((length(k5)==2)){
        dsmp[1:2] <- str[1,k5]
    }
    if((length(k5)<2) && (length(k5)!=0)){
        dsmp[1:length(prod)] <- rep(str[1,k5], length(dsmp))
    }
    rrf <- rep("", length(sub))
    if(length(k6)!=0){
        rrf[1:length(k6)] <- str[1,k6]
    }
    inh <- rep("", length(sub))
    if(length(k7)!=0){
        inh[1:length(k7)] <- str[1,k7]
    }
    rns1 <- c()
    rns2 <- c()
    for(i in 1:length(k1)){
        if(rrf[i]!=""){
            rrf[i] <- paste("(",rrf[i],")", sep="")
        }
        tmp1 <- c(inh[i],sub[i],rrf[i])
        rns1 <- c(rns1, unlist(tmp1))
        rns1 <- c(rns1,"+")
    }
    for(i in 1:length(k2)){
        if(drn[i]!=""){
            drn[i] <- paste("(",drn[i],")", sep="")
        }
        if(dsmp[i]!=""){
            dsmp[i] <- paste("(",dsmp[i],")", sep="")
        }
        tmp2 <- c(prod[i],drn[i],dsmp[i])
        rns2 <- c(rns2, unlist(tmp2))
        rns2 <- c(rns2,"+")
    }
    rns1 <- rns1[1:(length(rns1)-1)]
    rns2 <- rns2[1:(length(rns2)-1)]
    if((rns2[1]=="" )&&(rns2[2]=="")){rns2 <- ""}
    paste(paste(rns1, collapse=""), "--[",rc,"]-->", paste(rns2, collapse=""),";", collapse="", sep="")

}

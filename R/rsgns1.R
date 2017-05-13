rsgns.rn <- function(data, param,  sample=1,waitlist=NULL,  timeseries=TRUE, delay=FALSE, induce=FALSE, indpop=NULL){
 
 qind <- NULL;
	g <- NULL
   if(class(data)!="rsgns.reactions")
    {
        if(is.null(V(data@network)$name)){
            V(data@network)$name <- paste("N",c(1:vcount(data@network)),sep="")
            
        }
    
    g <- data@network
    rc <- c(data@rconst)
    
    if(length(rc)<6){
     rc <- rep(rc,6)
    }
   if(is.null(rc)||(length(rc)==0)){
     rc <- c(0.002,.005,.005,.005,.01,.2)
    }
  	
   
  rb_d <- p_d <- pro_d <-  c()

	  
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
			#dly.rn <- rep(dly.rn,3)
			for(i in (length(dly.rn)+1):3){
				dly.rn <- c(dly.rn,"")
				#dly.rn <- c(dly.rn)
				#print(dly.rn)
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
			#dly.samp <- rep((dly.samp),3)
		}
	}	
   }
else{
	dly.rn <- NULL
	dly.samp <- NULL
}
}
  if(class(data)=="rsgns.reactions"){
    	g <- data$network 
 }		
	params <- c()
    rnd <- paste( sample( 0:9, 20, replace=TRUE ), collapse="" )
    oup <- file.path(".",paste(".rsgns_output",rnd,sep=""))
    savef <- tempfile()
    savef="tmp2"
    
    
    inp <- file.path(".",paste(".rsgns_input",rnd,sep=""))
    
    
    params <- rbind(params, paste("time", param@time,";"))
    params <- rbind(params, paste("stop_time ", param@stop_time,";", sep=" ") )
    params <- rbind(params, paste("readout_interval ", param@readout_interval,";", sep=" ") )
    params <- rbind(params, paste("output_file ", oup,";", sep="") )
    #params <- rbind(params, paste("save_interval ", param@save_interval,";", sep=" ") )
    #params <- rbind(params, paste("save_index ", param@save_index,";", sep=" ") )
    #params <- rbind(params, paste("save_file ", savef,"_%%.sim;", sep="") )
    params <- rbind(params, paste("lua!{\n coop = 0", param@coop, ";\n}!"))
    #reactions <- getreactions(data)
    #reactions <- rbind(paste("reaction {\n"), reactions)
    #reactions <- rbind(reactions, "\n }")
    #print(dly.samp)
    #print("******")		
    rns <- NULL
    if(class(data)!="rsgns.reactions"){
        #rns <- getreactions(data, waitlist=waitlist, delay=delay, induce=induce,indpop=indpop)
        rns <- .get_inipop(g,rc,delay=delay,dly.rn=dly.rn, dly.samp=dly.samp[1:3], rn.rate.fn=rn.rate.fn, inhib=inhib, induce=induce, indpop=indpop, formula=TRUE)
	#print(rns)
	wl <- NULL
    if(!is.null(waitlist)){
    	typ <- waitlist@type
	mol <- waitlist@mol
	time <- waitlist@time
	nd <- waitlist@nodes
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
	if(typ[i]=="Protein"){
		typm <- c(typm,"P")
	}
	}
	if(length(time)==1){
		time <- rep(time, length(mol))
	}
	nd <- paste(typm, nd,sep="")
	wl <- cbind(mol,nd,time)
	wl  <- apply(wl, 1, function(x)paste("queue", paste("[",x[1],"]",x[2],"(",x[3],");", sep="") )) 
    }
	rns[["waitlist"]] <- wl		
	#print(class(rns))
	#data <- rns
    }
    if(class(data)=="rsgns.reactions"){
	wl <- data$waitlist 
	if(!is.null(wl)){
		wl  <- apply(wl, 1, function(x)paste("queue", paste("[",x[1],"]",x[2],"(",x[3],");", sep="") )) 
 	}
	if(data$formula){
	data[["waitlist"]] <- wl
	rns <- data	
	}	
	else{
	population <- apply(data[[1]],1,function(x)paste(x, collapse=" = "))
	trans_degradation <- sapply(data[[2]], function(x){if(!is.null(x)){apply(x,1,.pf)}})	
	binding_unbinding <- sapply(data[[3]], function(x){if(!is.null(x)){apply(x,1,.pf)}})	
	activation <- unlist(sapply(data[[4]], function(x){if(!is.null(x)){apply(x,1,.pf)}}))
	
	rns <- list(population=population, trans_degradation=trans_degradation,binding_unbinding=binding_unbinding,activation=activation,waitlist=wl)
	}
	}
    pop <- unlist(rns[["population"]])
    
    names(pop)=c()
    pop <- paste(pop,";")
    #print(names(rns))
    #bub <-   names(rns$reactions$binding_unbinding)
    #print(bub)	
    #binding <- rns$reactions$binding_unbinding[which(bub=="binding")]
    #unbinding <- rns$reactions$binding_unbinding[which(bub=="unbinding")]
    #names(binding) <- c()
    #names(unbinding) <- c()		
#	print(binding)
#	print(unbinding)

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
    rnsx <- rns[["binding_unbinding"]]
    rnsxn <- getbub(rnsx)

    reactns <- c(unlist(rns[["trans_degradation"]]), rnsxn, unlist(rns[["activation"]]),unlist(rns[["reactions"]]))
    names(reactns) <- c()
    #print(pop)	    
    rnstmp <- c("population {",pop,"}","reaction {",reactns,"}", qind)
    if(!is.null(rns$waitlist)){
        rnstmp <- c(rnstmp, rns$waitlist)
    }	 
    fl <- c(params,rnstmp)
    write.table(fl, file=inp, row.names=FALSE, col.names=FALSE, quote=FALSE)
    exp <- c()
    if(timeseries){
    out <- tryCatch({
        .rsgnscall(inp)
   }, warning=function(war){print("Interrupted by user")}, 
	error=function(err){print("error occured")},
	finally={
	if(file.exists(oup)){	
        	datx=read.table(oup,sep="\t",header=TRUE)
        	exp=datx[,-(1:4)]
	}
	if(file.exists(inp)){file.remove(inp)}
	if(file.exists(oup)){file.remove(oup)}
   })
    }
    else{
    		out <- tryCatch({
        	for(i in 1:sample){
			.rsgnscall(inp)
			if(file.exists(oup)){	
	   			datx=read.table(oup,sep="\t",header=TRUE)
           			expt=datx[nrow(datx),-(1:4)]
				if(class(data)!="rsgns.reactions"){
           			exp <- rbind(exp, (expt[paste("R",V(data@network)$name, sep="")]))
				}
				else{
				  exp <- rbind(exp, expt[data[[1]][,1]])	
				}
    			}			
		}

    	}, warning=function(war){print("Interrupted by user")}, 
        error=function(err){print("error occured")},
        finally={
        	if(file.exists(inp)){file.remove(inp)}
        	if(file.exists(oup)){file.remove(oup)}
  	} )	
           
    }

    if(file.exists(inp)){
	file.remove(inp)
    }
    if(file.exists(oup)){		
    	file.remove(oup)
    }	
    #print(rns)
    if(!timeseries){
        rownames(exp) <- c()
    }
    rns1 <- list()	
    rns1[["expression"]] <- t(exp)
    rns1[["reactions"]] <- rns	
    rns1[["network"]] <- g
    rns1[["waitlist"]] <- rns$waitlist	
    if(!is.null(qind)){	
    	rns1[["waitlist"]] <- qind
    }			
    #if(class(data)!="rsgns.reactions"){
    #    rns@network <- as.matrix(get.adjacency(data@network))
    #    rns <- .get_inipop(g)
    #}
    cat("\n \n Done.\n \n")
    
    rns1
}


.inpop <- function(P=NULL, N=NULL,OP=NULL,ind=TRUE, n1=1, n2=1, n3=1, n4=0){
    init <- c()
    pop1 <- pop2 <- tmp <- NULL
    rna <- paste("R",P,sep="")
    protein <- paste("P",P,sep="")
    pro <- paste("Pro",P,sep="")
    pop1 <- paste(pro,"=",n3)
    names(pop1) <- pro
    if(is.null(OP)){
        OP <- rep(1, length(N))
    }
    if(ind){
        pro <- paste("Pro",P,sep="")
        u_site <- paste("No",N,sep="")
    
	tmp  <- which(OP==2)
	if(length(tmp)>=1){
		u_site <- c(u_site, paste("No",N[tmp],".2",sep=""))
		
	}	
        pop1x <- paste(pro,u_site, sep=".")
        pop1 <- paste(pop1x,"=",n3)
        names(pop1) <- pop1x
    
        pop2x <- paste(pro,N, sep=".")
	if(length(tmp)>=1){
		pop2x <- c(pop2x, paste(pro,N[tmp],"2",sep="."))
		
	}	
        pop2 <- paste(pop2x, "=", n4)
        names(pop2) <- pop2x
    }
        pop3 <- paste(rna, "=", n1)
        names(pop3) <- rna
        pop4 <- paste(protein, "=", n2)
        names(pop4) <- protein
        pop5 <- pop6 <- c()
        if(is.null(pop5)){
            pop5 <- paste("P",N,sep="")
        }
        if(is.null(pop6)){
            pop6 <- paste("R",N,sep="")
        }
	if(length(tmp)>=1){
	    pop5 <- c(pop5, pop5[tmp])
	    pop6 <- c(pop6, pop6[tmp])
	    OP[tmp] <- 1
	    OP <- c(OP, rep(-1, length(tmp)))	
	   	
	}
        init <- list(prosites=pop1,promotors=pop2,rna=pop3,protein=pop4,pop5=pop5,pop6=pop6,op=OP)
        init

}

.get_inipop <- function(g,rc=NULL, delay=FALSE, dly.samp=NULL, dly.rn=NULL, rn.rate.fn=NULL, inhib=NULL, induce=FALSE, indpop=NULL, formula=TRUE){
    if(is.null(V(g)$name)){
        V(g)$name <- paste("G", c(1:vcount(g)),sep="")
    }
    #print(V(g)$name)
    ppop <- rpop <- NULL
    if(!is.null(V(g)$Ppop)){
    	ppop <- V(g)$Ppop
        ppop <- ppop[1:vcount(g)]
        ppop[is.na(ppop)] <- 1
        names(ppop) <- V(g)$name
	}	 	    
    if(!is.null(V(g)$Rpop)){
    	rpop <- V(g)$Rpop
        rpop <- rpop[1:vcount(g)]
        rpop[is.na(rpop)] <- 1
        names(rpop) <- V(g)$name
	}

    inh1 <- inh2 <- inh3  <- NULL
    if(!is.null(inhib)){
    	if(length(inhib)==1){inh1 <- inhib[1]}
    	if(length(inhib)==2){inh1 <- inhib[1]; inh2 <- inhib[2]}
    	if(length(inhib)==3){inh1 <- inhib[1]; inh2 <- inhib[2:3]}
    	if(length(inhib)==4){inh1 <- inhib[1]; inh2 <- inhib[2:4];}
    	if(length(inhib)>4){inh1 <- inhib[1]; inh2 <- inhib[2:4];inh3 <- inhib[5:length(inhib)]}
	#print("init inhb")
	#print(inh3)	
    	
    }

    rnl1 <- rnl2 <- rnl3 <- rnl4 <- NULL
    for(i in 1:length(rn.rate.fn)){
    	if(!is.null(rn.rate.fn[i][[1]])){
		if(i==1){rnl1 <- rn.rate.fn[i]}
		if(i==2){rnl2 <- rn.rate.fn[i]}
		if(i==3){rnl3 <- rn.rate.fn[i]}
		if(i==4){rnl4 <- rn.rate.fn[i]}
		
	}
     }		
	#print(rnl1)
	#print(rnl2)
	#print("#############")
    
    pop1 <- c()
    pop2 <- c()
    pop3 <- c()
    pop4 <- c()
    dg <- degree(g, mode="in")
    dg <- sort(dg, decreasing=TRUE)
    nm <- names(dg)
    bind <- activ <- dgrd <-c()
    rnrfn <- NULL
    for(i in 1:length(nm)){
        ind <- TRUE
        if(dg[i]==0){
            ind <- FALSE
        }
        nb <- neighbors(g,nm[i],mode="in")
        nb <- V(g)$name[nb]
        opx <- rbind(nb, rep(nm[i], length(nb)))
        opx <- get.edge.ids(g, (opx))
        opx <- E(g)$op[opx]
         

        xx <- NULL
        #print(ppop[nm[i]])
	if(is.null(ppop)&&is.null(rpop)){
        	xx <- .inpop(P=nm[i], N=nb,OP=opx)
	}
	if((!is.null(ppop))&&is.null(rpop)){
        	xx <- .inpop(P=nm[i], N=nb,OP=opx, ind=ind, n1=ppop[nm[i]])
	}
	if((is.null(ppop))&&(!is.null(rpop))){
        	xx <- .inpop(P=nm[i], N=nb,OP=opx, ind=ind, n2=rpop[nm[i]])
	}
	if((!is.null(ppop))&&(!is.null(rpop))){
        	xx <- .inpop(P=nm[i], N=nb,OP=opx, ind=ind, n1=ppop[nm[i]], n2=rpop[nm[i]])
	}
        #print(nb) 
	#xx <- .inpop(P=nm[i], N=nb,OP=opx, ind=ind)
        #print(xx)
	#print(opx)
	if(dg[i]!=0){
         bxx <- .conn_bind(xx, rc[4], rc[5], rn.rate.fn1=rnl2, rn.rate.fn2=rnl3, inhib=inh2, formula=formula)
         bind <- c(bind, bxx)
          
         cxx <- .conn_activ(xx,tr=rc[6],dly.rn=dly.rn[1], dly.samp=dly.samp[1], rn.rate.fn=rnl4, inhib=inh3, formula=formula)
         activ <- c(activ, cxx)

        }
        trxx <- .trans_deg(xx, rc[1], rc[2], rc[3], dly.rn= dly.rn[2:3], dly.samp=dly.samp[2:3], rn.rate.fn=rnl1, inhib=inh1, induce=induce, formula=formula)
        #print(trxx)
	names(trxx) <- c("translation","decayR","decayP")
	dgrd <- c(dgrd, trxx)
        
        pop1 <- c(pop1,xx[[1]] )
        pop2 <- c(pop2,xx[[2]] )
        pop3 <- c(pop3,xx[[3]] )
        pop4 <- c(pop4,xx[[4]] )


    }
    indmol <- NULL
    if(induce){
     if(!is.null(indpop)){
	indmol <- paste("ind = ",indpop)
     }	
     else{
	indmol <- paste("ind = ",10)
	}	
}
    if(formula){
	pop <- list(prosites=pop1,promotors=pop2,rna=pop3,protein=pop4, inducer=indmol)
    }	  
    else{
	k <- c(pop1,pop2,pop3,pop4, indmol)
	pop <- matrix(unlist(strsplit(k, " = ")), ncol=2, byrow=TRUE)
    }		
    #pop <- list(pop1, pop2, pop3, pop4, indmol)
   list(population=pop, trans_degradation=dgrd, binding_unbinding=bind, activation=activ)
 	
}

.conn_bind <- function(xx=NULL, br=NULL, ubr=NULL, rn.rate.fn1=NULL, rn.rate.fn2=NULL,inhib=NULL, formula=TRUE){
	in1 <- in2 <- NULL
	if(length(inhib)<=2){

		in1 <- inhib
	}
	if(length(inhib)>2){
		in1 <- inhib[1:2]
		in2 <- inhib[3]
	}

    sub1 <- names(xx[[1]])
    sub2 <- xx[[5]]
    pro1 <- names(xx[[2]])
    tmp1 <- c()
    tmp2 <- c()
    k <- 1
    if(!is.null(rn.rate.fn1)){
    	if(!is.null(rn.rate.fn1[[1]][3][[1]])){
		tmls <- list()
		k <- rn.rate.fn1[[1]][[3]]
		for(i in 1:length(k)){
		 if(k[i]==1){tmls[[i]] <- rn.rate.fn1[[1]][1:2]}
		 if(k[i]==0){tmls[[i]] <- list("","")}
		}	
		rn.rate.fn1 <- tmls
	}
    }
    		

    for(i in 1:length(xx[[1]])){
	if(formula){
        rn1 <- .set.reactions(substrate=c(sub1[i],sub2[i]), product=pro1[i],rconst=br, rn.rate.function=rn.rate.fn1, inhib=in1)
        tmp1 <- c(tmp1,rn1)
	}
	else{
		#rn.rate.fnx <- lapply(rn.rate.fn1, function(x)paste(unlist(x), collapse=":"))
		#print(rn.rate.fn1)
		rn.rate.fnx <- lapply(rn.rate.fn1, .gtf)
		#sapply(p, function(x)paste(x[[1]], paste(x[[2]], collapse=","), sep=":"))
		tmp_vec <- c(substrate=c(sub1[i], sub2[i]), product=pro1[i], rconstant=br, rratefunction=unlist(rn.rate.fnx),inhibitor=in1)
		tmp1 <- rbind(tmp1, tmp_vec)
		if(!is.null(tmp1)){
		rownames(tmp1) <- paste("bind", 1:nrow(tmp1),sep="")
		}
        }
	if(formula){
        rn2 <- .set.reactions(substrate=pro1[i], product=c(sub1[i],sub2[i]),rconst=ubr, rn.rate.function=rn.rate.fn2, inhib=in2) 
	tmp2 <- c(tmp2,rn2)
	}
	else{
		#rn.rate.fnx <- lapply(rn.rate.fn2, function(x)paste(unlist(x[1:2]), collapse=":", sep=""))
	
		rn.rate.fnx <- lapply(rn.rate.fn2, .gtf)
		tmp_vec <- c(substrate=pro1[i], product=c(sub1[i], sub2[i]), rconstant=ubr, rratefunction=unlist(rn.rate.fnx),inhibitor=in2)
		tmp2 <- rbind(tmp2, tmp_vec)
		if(!is.null(tmp2)){
		rownames(tmp2) <- paste("unbind", 1:nrow(tmp2), sep="")
		}
	}
 	}
    list(binding=tmp1,unbinding=tmp2)
    
}

.conn_activ <- function(xx, dly.samp=NULL, dly.rn=NULL, tr=20, rn.rate.fn=NULL, inhib=NULL, formula=TRUE){
   inh1 <- inh2 <- NULL
    if(length(inhib)==1){
     inh1 <- inhib[1]
     inh2 <- ""	
     }	
    if(length(inhib)==2){
     inh1 <- inhib[1]
     inh2 <- inhib[2]
     }
    #print("instart")	
    #print(inhib)	
    #print(c(inh1,inh2))		
    xxt <- rn.rate.fn[[1]][1:2]
    cnt <- xx$op
    cnt <- which(cnt==1)
    cntneg <- which(xx$op==-1)
    #print(dly.samp)	
    n <- length(cnt)
    rns <- c()
    tmls <- list() 
    flag <- 0
    if(!is.null(rn.rate.fn[[1]][3][[1]])){
              inds <- rn.rate.fn[[1]][[3]]
	      flag=1	 
    }
	#rn.rate.fn1[[1]][1:2]
for(i in (2^n-1):1){
	if(n==0){
		break
	}
	inhx1 <- inhx2 <- NULL	
        indx <- intToBits(i)
        indx <- as.integer(indx)
        indx1 <- indx[n:1]

        indx2 <- abs(indx1-1)
        # print(indx1)
        #print(indx2)
	if(flag){
	if((inds[1]==1)&&(inds[2]==0)){
		tmp <- c(rep(list(xxt), length(which(indx1==1))), rep(list(list("","")), c(length(which(indx1==0))+ length(cntneg))))
		rn.rate.fn <- tmp
	}	
	if((inds[1]==0)&&(inds[2]==1)){
		tmp <- c( rep(list(list("","")), length(which(indx1==1)) ), rep(list(xxt), length(which(indx1==0))+length(cntneg) ) )
		rn.rate.fn <- tmp
	}	
	if((inds[1]==1)&&(inds[2]==1)){
		tmp <- rep(list(xxt), c(length(which(indx1==1) )+length(which(indx1==0)+length(cntneg) )))
		rn.rate.fn <- tmp
	}	
	}
	if(!is.null(inh1)){
		inhx1 <- rep(inh1, length(which(indx1==1)))
	}
	if(!is.null(inh2)){
		inhx2 <- rep(inh2, length(which(indx2==1))+length(cntneg))
	}
	#print("####**")
	#print(indx1)
	#print(inhx1)
	#print(inhx2)
	#print("###")	
        mol <- c(names(xx[[2]])[which(indx1==1)], names(xx[[1]])[which(indx2==1)], names(xx[[1]])[cntneg])
        if(formula){
	rnx <- .set.reactions(substrate=mol, product=c(mol, names(xx[[3]])),rconst=tr, dly.rn=dly.rn, dly.samp=dly.samp, rn.rate.function=rn.rate.fn, inhib=c(inhx1,inhx2) )
        rns <- c(rns,rnx)
	}
	else{
		#rn.rate.fnx <- lapply(rn.rate.fn, function(x)paste(unlist(x), collapse=":"))
		rn.rate.fnx <- lapply(rn.rate.fn, .gtf)
		#dly.smpx <- lapply(dly.samp, function(x)paste(unlist(x), collapse=":"))
		dly.smpx <- lapply(dly.samp, .gtf)
		tmp_vec <- c(substrate=mol, product=c(mol, names(xx[[3]])), rconstant=tr,delayrn=unlist(dly.rn), dlysmp=unlist(dly.smpx), rratefunction=unlist(rn.rate.fnx),inhibitor=c(inhx1,inhx2))
		rns <- rbind(rns, tmp_vec)
	}
    
    }
    if(!formula){
    if(!is.null(rns)>0){	
    	rownames(rns) <- paste("reaction", 1:nrow(rns), sep="")	
    }
    rns <- list(rns)
    }		
    rns
}


## translation degreation pre gene ###


.trans_deg <- function(xx, tr=NULL, rdr=NULL, pdr=NULL,dly.rn=NULL, dly.samp=NULL, rn.rate.fn=NULL, inhib=NULL, induce=FALSE, formula=TRUE){
    
    rpa <- c()
    ra <- c()
    pa <- c()
    if(induce){
    rn.rate.fn <- c(list(list("","")),rn.rate.fn)		
    #if(formula){
	#rpa <- c(rpa,.set.reactions(substrate=names(xx[[3]]), product=c(names(xx[[3]]), names(xx[[4]])),rconst=tr, dly.rn=dly.rn, dly.samp=dly.samp, rn.rate.function=rn.rate.fn, inhib=inhib  ))
    #}
    if(formula){
    rpa <- c(rpa,.set.reactions(substrate=c(names(xx[[3]]), "ind"), product=c(names(xx[[3]]), names(xx[[4]])),rconst=tr, dly.rn=dly.rn, dly.samp=dly.samp, rn.rate.function=rn.rate.fn, inhib=c("",inhib)  ))
    }
   else{
	
		rnx <- lapply(rn.rate.fn, .gtf)
		dlnx <- lapply(dly.samp, .gtf)
		rpa <- rbind(rpa, c(substrate=c(names(xx[[3]]),"ind"), product=c(names(xx[[3]]), names(xx[[4]])), rconstant=tr, delayrn=dly.rn , dlysmp=unlist(dlnx), rratefunction=rnx, inhibitor=c("",inhib)))
    }	
    }
    else{
	if(formula){
		rpa <- c(rpa,.set.reactions(substrate=names(xx[[3]]), product=c(names(xx[[3]]), names(xx[[4]])),rconst=tr, dly.rn=dly.rn, dly.samp=dly.samp, rn.rate.function=rn.rate.fn, inhib=inhib  ))
	}
	else{
		#dlnx <- lapply(dly.samp, function(x)paste(unlist(x), collapse=":"))
		#rnx <-  lapply(rn.rate.fn, function(x)paste(unlist(x), collapse=":"))
		rnx <- lapply(rn.rate.fn, .gtf)
		dlnx <- lapply(dly.samp, .gtf)
		rpa <- rbind(rpa, c(substrate=names(xx[[3]]), product=c(names(xx[[3]]), names(xx[[4]])), rconstant=tr, delayrn=dly.rn , dlysmp=unlist(dlnx), rratefunction=rnx, inhibitor=inhib))
	}		
	}
    if(formula){	
    ra <- c(ra, .set.reactions(substrate=names(xx[[3]]),rconst=rdr) )
    pa <- c(pa, .set.reactions(substrate=names(xx[[4]]),rconst=pdr) )
    }
    else{
	ra <- rbind(ra,c(substrate=names(xx[[3]]), product=c("",""), rconstant=rdr))
	pa <- rbind(pa,c(substrate=names(xx[[4]]), product=c("",""), rconstant=pdr))

    }	
    list(rpa, ra, pa)
    
}
.gtf <- function(x){
    #print((x))
     rn <- ""
    if((x[[1]]=="")||is.null(x)){
        rn <- ""
    }
    else{
       if(!is.null(x[[2]])){	
       	rn <- paste(x[[1]], paste(x[[2]], collapse=","), sep=":")
       }
      if(length(x[[2]])==1){
       	rn <- paste(unlist(x), collapse=":")
      }  		
      	
     }
    rn
}


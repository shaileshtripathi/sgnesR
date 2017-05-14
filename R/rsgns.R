
setOldClass("igraph")
#rsgns.data <- 
setClass("rsgns.data", slots= c(network="igraph", rconst="numeric", submol= "numeric", promol="numeric", dly.rn="numeric", init_pop="numeric",
	pop_op="character", pop="numeric",dly.samp="list", rn.rate.function="list", inhib="character"))

#rsgns.param <- 
setClass("rsgns.param", slots= c(time = "numeric", stop_time="numeric", readout_interval = "numeric",  
	save_interval="numeric", coop="numeric"))

#rsgns.waitlist <- 
setClass("rsgns.waitlist", slots= c(nodes="numeric", time = "numeric", mol="numeric", type="character"))
#setClass("rsgns.reactions", slots= c(reactions="character", decay = "character", population="character", expression="matrix",network="matrix", waitlist="character"))

plot.sgnesR <- function(x,mol="R", method=c("boxplot","exp"), gn=NULL,ln=FALSE,...){
		#if(class(x)=="rsgns.reactions"){
		#	k <- (x@network)
		#	k <- Matrix::t(k)+k
		#	k1 <- which(k[upper.tri(k)]>=1)
		#	ex <- x@expression
		#	ex <- ex[rownames(k),]
		#	 
		#	cr <- cor(t(ex), method=method)
		#	k2 <- which(k[upper.tri(k)]==0)
		#	lst <- list(abs(cr[upper.tri(cr)][k1]), abs(cr[upper.tri(cr)][k2]))
		#	names(lst) <- c("Edges", "Non- edges")
		#	boxplot(lst, col=c("red","blue"), cex.axis=1.5, cex.label=2)
		#
		#}
		g <- x$network
		exp <- c()
		nm <- c()
		k <- grep(mol, rownames(x$expression))
		exp <- x$expression[k,]
		if(ln){
			exp <- log(exp)
		}
		#print(exp)	
		if(method[1]=="boxplot"){
			boxplot(t(exp))
		}
		if(method[1]=="exp"){
				
			if(is.null(gn)){
				gn <- which.max(degree(g, mode="in"))
				nb <- as.vector(neighbors(g,gn, mode="in"))
				nm <- paste(mol,V(g)$name[c(gn,nb)],sep="")
				
			}
			else{
				
				nm <- paste(mol,V(g)$name[c(gn)],sep="")
			}
			
			if(length(col)<length(nm)){
				col <- rep(col, length(nm))
			}
			
			if(length(lty)<length(nm)){
				lty <- rep(lty, length(nm))
			}
			if(length(nm)>1){
			plot(exp[nm[1], ], xlab="Time", ylab="Expression")

			for(i in 2:length(nm)){
				lines(exp[nm[i], ], type=type, cex=cex, pch=pch, col=col[i], lwd=lwd, lty=lty[i])
			}
			}
			else{
			plot(exp[nm[1], ],xlab="Time", ylab="Expression")
			}
		}

}

.getreactionsx <- function(data, waitlist, decay=TRUE, decayrate=1,  pop=NULL, pop.up=NULL, rev=FALSE, wait=FALSE){
	if(is.null(V(data@network)$name)){
			V(data@network)$name <- paste("N",c(1:vcount(data@network)),sep="")

	}
	setpop <- .set.population(V(data@network)$name, data@pop, data@pop_op)
	el <- get.edgelist(data@network)
	if(decay){
		if(length(data@decayrate)==0){
			data@decayrate <- rep(1, vcount(data@network))
		}
	}
	if(length(data@rconst)==1){
			data@rconst <- rep(data@rconst, ecount(data@network))
	}
	#if(data@rconst)
	rn <-c()
	for(i in 1:nrow(el)){
		rn <- c(rn, .set.reactions(el[i,1],el[i,2], data@rconst[i], data@submol[i], data@promol[i], data@dly.rn[i], data@dly.samp, data@rn.rate.function, data@inhib[i]) )
	}
	rn1<-c()
	if(decay){
		if(rev){
		for(i in 1:vcount(data@network)){
		 rn1 <- c(rn1, .set.reactions("",V(data@network)$name[i], data@decayrate[i]))
		}
		}
		else{
			for(i in 1:vcount(data@network)){
		    rn1 <- c(rn1, .set.reactions(V(data@network)$name[i],"", data@decayrate[i]))
		}

		}
	}
	else{rn1=""}
	sw <- ""
	if(wait){
		sw <-  .set.waitinglist(V(data@network)$name[waitlist@nodes], t=waitlist@time, equilibrium=waitlist@mol)
	}
	rnob <- new("rsgns.reactions", reactions=rn, decay=rn1, population=setpop, waitlist=sw)
	rnob
}

#rsgns <- function(data, param, waitlist, sample=1, decay=TRUE, wait=FALSE, rev=FALSE, timeseries=TRUE){
#	if(class(data)!="rsgns.reactions")
#	{
#	if(is.null(V(data@network)$name)){
#			V(data@network)$name <- paste("N",c(1:vcount(data@network)),sep="")
#
#	}
#	}	
#	params <- c()
#	rnd <- paste( sample( 0:9, 20, replace=TRUE ), collapse="" )
#	oup <- file.path(".",paste(".rsgns_output",rnd,sep=""))
#	#oup <- "tmp1.txt"
#	savef <- tempfile()
#	savef="tmp2"
#	inp <- file.path(".",paste(".rsgns_input",rnd,sep=""))
#	#inp <- "tmp.txt"
#
#	
#	
#	params <- rbind(params, paste("time", param@time,";"))
#	params <- rbind(params, paste("stop_time ", param@stop_time,";", sep=" ") )
#	params <- rbind(params, paste("readout_interval ", param@readout_interval,";", sep=" ") )
#	params <- rbind(params, paste("output_file ", oup,";", sep="") )
#	#params <- rbind(params, paste("save_interval ", param@save_interval,";", sep=" ") )
#	#params <- rbind(params, paste("save_index ", param@save_index,";", sep=" ") )
#	#params <- rbind(params, paste("save_file ", savef,"_%%.sim;", sep="") )
#	params <- rbind(params, paste("lua!{\n coop = ", param@coop, ";\n}!"))
#	#reactions <- getreactions(data)	
#	#reactions <- rbind(paste("reaction {\n"), reactions)
#	#reactions <- rbind(reactions, "\n }")
#	rns <- NULL
#	if(class(data)!="rsgns.reactions"){
#		rns <- getreactions(data, waitlist, decay=decay, rev=rev, wait=wait)
#	}
#	else{
#		rns <- data
#	}
#	rnstmp <- c( "reaction {",rns@reactions,rns@decay,"}", "population {",rns@population,"}")
#
#	
#  fl <- c(params,rnstmp,rns@waitlist)
#	write.table(fl, file=inp, row.names=FALSE, col.names=FALSE, quote=FALSE)
#    exp <- c()
#    if(timeseries){
#  	    .rsgnscall(inp)
#    	datx=read.table(oup,sep="\t",header=TRUE)
#    	exp=datx[,-(1:4)]
#    }
#    else{
#    for(i in 1:sample){
#    .rsgnscall(inp)
#	datx=read.table(oup,sep="\t",header=TRUE)
#    expt=datx[nrow(datx),-(1:4)]
#	exp <- rbind(exp, (expt[V(data@network)$name]))
#	}
#	}
#	#print(oup)
#	file.remove(inp)
#	file.remove(oup)
#	#print(rns)
#	if(!timeseries){
#	rownames(exp) <- c()
#	}
#	rns@expression <- t(exp)
#	if(class(data)!="rsgns.reactions"){
#		rns@network <- as.matrix(get.adjacency(data@network))
#	}
#	cat("\n \n Done.\n \n")
#	
#	rns
#}
#prepare data
#library(igraph)
#g <- barabasi.game(20, directed=FALSE)
#rsg <- new("rsgns.data",network=g, rconst=25, inhib=c(rep("*",vcount(g))), rn.rate.function=list("invhill", c(10,2)))
#rp <- new("rsgns.param", time=1, stop_time=100, readout_interval=10, save_interval=200, save_index=5,coop=10)
#rwl <- new("rsgns.waitlist",nodes=sample(vcount(g),8), time=sample(rp@stop_time,8), mol=20) 
#xx <- rsgns(rsg, rp, rwl, wait=TRUE,sample=20)

##############################################################




.rsgnscall <- function(inp=NULL){
	if(!is.null(inp)){
		tmp <- c("--include", inp, "--progress")
		k <- .Call("getArgs", as.numeric(length(tmp)), tmp)
		
	}
}

.onLoad <- function(libname, pkgname) {
  library.dynam("sgnesR", pkgname, libname, local=FALSE);
}

.onUnload <- function(libpath) {
  library.dynam.unload("sgnesR", libpath)
}



.setmol <- function(mol=NULL, mol.no= NULL, dly.rn=NULL, dly.samp=NULL, rn.rate.function=NULL, inhib=NULL){
	subm <- dlyrn <- param <-  inhibp <- c()
	rfn <- dlysmp <- list()
	if(!is.null(mol.no)){if(is.na(mol.no)){mol.no <- NULL}}
	if(!is.null(dly.rn)){if(is.na(dly.rn[1])){dly.rn=NULL}}
	#if(!is.null(dly.samp)){if(is.na(dly.samp)){dly.samp=NULL}}
	#if(!is.null(rn.rate.function)){if(is.na(rn.rate.function)){rn.rate.function=NULL}}
	if(is.null(mol)){
		mol <- ""
	}
	if(!is.null(inhib)){if(is.na(inhib[1])){inhib <- NULL}}

	if(is.null(mol.no)){
		subm <- rep("", length(mol))
	}
	else{
		if(length(mol.no)==length(mol)){
			subm <- mol.no
		}
		else{
			subm <- c(mol.no, rep("", (length(mol)-length(mol.no))))
		}
	}
	#print(!is.null(dly.rn) && is.null(dly.samp))
	if(!is.null(dly.rn)){
		if(length(dly.rn)==length(mol)){
			dlyrn <- dly.rn
		}
		if(length(dly.rn)<length(mol)){
			dlyrn <- rep(dly.rn, length(mol))
			dlyrn <- dlyrn[1:length(mol)]
			#dlyrn <- rep(list(dly.rn), length(mol))
			#dlyrn <- as.list(dly.rn)
		}
		#else{
		#	dlyrn <- c(dly.rn, rep("", (length(mol)-length(dly.rn))))
		#}
	}
	else{
		dlyrn <-  rep("", length(mol))
	}
	tmsmp <- list()
	if(!is.null(dly.samp) && (length(dly.samp)!=0)) {
			if(length(dly.samp)==0){
				dly.samp <- list(list("",""))
			}
			if(class(dly.samp) == "list"){
				#print(dly.samp)
				for(i in 1:length(dly.samp)){
					tmsmp[["function"]] <- dly.samp[[i]][[1]]
					tmsmp[["parameter"]] <- dly.samp[[i]][[2]]
					dlysmp[[i]] <-  tmsmp  
				}
			}
	}
	else{
		dlysmp[["function"]] <- ""
		dlysmp[["parameter"]] <- ""
	}
	
	#if(!is.null(rn.rate.function) && (length(rn.rate.function)!=0)){
	#	if(class(rn.rate.function) == "list"){
	#			rfn[["function"]] <- rn.rate.function[[1]]
	#			rfn[["parameter"]] <- rn.rate.function[[2]]
	#
	#	}
	#}
	#else{
	#	rfn[["function"]] <- ""
	#	rfn[["parameter"]] <- ""
	#}
#########################		

	rnsmp <- list()
	if(!is.null(rn.rate.function) && (length(rn.rate.function)!=0)) {
			if(class(rn.rate.function) == "list"){
				#print(dly.samp)
				for(i in 1:length(rn.rate.function)){
					rnsmp[["function"]] <- rn.rate.function[[i]][[1]]
					rnsmp[["parameter"]] <- rn.rate.function[[i]][[2]]
					rfn[[i]] <-  rnsmp  
				}
			}
	}
	else{
		rfn[["function"]] <- ""
		rfn[["parameter"]] <- ""
	}
	


##########################
	if(is.null(inhib)){
		inhibp<- rep("", length(mol))
	}
	else{
		if(length(mol)==length(inhib)){
			inhibp <- inhib
		}
		else{
			inhibp <- c(inhib, rep("", (length(mol)-length(inhib))))
		}
	}

	list(mol, subm, dlyrn, dlysmp, rfn, inhibp)

}

#########################################
#molsubs <- setmol(c("A","B","C"), 2,NULL,dly.samp=NULL, rn.rate.function=list(list("a",c(1:2)), list("x", c(3,5))))
#molsubs <- setmol(c("A","B","C"), 2,c(20,30),dly.samp=NULL)

.set.waitinglist <- function(compound=NULL, t=NULL, equilibrium=NULL){
	if(length(t)==0){
		t=0
	}
	if(length(equilibrium)==0){
		equlibrium=0
	}
	eq = equilibrium
	pr <- paste("queue", paste("[",eq,"]",compound,"(",t,");", sep="") )
	pr
}

.set.population <- function(mol=NULL, pop=NULL, pop.op=NULL){


	init_pop=NULL
	if(is.null(pop)||(length(pop)==0)){
		pop <- rep(0,length(mol))
	}
	if((length(pop))<(length(mol))) {pop <- c(pop, rep(0,(length(mol)-length(pop)) ))}

	if(is.null(pop.op) || (length(pop.op)==0)){
		m <- cbind(mol, pop)
		init_pop=apply(m,1,function(x)paste(x, collapse="="))
	}
	if(!is.null(pop.op) && (length(pop.op)>0) ){
		op_vec=NULL
		op_vec <- pop.op
		if(length(pop.op)==1){
			op_vec <- c(pop.op, rep("", (length(mol)-length(pop.op))) )
		}
		if(length(pop.op)>1 && (length(pop.op)<length(pop))){
			op_vec <- c(pop.op, rep("", (length(mol)-length(pop.op))))

		}
		m <- cbind(mol, op_vec, pop)
		init_pop=apply(m[,1:2],1,function(x)paste(x, collapse=" "))
		init_pop=apply(cbind(init_pop,m[,3]),1,function(x)paste(x, collapse="="))

	}
	init_pop <- paste(init_pop,";", sep="")
	init_pop

}

.getreactionvec <- function(rnfn){
	rnfnvec <- c()
	#print("#######")
	#print(rnfn)
	for(i in 1:length(rnfn)){
		if((rnfn[[i]][[1]]!="") && (rnfn[[i]][[2]]!="")){
		rnfnvec <- c(rnfnvec, paste("(",rnfn[[i]][[1]], ":", paste(rnfn[[i]][[2]],collapse=","),")", sep=""))
		}
		else{
			rnfnvec <- c(rnfnvec, "")
		}
	}
	#print(rnfnvec)
	#print("abc")	
	#print("#######")
	rnfnvec
}
.getreactionlist <- function(mol=NULL, rnfn=NULL){
	xxtmp <- list("","")
	rn <- list()
	rnls <- c()
	if(class(rnfn[[1]])!="list"){
		rnfn <- list(rnfn)
		if(length(mol)>1){
		rnfn <- c(rnfn, rep(list(xxtmp), (length(mol)-1)))
		}
	}
	if((class(rnfn)=="list")&&(length(rnfn)< length(mol))){
		for(i in 1:length(mol)){
			#print(i)
			if(!is.null(rnfn[i][[1]])){
				rn[["function"]] <- rnfn[[i]][[1]]
				rn[["parameter"]] <- rnfn[[i]][[2]]
			}
			else{
				#print("yes")
				rn[["function"]] <- rnls[[(i-1)]][[1]]
                                rn[["parameter"]] <- rnls[[(i-1)]][[2]]

			}
		rnls <- c(rnls, list(rn))
		#print(rnls)
		}
	
	}
	else{
	for(i in 1:length(mol)){
		tmp <- list(list())
		#rnls[[i]][[1]] <- rnfn[[i]][[1]]
		#rnls[[i]][[2]] <- rnfn[[i]][[2]]
		tmp[[1]][["function"]] <- rnfn[[i]][[1]]
		tmp[[1]][["parameter"]] <- rnfn[[i]][[2]]
		rnls <- c(rnls, tmp) 
	}
	}
	rnls

}
.getreaction <- function(mols){
	sm <- mols[[1]]
	mol <- mols[[2]]
	del <- mols[[3]]
	delvec <- c()
	delfn <- mols[[4]]
	for(i in 1:length(del)){
		if(del[i]!=""){
			delvec <- c(delvec, paste("(", del[i],")", sep=""))
		}
		else{
			delvec <- c(delvec,"")
		}
	}
	rnfn <- mols[[5]]
	inhb <- mols[[6]]
	op <-  c(rep("+", (length(sm)-1)),"")

	rnfnvec <- .getreactionvec(.getreactionlist(mol, rnfn))
	#print("##")
	#print(delfn)
	#print("##")
	delfnvec <- .getreactionvec(.getreactionlist(mol, delfn))

	#print(delfnvec)	
	lhs <- paste(inhb,mol, sm, delvec, delfnvec,  rnfnvec,op, collapse="", sep="" )
	#print(lhs)
	lhs
}


.set.reactions <- function(substrate=NULL, product=NULL, rconst=NULL, sub.mol=NULL, pro.mol=NULL, dly.rn=NULL, dly.samp=NULL, rn.rate.function=NULL, inhib=NULL){
	
	molsubs <- .setmol(mol=substrate, mol.no=sub.mol, dly.rn=NULL, rn.rate.function=rn.rate.function, inhib=inhib)
	molprod <- .setmol(mol=product, mol.no=pro.mol, dly.rn=dly.rn, dly.samp=dly.samp)
	if(!is.null(rconst)){if(is.na(rconst)){rconst=""}}
	if(is.null(rconst)){rconst <- 1}
	tmp1 <- .getreaction(molsubs)
	#print("mprd")
	#print(molprod)
	tmp2 <- .getreaction(molprod)
	#print(molprod)
	if(!is.null(rconst)){if(is.na(rconst)){rconst=""}}

	paste(tmp1,"--","[",rconst,"]", "-->",tmp2, ";", sep="")
	

}






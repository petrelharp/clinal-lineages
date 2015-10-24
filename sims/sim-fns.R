
getRuns = function(index){
    #This returns an array in which rows are the the beggingns and end of a TRUE/FALSE vector
    starts  = which(c(0,index)[-1]==1 & c(0,index)[-length(index)]!=1)
    stops   = which(c(index,0)[-1]!=1 & c(index,0)[-(1+length(index))]==1)
    if(index[length(index)] == 1 & index[(length(index)-1)] !=1){starts = c(starts,length(index))}
    return (cbind(starts,stops))
}


namesList = function(these.names){
    # this is a way to get a list where each elelment has a name equalto its value, so that when I lapply, I get back a [correctly] named list
    all.things = as.list(these.names)
    names(all.things)=these.names
    return(all.things)
}

#get breaks between species' ancestry
spBreaks = function(sim){
    ind.ancest = lapply( sim$inds,function(IND){ 
            lapply(IND,function(CHR){
               lapply(CHR,function(PAR){
                    if(nrow(PAR)==1){return(data.frame(starts=0,stops=1,sp2 = PAR$ancest > length(sim[[1]])))}
                    sp1 = getRuns(as.numeric( PAR$ancest <= length(sim[[1]])))
                    sp2 = getRuns(as.numeric( PAR$ancest > length(sim[[1]])))
                    tmp = rbind(data.frame(sp1,sp2=rep(FALSE,nrow(sp1))),data.frame(sp2,sp2=rep(TRUE,nrow(sp2))))
                    tmp = tmp[order(tmp$starts),]
                    tmp$starts = PAR$pos[tmp$starts]
                    tmp$stops = c(PAR$pos[-1],1)[tmp$stops]
                    return(tmp)
                })
            })
        })
    ancest.prop = t(sapply(ind.ancest,function(IND){
        tmp = t(sapply(IND,function(CHR){ c(X1=with(CHR$X1,sum(sp2*(stops-starts))), X2 = with(CHR$X2,sum(sp2*(stops-starts))))  }))
        return(colSums(tmp)/nrow(tmp))
    }))		
    return(list(ind.ancest = ind.ancest, ancest.prop = ancest.prop))
}


initializeInds = function(arr,n.chr=n.chr){	
    return(lapply(seq_along(arr[,1]),function(IND){
        X=replicate(n.chr,list( X1=data.frame(pos=0,ancest = arr[IND,"X1"]), 
            X2=data.frame(pos=0,ancest = arr[IND,"X2"])  ),simplify=FALSE)
        names(X) = paste("chr",c(1:n.chr),sep="")
        return(X)
    }))
}


#######  GOING THROUGH MEIOSIS    #######
meiosis = function(CHR,r=1){
    start = sample(c(0,1),1)
    breaks = sort(runif(rpois(1,1)))
    A= breaks[seq_along(breaks)%%2==start]
    B= breaks[seq_along(breaks)%%2!=start]
    X = rbind(
        do.call(rbind,lapply(A,function(a){
            with(CHR[[1]],c(pos=a,ancest=ancest[max(which(pos<a))]))})),
            CHR[[1]][(findInterval(CHR[[1]][,"pos"],c(-1,breaks,1.2))%%2)!=start,],
            do.call(rbind,lapply(B,function(b){with(CHR[[2]],c(pos=b,ancest=ancest[max(which(pos<b)) ]))})),
            CHR[[2]][(findInterval(CHR[[2]][,"pos"],c(-1,breaks,1.2))%%2)==start,]
        )
    return(		X[order(X$pos),]		)
}
#######  GENOTYPING INDIVIDUALS   #######
genotype = function(IND.CH,L,sp){
    do.call(cbind,lapply(IND.CH,function(PAR){	
        ancest.chunk = with(PAR,ancest[findInterval( L, c(pos,1) )]  )
        sp[ancest.chunk]  
    }))
} 	
genoAll = function(inds,QTL,sp.ids){
    tmp.ch = namesList(names(QTL))
    sp=rep(sp.ids,each=2)
    lapply(inds,function(IND){
        lapply(tmp.ch,function(CHR){
            genotype(IND[[CHR]],QTL[[CHR]]$pos,sp)
        })
    })
}



#######   CALCULATING FITNESS   #######

#######   ASSORTATIVE MATING    #######

#######    MAKING A NEW GEN     #######

makeNewGen=function(inds,QTL,sp.ids,demes,edge="M/2",sigma,max.deme=max(demes),min.deme=min(demes)){
    pop.sizes <- table(demes) # number of individuals in each deme
    chrnames <- namesList( names(inds[[1]]) )
    #genotyping
    if(!is.null(QTL)){ 
        genos = lapply(genoAll(inds=inds,QTL=qtl,sp.ids=sp.ids) ,function(IND){do.call(rbind,IND)})  
        QTL=do.call(rbind,QTL)
        #make fitness from genotypes
        w_abs = sapply(genos,function(IND){
            #NOTE I WILL CHANGE THIS IF THERE ARE EPISTSTIC LOCI
            loc.vals = sapply(seq_along(IND[,1]),function(LOCUS){
                if(sum(QTL[LOCUS,"traits",drop=FALSE] == "additive_fitness") == 1){
                    return(  1+as.numeric(IND[LOCUS,] == "A") * QTL[LOCUS,"s"]  )
                }	
                if(sum(QTL[LOCUS,"traits",drop=FALSE] == "underdominant") == 1){
                    return( 1+QTL[LOCUS,"s"] *(-1)^(sum(IND[LOCUS,] == "A")) )
                }
            })
            return(prod(unlist(loc.vals)))
        })
    }
    if(is.null(QTL)){ w_abs = rep(1,length(inds)) }
    rel.w = w_abs/mean(w_abs)
    inds.df = data.frame(w=rel.w,deme=demes,edge=ifelse(demes%in%range(demes),2,1))
    # sample the displacements of individuals
    sample_move = floor(rnorm(NROW(inds.df),0.5,sigma))
    # new locations of individuals
    inds.df$new.deme = pmin(max.deme, pmax(min.deme, inds.df$deme + sample_move ) )
    inds.df$id = seq_along(inds.df$new.deme)
    new.inds.df = do.call(rbind,lapply(  unique(inds.df$deme), function(d){	# iterate over demes
                these.parents = which(inds.df$new.deme==d)
                # HACK: if there's no available parents, pick ome
                while (length(these.parents)==0) {
                    parent.probs = inds.df$w*exp(-(inds.df$deme-d)^2/sigma)
                    these.parents = sample.int(nrow(inds.df),1,prob=parent.probs)
                }
                if (length(these.parents)>1) {  # this is needed because of sample()'s bad behavior
                    #find parents
                    a = with(inds.df[these.parents,], sample(id,size=pop.sizes[d], prob=w, replace=TRUE))
                    # MODIFY B|A TO ALLOW FOR ASSORTATIVE MATING & eg no selfing
                    b = with(inds.df[these.parents,], sample(id,size=pop.sizes[d], prob=w, replace=TRUE))
                } else {
                    a = b = rep(inds.df[these.parents,"id"][1],pop.sizes[d])  # will be NA if no available parents
                }
                return(cbind(a,b))
    }))
    lapply(1:nrow(new.inds.df),function(k){
       if (any(is.na(new.inds.df[k,]))) {
           return( list( NA ) )
       } else {
           sperms <- lapply( inds[[new.inds.df[k,1]]], meiosis )
           eggs <- lapply( inds[[new.inds.df[k,2]]], meiosis )
           return( lapply( chrnames, function (chrom) {
                  list( X1=sperms[[chrom]], X2=eggs[[chrom]] )
                } ) )
       }
    })
}

#	calculating genome-wide ancetry proportions for each ind
getAlphas = function(inds,sp.ids){
    sp=rep(sp.ids,each=2)
    do.call(rbind,lapply(inds,function(IND){
        X=do.call(rbind,lapply(IND,function(CHR){
            do.call(rbind,lapply(CHR,function(PAR){
                sp.origin=factor(sp[PAR$ancest],levels=unique(sp))
                tapply(diff(c(PAR$pos,1)),sp.origin,sum)
            }))
        }))
        X[is.na(X)]=0
        colMeans(X)  
    }))
}

sim.gens = function(
                    n,
                    sp.ids,
                    deme,
                    n.gen,
                    n.chr,
                    QTL, 
                    sigma,
                    SAVE=NULL,
                    RETURN=TRUE,
                    quiet=FALSE
            ){
    initial.inds.df = data.frame(matrix(1:(2*(n)),ncol=2,byrow=TRUE),sp.ids=sp.ids,deme=deme)
    initial.inds = initializeInds(arr=initial.inds.df,n.chr=n.chr)
    inds = initial.inds
    inds.demes = initial.inds.df$deme
    for(i in 1:n.gen){
        inds = makeNewGen(inds,QTL=qtl,sp.ids=sp.ids,demes=inds.demes,sigma=sigma)  # add in selection and assortative mating
        #plot(getAlphas(inds,sp.ids)[,1])
        if (!quiet) { print(c(i,date())) }
        if(!is.null(SAVE)){ if(i%in%SAVE$gen){  
                sim = list(inds=inds, sp.ids=sp.ids,deme=deme,pars = c(sigma=sigma,n.gen=n.gen),QTL=QTL )
                save(sim,file=sprintf("%s_t%i_g%i.Robj",SAVE$name,SAVE$try,i))
        }}
    }
    if(RETURN){return(list(inds=inds, sp.ids=sp.ids,deme=deme,pars = c(n.gen=n.gen),QTL=QTL ))}
    return(NULL)
}		


# Simulate a hybrid zone with a total of nind individuals
# across n.demes, for ngens,
# with the left half called 'A' and the right half called 'B'.
sim.zone <- function (n.ind, n.deme, ...) {
	n.ind <- ceiling(n.ind/n.deme)*n.deme
	sim.gens( n=n.ind, 
             sp.ids = rep(c("A","B"),each=n.ind/2), 
             deme=rep(1:n.deme,each=floor(n.ind/n.deme)),
             n.chr=1,
             QTL = qtl,...)
}


#GENOTYPE  [ASSUMES A ONE CHR GENOME WHERE WE HAVE PHASED DIPLOID DATA]
geno.pop <- function (ind.ancest, loci) { 
    sapply(ind.ancest, geno.ind, loci) 
}

geno.ind <- function(IND,loci) {
    # Each individual is a list of chromosomes; each chromosome is a two-column data frame of starts, stops, and ancestral identities
    # Returns diploid genotypes;
    #  as a vector, with all chromosomes squooshed together
    if (length(IND)!=length(loci)) { stop("Number of chromsomes not matching.") }
    unlist( lapply( seq_along(IND), function (CHR) {
            rowSums( sapply(IND[[CHR]],function(HAP){
                   HAP$sp2[ findInterval( loci[[CHR]], c(HAP$starts,1), rightmost.closed=TRUE ) ]
              }) )
        }) )
}


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
    ind.ancest = lapply( sim$inds,function(IND){ lapply(IND,function(CHR){lapply(CHR,function(PAR){
            if(nrow(PAR)==1){return(data.frame(starts=0,stops=1,sp2 = PAR$ancest > length(sim[[1]])))}
            sp1 = getRuns(as.numeric( PAR$ancest <= length(sim[[1]])))
            sp2 = getRuns(as.numeric( PAR$ancest > length(sim[[1]])))
            tmp = rbind(data.frame(sp1,sp2=rep(F,nrow(sp1))),data.frame(sp2,sp2=rep(T,nrow(sp2))))
            tmp = tmp[order(tmp$starts),]
            tmp$starts = PAR$pos[tmp$starts]
            tmp$stops = c(PAR$pos[-1],1)[tmp$stops]
            return(tmp)
    })})})
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
    ##Determine if an individual moves up/down a deme by binomial sampling. 
    ##The prob vector in the first rbinom determines whether an individual is in an edge deme or not. This is done because an individual in the leftmost deme cannot be sent left, and vice versa.
    ##The first rbinom gives the probability of moving right, if a move is made. The sign then just indicates which direction the move occurs. The second rbinom determines whether or not an individual moves.
    #move = sign(rbinom(nrow(inds.df),1,prob=with(inds.df,ifelse(deme==min(deme),1,ifelse(deme!=max(deme),1/2,0))))-.5) * with(inds.df, rbinom(nrow(inds.df),1, c(mig/edge)))
        #MOD (Alisa): Change migration so that we can give a vector or migration probabilities e.g. c(+1 deme, +2 demes)
        #Decide whether you stay or move n steps (MAKE SURE mig vector DOES NOT include probability that you stay!)
        #JUMP_SIZE = apply(inds.df,1,function(X){sample(0:(length(mig)),1,prob=c(1-sum(mig/as.numeric(X[3])),mig/as.numeric(X[3])))})
        #move = sign(rbinom(nrow(inds.df),1,prob=with(inds.df,ifelse(deme==min(deme),1,ifelse(deme!=max(deme),1/2,0))))-.5) * JUMP_SIZE

        #NEW JUMP_SIZE: sample from Gaussian, send to closest deme. 
        sample_move = floor(rnorm(NROW(inds.df),0.5,sigma))
        
    #new locations of individuals
    inds.df$new.deme = pmin(max.deme, pmax(min.deme, inds.df$deme + move ) )
    inds.df$id = seq_along(inds.df$new.deme)
    new.inds.df = do.call(rbind,lapply(  unique(inds.df$deme)  , function(d){			
        #find parents
        a = with(inds.df[inds.df$new.deme==d,], sample(id,size=length(which(demes==d)), prob=w/sum(w) , replace=T))
        # MODIFY B|A TO ALLOW FOR ASSORTATIVE MATING & eg no selfing
        b = with(inds.df[inds.df$new.deme==d,], sample(id, size=length(which(demes==d)), prob=w/sum(w) , replace=T))
        return(cbind(a,b))
    }))
    lapply(seq_along(new.inds.df[,1]),function(NEW){
        NEW = lapply(list(X1 = inds[[new.inds.df[NEW,1]]],X2 = inds[[new.inds.df[NEW,2]]]),function(PAR){  lapply(PAR,meiosis)   })
        lapply(namesList(names(NEW[[1]])),function(C){lapply(NEW,function(PARNTL){PARNTL[[C]]})})
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

sim.gens = function(n,sp.ids,deme,n.gen,n.chr,QTL, sigma,SAVE=NULL,RETURN=TRUE){
    initial.inds.df = data.frame(matrix(1:(2*(n)),ncol=2,byrow=TRUE),sp.ids=sp.ids,deme=deme)
    initial.inds = initializeInds(arr=initial.inds.df,n.chr=n.chr)
    inds = initial.inds
    inds.demes = initial.inds.df$deme
    for(i in 1:n.gen){
        inds = makeNewGen(inds,QTL=qtl,sp.ids=sp.ids,demes=inds.demes,sigma=sigma)  # add in selection and assortative mating
        #plot(getAlphas(inds,sp.ids)[,1])
        print(c(i,date()))
        if(!is.null(SAVE)){ if(i%in%SAVE$gen){  
                sim = list(inds=inds, sp.ids=sp.ids,deme=deme,pars = c(n.gen=n.gen),QTL=QTL )
                save(sim,file=sprintf("%s_t%i_g%i.Robj",SAVE$name,SAVE$try,i))
        }}
    }
    if(RETURN){return(list(inds=inds, sp.ids=sp.ids,deme=deme,pars = c(n.gen=n.gen),QTL=QTL ))}
    return(NULL)
}		


# Simulate a hybrid zone with a total of nind individuals
# across ndemes, for ngens,
# with the left half called 'A' and the right half called 'B'.
sim.zone <- function (nind, ndeme, ngens,sigma) {
	nind <- ceiling(nind/ndeme)*ndeme
	sim.gens( n=nind, sp.ids = rep(c("A","B"),each=nind/2),deme = rep(1:ndeme,each=floor(nind/ndeme)),n.chr=1,n.gen=ngens, QTL = qtl,sigma=sigma)
}


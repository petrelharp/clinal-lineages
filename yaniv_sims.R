#setwd("~/Desktop/swarmz")
#library(multicore)
# MAKE STARTING SAMPLES


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
				tmp$starts = PAR$bp[tmp$starts]
				tmp$stops = c(PAR$bp[-1],1)[tmp$stops]
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
			X=replicate(n.chr,list( X1=data.frame(bp=0,ancest = arr[IND,"X1"]), 
				X2=data.frame(bp=0,ancest = arr[IND,"X2"])  ),simplify=FALSE)
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
				with(CHR[[1]],c(bp=a,ancest=ancest[max(which(bp<a))]))})),
				CHR[[1]][(findInterval(CHR[[1]][,"bp"],c(-1,breaks,1.2))%%2)!=start,],
				do.call(rbind,lapply(B,function(b){with(CHR[[2]],c(bp=b,ancest=ancest[max(which(bp<b)) ]))})),
				CHR[[2]][(findInterval(CHR[[2]][,"bp"],c(-1,breaks,1.2))%%2)==start,]
			)
		return(		X[order(X$bp),]		)
	}
#######  GENOTYPING INDIVIDUALS   #######
	genotype = function(IND.CH,L,sp){
		do.call(cbind,lapply(IND.CH,function(PAR){	
			ancest.chunk = with(PAR,ancest[findInterval( L, c(bp,1) )]  )
			sp[ancest.chunk]  
		}))
	} 	
	genoAll = function(inds,QTL,sp.ids){
		tmp.ch = namesList(names(QTL))
		sp=rep(sp.ids,each=2)
		lapply(inds,function(IND){
			lapply(tmp.ch,function(CHR){
				genotype(IND[[CHR]],QTL[[CHR]]$bp,sp)
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
						return(  1+as.numeric(IND[LOCUS,] == "A") * QTL[LOCUS,"vals_spA"]  )
					}	
					if(sum(QTL[LOCUS,"traits",drop=FALSE] == "underdominant") == 1){
						return( 1+QTL[LOCUS,"vals_spA"] *(-1)^(sum(IND[LOCUS,] == "A")) )
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
					tapply(diff(c(PAR$bp,1)),sp.origin,sum)
				}))
			}))
			X[is.na(X)]=0
			colMeans(X)  
		}))
	}

	forqs = function(n,sp.ids,deme,n.gen,n.chr,QTL, sigma,SAVE=NULL,RETURN=TRUE){
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






#sims.s0.5 = lapply(c(5,10),function(Z){forqs(n=1000, sp.ids = rep(c("A","B"),each=500),deme = rep(1:20,each=50),n.chr=1,n.gen=Z, QTL = qtl ,mig = MIG)})


qtl=list(data.frame(traits=c("underdominant"),vals_spA = c(0.1),bp=c(0.5)))
names(qtl) = c("chr1")

sim.demes <- function (nind, ndeme, ngens,sigma) {
	nind <- ceiling(nind/ndeme)*ndeme
	forqs(n=nind, sp.ids = rep(c("A","B"),each=nind/2),deme = rep(1:ndeme,each=floor(nind/ndeme)),n.chr=1,n.gen=ngens, QTL = qtl,sigma=sigma)
}


#sigma = n*sqrt(s)/2/sqrt(tau). Define the width of the range by the number of demes. sigma/sqrt(s) is the width of the cline. sqrt(tau) is rate at which cline width is increasing (should also depend on rqtl=list(data.frame(traits=c("underdominant"),vals_spA = c(0.1),bp=c(0.5)))

sims.s0.1 = lapply(c(50),sim.demes,nind=2000,ndeme=200,sigma=2.5)

sims.sums = lapply(sims.s0.1,spBreaks)
#save(sims.sums,file="~/Projects/HybridZones/ClineProjects/sims_sums_s0.1.Robj")




loci = seq(0.5,1,0.1)

#MAKE file.snp
    info = data.frame(cbind(1,round(loci*1000000),0,loci))
    colnames(info) = NULL
    #write.csv(info,row.names=F,file="file.snp")


#GENOTYPE  [ASSUMES A ONE CHR GENOME WHERE WE HAVE PHASED DIPLOID DATA]
    genoInd = function(IND,loci){sapply(IND[[1]],function(CHR){as.numeric(CHR$sp2[as.numeric(cut(loci,c(CHR$starts,1)))])})}

#GENOTYPE ALL MY INDS
my.genos = lapply(sims.sums,function(Z){data.frame(do.call(cbind,lapply(Z$ind.ancest, genoInd, loci)))})
   #   colnames(my.genos) = NULL

#MAKE OUR FILE [NOTE THIS TAKES A LONG TIME AND IT MIGHT BE SMARTER TO WRITE THE THING ABOVE TO FILE ONE LINE AT A TIME]
   # write.table(my.genos, file = "file.geno",row.names = FALSE,sep="")


genotypes = lapply(my.genos,function(Z){(Z[,seq(1,ncol(Z),2)]+ Z[,seq(2,ncol(Z),2)])/2})

inds.per.deme = 10
freqs = lapply(genotypes,function(X){apply(X,1,function(Z){tapply(Z,cut(1:ncol(X),breaks=seq(0,ncol(X),inds.per.deme)),mean)})})

######COMPARE TO THEORY:

xx <- seq(-5,5,length.out=500)
plot(-xx, pcline(xx,r=1), type='l', xlab="distance from cline center", ylab="prob of inherited from the left" )
rr=(loci-0.5)*SIGMA/sqrt(0.1)
#rr <- seq(0,2,length.out=11)
# pp <- lapply(c(5,10,50),function(T){sapply( rr, function (r) { pcline(xx,r=r,sigma=2.5,tau=T) } )})

SIGMA = sqrt(sum(MIG*(1:length(MIG))^2))
xx <- (1:ndemes)-0.5-ndemes/2
pp <- pcline(xx,r=rr,tau=50,sigma=SIGMA)

par(mfrow=c(1,3))

for(TIME in 1:3){
matplot(-xx,pp[[TIME]],type='l',lty=1,col=rainbow(length(rr)))
matpoints(freqs[[TIME]],x=(-100:99)+0.5,lty=2,col=rainbow(length(rr)),type="l")
}

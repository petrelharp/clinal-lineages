#setwd("~/Desktop/swarmz")
#library(multicore)
# MAKE STARTING SAMPLES

source("sim-fns.R")

#PARAMS for simulation:

SIGMA = 1
ndemes = 50
deme_size = 500
ninds = ndemes*deme_size
S = 0.5


zone_age = c(50)

#PARAMS for parsing:
loci = list(seq(0.5,1,0.05))



#DEFINE selection against genotypes:
qtl=list(chr1=data.frame(traits=c("underdominant"), s = S,pos=c(0.5)))

#SIMULATE the populations (NB: Can't make per-deme population size too small, otherwise deme may end up empty and throws up an error.)
sims.s0.1 = lapply(zone_age,sim.zone,n.ind=ninds,n.deme=ndemes,sigma=SIGMA)
	#output is: a list with hierarchy: [[zone_age]][[c("inds","sp.inds","deme","pars","QTL")]]
		#where:
			#"inds" is a list of length ninds, such that [[ind]][[chr]][[X1/X2]][[data.frame where col1 = coordinates of rec., col2 = ancestor ID]]
			#"sp.inds" is a character vector of length ninds, giving the ancestry given an ancestor ID
			#"deme" is an integer vector of length ninds, giving the deme number of an ancestor, given their ancestor ID
			#"pars" is an integer that gives the number of generatsion of the simulation
			#"QTL" returns the QTL vector used in the simulation

			
#PARSE the simulated chromosomes
sims.sums = lapply(sims.s0.1,spBreaks)
	#output is:  a list: [[zone_age]][[c("ind.ancestry","ancest.prop")]]
		#"ind.ancestry" is a list of length ninds: [[ind][[chr]]][[X1/X2]][[data.frame where col1 = (chunk)starts, col2 = (chunk)stops, col3 = sp2(T/F)]]
		#"ancest.prop" is a ninds*2 matrix giving the ancestry proportion of each individual on each of their two copies of the chromosome.
		

#save(sims.sums,file="~/Projects/HybridZones/ClineProjects/sims_sums_s0.1.Robj")


######################
######################


###NOW get frequencies at each site. 

#GENOTYPE ALL MY INDS
my.genos = lapply(sims.sums,function(Z){data.frame(do.call(cbind,lapply(Z$ind.ancest, geno.ind, loci)))})
	#my.genos is a list of nloci*ninds dataframes where each entry is the number of ancestry B alleles in an individual at the given locus.
freqs = lapply(my.genos,function(X){apply(X,1,function(Z){tapply(Z,cut(1:ncol(X),breaks=seq(0,ncol(X),deme_size)),mean)/2})})

matplot(xx,freqs[[1]],col=rainbow(20),lty=1,type="l", main="50 generations",ylab="freq",xlab="deme")
legend('bottomright',col=rainbow(20),lty=1,legend=loci[[1]],cex=0.75)
matpoints(-xx,pp,col=rainbow(20),lty=3,type="l")

###########################
###########################
#To get distribution of chunks of ancestry B to the (wlog) right of the selected locus:
focal_sites = lapply(qtl,function(Z){Z$pos})

get.interval.size = function(IND_DATA=sims.sums[[1]]$ind.ancest[[1]],CHR=1,POS=0.5,ancA = TRUE){
	#return the interval containing focal site for a individual	
	chunk = sapply(IND_DATA[[CHR]],function(X){diff(as.numeric(X[which(X$starts<POS & X$stops>POS),1:2]))})	
	identity = sapply(IND_DATA[[CHR]],function(X){X[which(X$starts<POS & X$stops>POS),3]})
	replace(chunk,which(identity==ancA),0)
}


get.deme.chunks = function(IND_DATA=sims.sums[[1]]$ind.ancest, DEME = 1, CHR=1,POS=0.5,ancA = TRUE){
	#get distribution of chunks within a deme
	INDS = IND_DATA[which(sims.s0.1[[1]]$deme==DEME)]
	intervals = do.call(rbind,lapply(INDS,get.interval.size,CHR=CHR,POS=POS,ancA=ancA)) 	
	return(intervals)
}

#output is a nind*2 matrix, where each column is a chromosome. 
#EXAMPLE: testB = lapply(1:50,function(X){get.deme.chunks(DEME=X,ancA=F)})
#D = 23; hist(testB[[D]][which(testB[[D]]<1 & testB[[D]]>0)], col="black",breaks=seq(0,1,0.05))

add.alpha <- function(col, alpha=1){
if(missing(col))
stop("Please provide a vector of colours.")
apply(sapply(col, col2rgb)/255, 2, 
function(x) 
rgb(x[1], x[2], x[3], alpha=alpha)) 

}

transparent_rainbow = add.alpha(rainbow(ndemes),0.75)

D = 1; hist(testB[[D]], border="black",breaks=seq(0,1,0.05),ylim=c(0,100),xlab="length",main="Distribution of tracts")

for(D in seq(3,ndemes,5)){
	hist(testB[[D]], col=transparent_rainbow[D],border=NA,breaks=seq(0,1,0.05),ylim=c(0,50),add=T)	
}

######################
######################
#Get LD between two specified loci:

get.genotype = function(IND_DATA=sims.sums[[1]]$ind.ancest[[1]],CHR=1,POS1=0.5,POS2=0.5){
	site1 = sapply(IND_DATA[[CHR]],function(X){X[which(X$starts<POS1 & X$stops>POS1),3]})	
	site2 = sapply(IND_DATA[[CHR]],function(X){X[which(X$starts<POS2 & X$stops>POS2),3]})
	return(rbind(site1,site2))
}

get.LD = function(IND_DATA= sims.sums[[1]]$ind.ancest,DEME=1,CHR=1,POS1 = 0.5, POS2=0.1){
	INDS = IND_DATA[which(sims.s0.1[[1]]$deme==DEME)]
	genotypes = do.call(cbind,lapply(INDS,get.genotype,CHR=CHR,POS1=POS1,POS2=POS2))
	freqB1 = length(which(genotypes[1,]))/length(genotypes[1,])
	freqB2 = length(which(genotypes[2,]))/length(genotypes[2,])
	freqB12 = length(which(genotypes[1,] & genotypes[2,]))/length(genotypes[1,])
	LD = freqB12 - freqB2*freqB1
	return(LD)
}

#EXAMPLE: LD_by_deme = sapply(1:ndemes,get.LD,POS1=0.5,POS2=0.6,CHR=1,IND_DATA=sims.sums[[1]]$ind.ancest)

LD_matrix = do.call(cbind, lapply(seq(0.51,0.99,0.05),function(Z){sapply(1:ndemes,get.LD,POS1=0.5,POS2=Z,CHR=1,IND_DATA=sims.sums[[1]]$ind.ancest)}))
matplot(xx,LD_matrix,type="l",col=rainbow(20),xlim=c(-10,10),lty=1)

######COMPARE TO THEORY:

zeta <- function(x,r){
	2*pnorm(abs(x)) - 1 + 2 * exp(r*abs(x) + r^2/2 + pnorm(abs(x)+r,lower.tail=FALSE,log.p=TRUE))	
}
#zeta 

pcline <- function(x,r,tau=zone_age,sigma=SIGMA){
	(1/2)*(1-sign(x)*zeta(x=x/sqrt(tau*sigma^2),r=r*tau))	
}
#

#xx <- seq(-5,5,length.out=500)
#plot(-xx, pcline(xx,r=1,sigma=SIGMA), type='l', xlab="distance from cline center", ylab="prob of inherited from the left" )
rr=(loci[[1]]-0.5)*SIGMA/sqrt(S)
#rr <- seq(0,2,length.out=11)
# pp <- lapply(c(5,10,50),function(T){sapply( rr, function (r) { pcline(xx,r=r,sigma=2.5,tau=T) } )})

xx <- (1:ndemes)-0.5-ndemes/2
pp <- do.call(cbind,lapply(rr,function(R){pcline(xx,r=R,tau=zone_age,sigma=SIGMA)}))

par(mfrow=c(1,3))

for(TIME in 1:3){
matplot(-xx,pp[[TIME]],type='l',lty=1,col=rainbow(length(rr)))
matpoints(freqs[[TIME]],x=(-100:99)+0.5,lty=2,col=rainbow(length(rr)),type="l")
}

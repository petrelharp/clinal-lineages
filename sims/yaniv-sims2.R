#setwd("~/Desktop/swarmz")
#library(multicore)
# MAKE STARTING SAMPLES

source("sim-fns.R")

add.alpha <- function(col, alpha=1){
if(missing(col))
stop("Please provide a vector of colours.")
apply(sapply(col, col2rgb)/255, 2, 
function(x) 
rgb(x[1], x[2], x[3], alpha=alpha)) 

}


#PARAMS for simulation:

SIGMA = 1
ndemes = 50
deme_size = 500
ninds = ndemes*deme_size
S = 0.1
date = Sys.Date()

transparent_rainbow = adjustcolor(rainbow(ndemes),0.7)

zone_age = c(1000)

outfile = sprintf("simulation_SIGMA%s_Ninds%s_ndemes%s_s%s_tau%s_run%s",SIGMA,ninds,ndemes,S,zone_age,date)

#PARAMS for parsing:
loci = list(seq(0.5,1,0.05),seq(0.5,1,0.05))
xx <- (1:ndemes)-0.5-ndemes/2



#DEFINE selection against genotypes:
qtl=list(chr1=data.frame(traits=c("underdominant"), s = S,pos=c(0.5)),chr2=data.frame(traits=c("underdominant"), s = 0,pos=c(0.5)))

#SIMULATE the populations (NB: Can't make per-deme population size too small, otherwise deme may end up empty and throws up an error.)
sims.s0.1 = lapply(zone_age,sim.zone,n.ind=ninds,n.deme=ndemes,sigma=SIGMA)
	#output is: a list with hierarchy: [[zone_age]][[c("inds","sp.inds","deme","pars","QTL")]]
		#where:
			#"inds" is a list of length ninds, such that [[ind]][[chr]][[X1/X2]][[data.frame where col1 = coordinates of rec., col2 = ancestor ID]]
			#"sp.inds" is a character vector of length ninds, giving the ancestry given an ancestor ID
			#"deme" is an integer vector of length ninds, giving the deme number of an ancestor, given their ancestor ID
			#"pars" is an integer that gives the number of generatsion of the simulation
			#"QTL" returns the QTL vector used in the simulation

save(sims.s0.1,file=paste(outfile,"_rawoutput.Robj",sep=""))
			
#PARSE the simulated chromosomes
sims.sums = lapply(sims.s0.1,spBreaks)
	#output is:  a list: [[zone_age]][[c("ind.ancestry","ancest.prop")]]
		#"ind.ancestry" is a list of length ninds: [[ind][[chr]]][[X1/X2]][[data.frame where col1 = (chunk)starts, col2 = (chunk)stops, col3 = sp2(T/F)]]
		#"ancest.prop" is a ninds*2 matrix giving the ancestry proportion of each individual on each of their two copies of the chromosome.
		

save(sims.sums,file=paste(outfile,"_simsums.Robj",sep=""))


######################
######################


###NOW get frequencies at each site. 

#GENOTYPE ALL MY INDS
my.genos.list = lapply(sims.sums,function(Z){lapply(Z$ind.ancest, geno.ind, loci)})
my.genos = lapply(my.genos.list,function(AGE){lapply(1:length(loci),function(CHR){data.frame(do.call(cbind,lapply(AGE,function(IND){IND[[CHR]]})))})})
	#my.genos is a list of nloci*ninds dataframes where each entry is the number of ancestry B alleles in an individual at the given locus. [age][chr][data.frame]
freqs = lapply(my.genos,function(X){lapply(X,function(CHR){apply(CHR,1,function(Z){tapply(Z,cut(1:ncol(CHR),breaks=seq(0,ncol(CHR),deme_size)),mean)/2})})})


#pdf(file="freqplot.pdf")
matplot(xx,freqs[[1]][[1]],col=rainbow(20),lty=1,type="l", main=paste(zone_age,"generations"),ylab="freq",xlab="deme")
matpoints(xx,freqs[[1]][[2]],col="grey",lty=1,type="l", main=paste(zone_age,"generations","neutral"),ylab="freq",xlab="deme")

legend('bottomright',col=rainbow(20),lty=1,legend=loci[[1]],cex=0.75)

#dev.off()

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
#EXAMPLE: testB = lapply(1:50,function(X){get.deme.chunks(DEME=X,ancA=FALSE)})
#D = 23; hist(testB[[D]][which(testB[[D]]<1 & testB[[D]]>0)], col="black",breaks=seq(0,1,0.05))

testB = lapply(1:ndemes,function(X){get.deme.chunks(DEME=X,ancA=FALSE)})
testB_far = lapply(1:ndemes,function(X){get.deme.chunks(DEME=X,ancA=FALSE,POS=0.01)})
testB_unlinked = lapply(1:ndemes,function(X){get.deme.chunks(DEME=X,ancA=FALSE,POS=0.5,CHR=2)})

chunks = list(selected=testB,distant=testB_far,unlinked=testB_unlinked)
save(chunks,file = paste(outfile,"_simsums_chunks.Robj",sep=""))

empty_deme = rep(0,2*deme_size)

chunks_ecdf = lapply(chunks,function(X){
	lapply(X,function(D){
		segregating =which(D<1 & D>0)
		if(length(segregating)>0){
		return(ecdf(D[segregating]))}else return(ecdf(D))	
	})})
	
pdf(file = paste(outfile,"_chunks_ecdf.pdf",sep=""))
for(type in names(chunks_ecdf)){
	plot(1, type="n",ylim=c(0,1),xlim=c(0,0.2),main=type,ylab="1-ecdf(x)",xlab="x=chunk length")	
	legend("topright",legend=seq(-24.5,24.5,5),col=transparent_rainbow[seq(1,50,5)],lty=1)
	for(i in rev(1:ndemes)){points(1-chunks_ecdf[[type]][[i]](seq(0,1,0.005))~seq(0,1,0.005),ylim=c(0,1),type="l",col=transparent_rainbow[i])}
}
dev.off()

pdf(file = paste(outfile,"_chunks_density.pdf",sep=""))
for(type in names(chunks)){
	D = 1; hist(chunks[[type]][[1]], border="white",breaks=seq(0,1,0.05),ylim=c(0,200),xlab="length",main=sprintf("Distribution of %s tracts",type),xlim=c(0,0.2))
	for(D in rev(seq(3,ndemes,1))){
		relevant_chunks = chunks[[type]][[D]][which(chunks[[type]][[D]]>0 & chunks[[type]][[D]]<1)]
	#hist(relevant_chunks, col=transparent_rainbow[D],border=NA,breaks=seq(0,1,0.001),ylim=c(0,50),add=T)	}
	if(length(relevant_chunks)>1){points(density(relevant_chunks), col=transparent_rainbow[D],border=NA,breaks=seq(0,1,0.001),type="l")}}
}
dev.off()

#pdf(file = paste(outfile,"_selected_chunks.pdf",sep=""))

#D = 1; hist(chunks[[1]][[1]], border="white",breaks=seq(0,1,0.05),ylim=c(0,200),xlab="length",main="Distribution of selected tracts",xlim=c(0,0.2))
#for(D in rev(seq(3,ndemes,5))){
#	hist(chunks[[1]][[D]][which(chunks[[1]][[D]]>0 & chunks[[1]][[D]]<1)], col=transparent_rainbow[D],border=NA,breaks=seq(0,1,0.01),ylim=c(0,50),add=T)	}
#dev.off()

#pdf(file = paste(outfile,"_faraway_chunks.pdf",sep=""))
#D = 1; hist(chunks[[2]][[1]], border="white",breaks=seq(0,1,0.05),ylim=c(0,100),xlab="length",main="Distribution of distant tracts")
#for(D in rev(seq(3,ndemes,5))){
#	hist(chunks[[2]][[D]][which(chunks[[2]][[D]]>0 & chunks[[2]][[D]]<1)], col=transparent_rainbow[D],border=NA,breaks=seq(0,1,0.01),ylim=c(0,50),add=T)	}
#dev.off()

#pdf(file = paste(outfile,"_unlinked_chunks.pdf",sep=""))
#D = 1; hist(chunks[[3]][[D]], border="white",breaks=seq(0,1,0.05),ylim=c(0,100),xlab="length",main="Distribution of unlinked tracts")
#for(D in rev(seq(3,ndemes,5))){
#	hist(chunks[[3]][[D]][which(chunks[[3]][[D]]>0 & chunks[[3]][[D]]<1)], col=transparent_rainbow[D],border=NA,breaks=seq(0,1,0.01),ylim=c(0,50),add=T)	}
#dev.off()

harmmean_sel=1/sapply(testB,function(Z){mean(sapply(Z[which(Z>0 & Z<1)],function(X){1/X}))})
harmmean_unlinked=1/sapply(testB_unlinked,function(Z){mean(sapply(Z[which(Z>0&Z<1)],function(X){1/X}))})
harmmean_far = 1/sapply(testB_far,function(Z){mean(sapply(Z[which(Z>0&Z<1)],function(X){1/X}))})


pdf(file = paste(outfile,"_chunks_hmean.pdf",sep=""))

plot(harmmean_unlinked,main="harmonic mean",ylim=c(0,1),xlab="deme",ylab="harmonic mean")
points(harmmean_sel,col="red")
points(harmmean_far,col="blue")

legend('topleft',legend=c("unlinked","sel.","far."),col=c("black","red","blue"),pch=1)
dev.off()


geomean_unlinked=sapply(testB_unlinked,function(Z){prod(Z+1)^(1/length(Z))})
geomean_sel=sapply(testB,function(Z){prod(Z+1)^(1/length(Z))})
geomean_far=sapply(testB_far,function(Z){prod(Z+1)^(1/length(Z))})

pdf(file = paste(outfile,"_chunks_gmean.pdf",sep=""))

plot(geomean_unlinked-1,main="geometric mean",ylim=c(0,1),xlab="deme")
points(geomean_sel-1,col="red")
points(geomean_far-1,col="blue")
legend('topleft',legend=c("unlinked","sel.","far."),col=c("black","red","blue"),pch=1)

dev.off()
#Get distribution of chunk length around selected locus. 

pdf(file = paste(outfile,"_chunks_mean.pdf",sep=""))
mean_unlinked=sapply(testB_unlinked,mean)
mean_sel=sapply(testB,mean)
mean_far =sapply(testB_far,mean)

plot(mean_unlinked,main="mean",ylim=c(0,1),xlab="deme")
points(mean_sel,col="red")
points(mean_far,col="blue")

legend('topleft',legend=c("neu.","sel.","far."),col=c("black","red","blue"),pch=1)
dev.off()

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

pdf(file = paste(outfile,"_LD.pdf",sep=""))
matplot(xx,LD_matrix,type="l",col=rainbow(20),xlim=c(-10,10),lty=1)
dev.off()

######COMPARE TO THEORY:

#zeta <- function(x,r){
#2*pnorm(abs(x)) - 1 + 2 * exp(r*abs(x) + r^2/2 + pnorm(abs(x)+r,lower.tail=FALSE,log.p=TRUE))	
#}
##zeta 
#
#pcline <- function(x,r,tau=zone_age,sigma=SIGMA){
#(1/2)*(1-sign(x)*zeta(x=x/sqrt(tau*sigma^2),r=r*tau))	
#}
##
#
##xx <- seq(-5,5,length.out=500)
##plot(-xx, pcline(xx,r=1,sigma=SIGMA), type='l', xlab="distance from cline center", ylab="prob of inherited from the left" )
#rr=(loci[[1]]-0.5)*SIGMA/sqrt(S)
##rr <- seq(0,2,length.out=11)
## pp <- lapply(c(5,10,50),function(T){sapply( rr, function (r) { pcline(xx,r=r,sigma=2.5,tau=T) } )})
#
#xx <- (1:ndemes)-0.5-ndemes/2
#pp <- do.call(cbind,lapply(rr,function(R){pcline(xx,r=R,tau=zone_age,sigma=SIGMA)}))
#
#par(mfrow=c(1,3))
#
#for(TIME in 1:3){
#matplot(-xx,pp[[TIME]],type='l',lty=1,col=rainbow(length(rr)))
#matpoints(freqs[[TIME]],x=(-100:99)+0.5,lty=2,col=rainbow(length(rr)),type="l")
#}

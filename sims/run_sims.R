#!/usr/bin/Rscript

source("sim-fns.R")

#PARAMS for simulation:

params = list(
    SIGMA = 3,
    ndemes = 100,
    deme_size = 50,
    S = 0.05,
    zone_age = 10
)
params$ninds = params$ndemes * params$deme_size

run.id = sprintf("%0.0f",1e6*runif(1))  # unique run id
outfile = with( list2env(params), sprintf("simulation_SIGMA%s_Ninds%s_ndemes%s_s%s_dir/tau%s_dir/results_runid_%s",SIGMA,ninds,ndemes,S,zone_age,run.id) )
dir.create(dirname(outfile),recursive=TRUE)

#PARAMS for parsing:
loci = list(seq(0.5,1,0.05),seq(0.5,1,0.05))
xx = (1:params$ndemes)-0.5-params$ndemes/2


#DEFINE selection against genotypes:
qtl=list(
        chr1=data.frame(
                traits=c("underdominant"), 
                s=params$S,
                pos=c(0.5)),
        chr2=data.frame(
                traits=c("underdominant"), 
                s = 0,pos=c(0.5))
        )

# save out parameters
save( params, outfile, loci, xx, qtl, file=paste0(outfile,"_params.Robj") )

#SIMULATE the populations
sims.full = lapply(params$zone_age, sim.zone, n.ind=params$ninds, n.deme=params$ndemes, sigma=params$SIGMA)
	#output is: a list with hierarchy: [[zone_age]][[c("inds","sp.inds","deme","pars","QTL")]]
		#where:
			#"inds" is a list of length ninds, such that [[ind]][[chr]][[X1/X2]][[data.frame where col1 = coordinates of rec., col2 = ancestor ID]]
			#"sp.inds" is a character vector of length ninds, giving the ancestry given an ancestor ID
			#"deme" is an integer vector of length ninds, giving the deme number of an ancestor, given their ancestor ID
			#"pars" is an integer that gives the number of generatsion of the simulation
			#"QTL" returns the QTL vector used in the simulation

save(sims.full,file=paste(outfile,"_rawoutput.Robj",sep=""))
			
#PARSE the simulated chromosomes
sims.sums = lapply(sims.full,spBreaks)
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
freqs = lapply(my.genos,function(X){lapply(X,function(CHR){apply(CHR,1,function(Z){tapply(Z,cut(1:ncol(CHR),breaks=seq(0,ncol(CHR),params$deme_size)),mean)/2})})})


pdf(file=paste0(outfile,"_freqplot.pdf"))
matplot(xx,freqs[[1]][[1]],col=rainbow(20),lty=1,type="l", main=paste(params$zone_age,"generations"),ylab="freq",xlab="deme")
matpoints(xx,freqs[[1]][[2]],col="grey",lty=1,type="l", main=paste(params$zone_age,"generations","neutral"),ylab="freq",xlab="deme")

legend('bottomright',col=rainbow(20),lty=1,legend=loci[[1]],cex=0.75)

dev.off()

###########################
###########################
# To get distribution of chunks of ancestry B surrounding the chosen locus:


focal_sites = lapply(qtl,function(Z){Z$pos})

get.interval.size = function(IND_DATA=sims.sums[[1]]$ind.ancest[[1]],CHR=1,POS=0.5,ancA = TRUE){
	#return the interval containing focal site for a individual	
	chunk = sapply(IND_DATA[[CHR]],function(X){diff(as.numeric(X[which(X$starts<POS & X$stops>POS),1:2]))})	
	identity = sapply(IND_DATA[[CHR]],function(X){X[which(X$starts<POS & X$stops>POS),3]})
	replace(chunk,which(identity==ancA),0)
}


get.deme.chunks = function(IND_DATA=sims.sums[[1]]$ind.ancest, DEME = 1, CHR=1,POS=0.5,ancA = TRUE){
	#get distribution of chunks within a deme
	INDS = IND_DATA[which(sims.full[[1]]$deme==DEME)]
	intervals = do.call(rbind,lapply(INDS,get.interval.size,CHR=CHR,POS=POS,ancA=ancA)) 	
	return(intervals)
}

#output is a nind*2 matrix, where each column is a chromosome. 
#EXAMPLE: testB = lapply(1:50,function(X){get.deme.chunks(DEME=X,ancA=FALSE)})
#D = 23; hist(testB[[D]][which(testB[[D]]<1 & testB[[D]]>0)], col="black",breaks=seq(0,1,0.05))

chunks.rvals = list( 
        selected=c(chr=1,pos=0.5),
        near0.01=c(chr=1,pos=0.51),
        near0.05=c(chr=1,pos=0.55),
        near0.1=c(chr=1,pos=0.6),
        near0.2=c(chr=1,pos=0.7),
        near0.4=c(chr=1,pos=0.9),
        unlinked=c(chr=2,pos=0.5),
        unlinked0.4=c(chr=2,pos=0.9) )
chunks <- lapply( chunks.rvals, function (xx) {
            lapply(1:params$ndemes, function (X) { get.deme.chunks(DEME=X,ancA=FALSE,POS=xx[2],CHR=xx[1]) } )
        } )

# testB = lapply(1:params$ndemes,function(X){get.deme.chunks(DEME=X,ancA=FALSE)})
# testB_far = lapply(1:params$ndemes,function(X){get.deme.chunks(DEME=X,ancA=FALSE,POS=0.01)})
# testB_unlinked = lapply(1:params$ndemes,function(X){get.deme.chunks(DEME=X,ancA=FALSE,POS=0.5,CHR=2)})
# chunks = list(selected=testB,distant=testB_far,unlinked=testB_unlinked)

save(chunks, chunks.rvals, file = paste(outfile,"_simsums_chunks.Robj",sep=""))

empty_deme = rep(0,2*params$deme_size)

transparent_rainbow = adjustcolor(rainbow(params$ndemes),0.7)

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
	for(i in rev(1:params$ndemes)){points(1-chunks_ecdf[[type]][[i]](seq(0,1,0.005))~seq(0,1,0.005),ylim=c(0,1),type="l",col=transparent_rainbow[i])}
}
dev.off()

pdf(file = paste(outfile,"_chunks_density.pdf",sep=""))
for(type in names(chunks)){
	D = 1; hist(chunks[[type]][[1]], border="white",breaks=seq(0,1,0.05),ylim=c(0,200),xlab="length",main=sprintf("Distribution of %s tracts",type),xlim=c(0,0.2))
	for(D in rev(seq(3,params$ndemes,1))){
		relevant_chunks = chunks[[type]][[D]][which(chunks[[type]][[D]]>0 & chunks[[type]][[D]]<1)]
	#hist(relevant_chunks, col=transparent_rainbow[D],border=NA,breaks=seq(0,1,0.001),ylim=c(0,50),add=T)	}
	if(length(relevant_chunks)>1){points(density(relevant_chunks), col=transparent_rainbow[D],border=NA,breaks=seq(0,1,0.001),type="l")}}
}
dev.off()

#Get distribution of chunk length around selected locus. 

pdf(file = paste(outfile,"_chunks_mean.pdf",sep=""))

mean_lens = sapply( chunks, sapply, mean )
mean_cond_lens = sapply( chunks, sapply, function (x) { if (any(x>0)) { mean(x[x>0]) } else { NA } } )

matplot(mean_lens, type='l', lty=sapply(chunks.rvals,"[",1),
        main=outfile,
        xlab="spatial position", ylab="mean B haplotype length (including zeros)")
legend("topleft", lty=sapply(chunks.rvals,"[",1), col=1:6,
        legend=sapply( chunks.rvals, function (x) { sprintf("chr %d: %0.3f",x[1],x[2]) } ) )

matplot(mean_cond_lens, type='l', lty=sapply(chunks.rvals,"[",1),
        main=outfile,
        xlab="spatial position", ylab="mean conditional B haplotype length")
legend("topleft", lty=sapply(chunks.rvals,"[",1), col=1:6,
        legend=sapply( chunks.rvals, function (x) { sprintf("chr %d: %0.3f",x[1],x[2]) } ) )

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
	INDS = IND_DATA[which(sims.full[[1]]$deme==DEME)]
	genotypes = do.call(cbind,lapply(INDS,get.genotype,CHR=CHR,POS1=POS1,POS2=POS2))
	freqB1 = length(which(genotypes[1,]))/length(genotypes[1,])
	freqB2 = length(which(genotypes[2,]))/length(genotypes[2,])
	freqB12 = length(which(genotypes[1,] & genotypes[2,]))/length(genotypes[1,])
	LD = freqB12 - freqB2*freqB1
	return(LD)
}

#EXAMPLE: LD_by_deme = sapply(1:params$ndemes,get.LD,POS1=0.5,POS2=0.6,CHR=1,IND_DATA=sims.sums[[1]]$ind.ancest)

LD_matrix = do.call(cbind, lapply(seq(0.51,0.99,0.05),function(Z){sapply(1:params$ndemes,get.LD,POS1=0.5,POS2=Z,CHR=1,IND_DATA=sims.sums[[1]]$ind.ancest)}))

pdf(file = paste(outfile,"_LD.pdf",sep=""))
matplot(xx,LD_matrix,type="l",col=rainbow(20),xlim=c(-10,10),lty=1)
dev.off()




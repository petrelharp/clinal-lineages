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

#########
#########
#########
## stuff from potential_statistics: haplotype lengths, etcetera.

chromosome = "chr1"
plot_start = 0.48
plot_stop=0.52
increments = 0.001

positions = seq(plot_start,plot_stop,increments)

outstring = with( list2env(params), sprintf("simulation_SIGMA%s_Ninds%s_ndemes%s_s%s_dir/tau%s_dir/results_runid_%s_chr%s_start%g_stop%g_by%g",
                                            SIGMA, ninds, ndemes, S, zone_age, run.id,  
                                            chromosome, plot_start, plot_stop, increments ) )

deme_ID = rep( 1:params$ndemes, each=params$deme_size )


intervalSizes = lapply(positions,function(POS){
                do.call(rbind,lapply(sims.sums[[1]]$ind.ancest,function(IND){
                    get.interval.size(IND_DATA=IND,CHR=chromosome,POS=POS,ancA=TRUE)	
                }))
            })

intervalSizes_allAncs = lapply(positions,function(POS){
            do.call(rbind,lapply(sims.sums[[1]]$ind.ancest,function(IND){
				get.interval.size(IND_DATA=IND,CHR=chromosome,POS=POS,ancA=TRUE,restrict.anc=FALSE)
            }))
        })

#[[POS]][ind,chr]

# from deSolve:::drawlegend
.drawlegend <- function ( zlim=c(0,1), colfn=heat.colors, horizontal=TRUE ) {
    iy <- 1
    minz <- zlim[1]
    maxz <- zlim[2]
    binwidth <- (maxz - minz)/64
    ix <- seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
    iz <- matrix(ix, nrow = length(ix), ncol = length(iy))
    if (horizontal) {
        opar <- par(mar=c(par("mar"),4,1)[c(6,2,5,4)])
        on.exit(par(opar),add=TRUE)
        image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = colfn(64))
        do.call("axis", list(side = 3, mgp = c(3, 1, 0), las = 1))
    } else {
        opar <- par(mar=c(par("mar"),3,0.5)[c(1,5,3,6)])
        on.exit(par(opar),add=TRUE)
        image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = colfn(64))
        do.call("axis", list(side = 2, mgp = c(3, 1, 0), las = 2))
    }
}


##
our_image <- function (z, main='', breaks=pretty(z,n=25), col=heat.colors(length(breaks)-1)) {
    # do an image plot with a scale bar
    layout(c(1,2),heights=c(1,4))
    on.exit(layout(1),add=TRUE)
    opar <- par(mar=c(0.5,3.5,4,1))
    on.exit(par(opar),add=TRUE)
    # do legend
        .drawlegend(zlim=range(breaks))
        mtext(main, side=3, line=2.5, font=2)
    # do plot
        par(mar=c(3.5,3.5,0,1),mgp=c(2.5,1,0))
        image(x=(positions-0.5)*100, y=xx, z=z, main="", breaks=breaks, col=col,
            xlab="Distance from selected locus (cM)",
            ylab="Distance from zone center" )
}

.spatial.legend <- function () {
    # do the legend for a matplot where one line correponds to a location
    these.lines <- unique(floor(seq(1,length(xx),length.out=10)))
    legend("topright",
           title="distance from center:",
           legend = xx[these.lines],
           cex=0.7,col=rainbow(length(xx))[these.lines],lty=1,xjust=1)
}

## from make_length_along_chr.R

pdf(height=4, width=6.25,file=paste(outstring,"blocksAlongChromNoConditioning.pdf",sep="_"))
par(mar=c(3.5,3.5,0.5,0.5))
mean_per_deme = do.call(rbind, lapply(intervalSizes_allAncs,function(P){tapply(1:length(deme_ID),deme_ID,function(Z){mean(P[Z,])})}))
plot(mean_per_deme[,1]/mean(mean_per_deme[,1],na.rm=T),x=positions,ylim=c(0.5,3),type="l",lty=1,col="white",main="",cex.main=0.8,xlab="",ylab="",yaxt="n",xaxt="n")
for(i in 1:50){points(mean_per_deme[,i]/mean(mean_per_deme[,i],na.rm=T),x=positions,type="l",col=rainbow(50)[i])}
axis(1,at=seq(0.48,0.52,0.01),labels=seq(-0.02,0.02,0.01),cex.axis=0.9)
axis(2,at=seq(0.5,3,0.5),cex.axis=0.9,las=2)
mtext("Distance from selected locus (M)",side=1,line=2.5)
mtext("Normalized mean block length",side=2,line=2.5)
.spatial.legend()
dev.off()

pdf(height=4, width=6.25,file=paste(outstring,"blocksAlongChromAncBConditioning.pdf",sep="_"))
par(mar=c(3.5,3.5,0.5,0.5))
mean_per_deme_AncB = do.call(rbind, lapply(intervalSizes,function(P){tapply(1:length(deme_ID),deme_ID,function(Z){all_sites = P[Z,]; mean(all_sites[which(all_sites>0)])})}))
plot(mean_per_deme_AncB[,1]/mean(mean_per_deme_AncB[,1],na.rm=T),x=positions,ylim=c(0.5,3),type="l",lty=1,col="white",main="",cex.main=0.8,xlab="",ylab="",yaxt="n",xaxt="n")
for(i in 1:50){points(mean_per_deme_AncB[,i]/mean(mean_per_deme_AncB[,i],na.rm=T),x=positions,type="l",col=rainbow(50)[i])}
axis(1,at=seq(0.48,0.52,0.01),labels=seq(-0.02,0.02,0.01),cex.axis=0.9)
axis(2,at=seq(0.5,3,0.5),cex.axis=0.9,las=2)
mtext("Distance from selected locus (M)",side=1,line=2.5)
mtext("Normalized mean block length",side=2,line=2.5)
.spatial.legend()
dev.off()

pdf(height=4, width=6.25,file=paste(outstring,"blocksAlongChromHeatmap.pdf",sep="_"))
our_image(mean_per_deme)
dev.off()

pdf(height=4, width=6.25,file=paste(outstring,"blocksAlongChromHeatmapAncBConditioning.pdf",sep="_"))
our_image(mean_per_deme_AncB)
dev.off()

pdf(height=4, width=6.25,file=paste(outstring,"adjacentBlocksAlongChromNoConditioning.pdf",sep="_"))
par(mar=c(3.5,3.5,0.5,0.5))
flanking.blocks.deme.matrix = do.call(rbind,lapply(1:params$ndemes,function(DEME){sapply(flanking.blocks.by.ind,function(POS){mean(POS[[DEME]])})}))
matplot(t(flanking.blocks.deme.matrix/apply(flanking.blocks.deme.matrix,1,mean))[,params$ndemes:1],xaxt="n",yaxt="n",xlab="",ylab="",lty=1,col=rainbow(params$ndemes),type="l",x=positions,ylim=c(0.5,3))
axis(1,at=seq(0.48,0.52,0.01),labels=seq(-0.02,0.02,0.01),cex.axis=0.9)
axis(2,at=seq(0.5,3,0.5),cex.axis=0.9,las=2)
mtext("Distance from selected locus (M)",side=1,line=2.5)
mtext("Normalized mean block length",side=2,line=2.5)
.spatial.legend()
dev.off()

## flanking blocks
flanking.blocks.deme.matrix = do.call(rbind,lapply(1:params$ndemes,function(DEME){sapply(flanking.blocks.by.ind,function(POS){mean(POS[[DEME]])})}))

#For AncB inds:
get.flanking.blocks.AncB = function(IND_DATA=sims.sums[[1]]$ind.ancest[[1]],CHR=chromosome,POS=0.5,ancB=TRUE){
	focal_chunks = do.call(rbind,lapply(IND_DATA[[CHR]],function(X){
		FOCUS = which(X$starts<=POS & X$stops>=POS);
		if(X[FOCUS,]$sp2){
		if(max(FOCUS-1,1)!=min(FOCUS+1,nrow(X))){return(X[c(max(FOCUS-1,1),min(FOCUS+1,nrow(X))),])}else return(X[max(FOCUS-1,1),])}else return()}
		))
	focal_chunks$stops-focal_chunks$starts
	}

flanking.blocks.by.ind.AncB = lapply(positions,function(POS){
	tapply(1:length(deme_ID),deme_ID,function(DEME){do.call(c,lapply(DEME,function(IND){
			get.flanking.blocks.AncB(IND_DATA=sims.sums[[1]]$ind.ancest[[IND]],POS=POS)}))})})

flanking.blocks.deme.matrix.AncB = do.call(rbind,lapply(1:params$ndemes,function(DEME){sapply(flanking.blocks.by.ind.AncB,function(POS){mean(POS[[DEME]])})}))


pdf(height=4, width=6.25,file=paste(outstring,"adjacentBlocksAlongChromAncBConditioning.pdf",sep="_"))
par(mar=c(3.5,3.5,0.5,0.5))
plot(flanking.blocks.deme.matrix.AncB[1,]/mean(flanking.blocks.deme.matrix.AncB[1,],na.rm=T),
     x=positions,ylim=c(0,3),type="l",lty=1,col="white",main="",cex.main=0.8,xlab="",ylab="",yaxt="n",xaxt="n")
for(i in params$ndemes:1){points(flanking.blocks.deme.matrix.AncB[i,]/mean(flanking.blocks.deme.matrix.AncB[i,],na.rm=T),x=positions,type="l",col=rainbow(params$ndemes)[params$ndemes-i])}
axis(1,at=seq(0.48,0.52,0.01),labels=seq(-0.02,0.02,0.01),cex.axis=0.9)
axis(2,at=seq(0.5,3,0.5),cex.axis=0.9,las=2)
mtext("Distance from selected locus (M)",side=1,line=2.5)
mtext("Normalized mean block length",side=2,line=2.5)
.spatial.legend()
dev.off()

pdf(height=4, width=6.25,file=paste(outstring,"ratioAdjacentBlocksAlongChromNoConditioning.pdf",sep="_"))
par(mar=c(3.5,3.5,0.5,0.5))
ratio.blocks = mean_per_deme/t(flanking.blocks.deme.matrix[params$ndemes:1,])
plot(ratio.blocks[,1]/mean(ratio.blocks[,1],na.rm=T),x=positions,ylim=c(0.5,3),type="l",lty=1,col="white",main="",cex.main=0.8,xlab="",ylab="",yaxt="n",xaxt="n")
for(i in 1:params$ndemes){points(ratio.blocks[,i]/mean(ratio.blocks[,i],na.rm=T),x=positions,type="l",col=rainbow(params$ndemes)[i])}
axis(1,at=seq(0.48,0.52,0.01),labels=seq(-0.02,0.02,0.01),cex.axis=0.9)
axis(2,at=seq(0.5,3,0.5),cex.axis=0.9,las=2)

mtext("Distance from selected locus (M)",side=1,line=2.5)
mtext("Normalized mean block length",side=2,line=2.5)
.spatial.legend()
dev.off()



pdf(height=4, width=6.25,file=paste(outstring,"ratioAdjacentBlocksAlongChromHeatmapNoConditioning.pdf",sep="_"))
our_image(ratio.blocks[,params$ndemes:1])
dev.off()


pdf(height=4, width=6.25,file=paste(outstring,"ratioAdjacentBlocksAlongChromAncBConditioning.pdf",sep="_"))
par(mar=c(3.5,3.5,0.5,0.5))
ratio.blocks.AncB = mean_per_deme_AncB/t(flanking.blocks.deme.matrix.AncB[params$ndemes:1,])
plot(ratio.blocks.AncB[,1]/mean(ratio.blocks.AncB[,1],na.rm=T),x=positions,ylim=c(0.0,3),type="l",lty=1,col="white",main="",cex.main=0.8,xlab="",ylab="",yaxt="n",xaxt="n")
for(i in 1:params$ndemes){points(ratio.blocks.AncB[,i]/mean(ratio.blocks.AncB[,i],na.rm=T),x=positions,type="l",col=rainbow(params$ndemes)[i])}
axis(1,at=seq(0.48,0.52,0.01),labels=seq(-0.02,0.02,0.01),cex.axis=0.9)
axis(2,at=seq(0.5,3,0.5),cex.axis=0.9,las=2)

mtext("Distance from selected locus (M)",side=1,line=2.5)
mtext("Normalized mean block length",side=2,line=2.5)
.spatial.legend()
dev.off()

pdf(height=4, width=6.25,file=paste(outstring,"ratioAdjacentBlocksAlongChromHeatmapAncBConditioning.pdf",sep="_"))
our_image(ratio.blocks.AncB[,params$ndemes:1])
dev.off()


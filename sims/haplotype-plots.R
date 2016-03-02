# needs to have defined:
#   params
# and have loaded:
#   sims.sums
# and source sim-fns.R

# this could be useful for that:
get_simparams <- function (fname) {
    strings <- c( SIGMA="SIGMA", ninds="Ninds", ndemes="ndemes", S="s", tau="tau", run.id="runid_" )
    out <- lapply( strings, function (x) { as.numeric(gsub(sprintf(".*[/_]%s([0-9.]+)_.*",x),"\\1",fname)) } )
    out$zone_age <- out$tau
    out$deme_size <- out$ninds/out$ndemes
    return(out)
}

if (!exists("xx")) { xx = (1:params$ndemes)-0.5-params$ndemes/2 }


#########
#########
#########
## stuff from make_length_along_chr.R

chromosome = "chr1"
plot_start = 0.48
plot_stop=0.52
increments = 0.001

positions = seq(plot_start,plot_stop,increments)

outstring = with( list2env(params), sprintf("simulation_SIGMA%s_Ninds%s_ndemes%s_s%s_dir/tau%s_dir/results_runid_%s_%s_start%g_stop%g_by%g",
                                            SIGMA, ninds, ndemes, S, zone_age, run.id,  
                                            chromosome, plot_start, plot_stop, increments ) )

deme_ID = rep( 1:params$ndemes, each=params$deme_size )


# returns a list, corresponding to positions along the genome,
# with the [[k]]th element a (ninds x ploidy) matrix, one column per chromosome,
# giving the lengths of segments of B ancestry at the k-th position (length is zero if ancestry is A)
if (!file.exists(paste(outstring,"intervalSizes.Robj",sep="_"))) {
    intervalSizes = lapply(positions,function(POS){
                    do.call(rbind,lapply(sims.sums[[1]]$ind.ancest,function(IND){
                        get.interval.size(IND_DATA=IND,CHR=chromosome,POS=POS,ancA=TRUE)	
                    }))
                })
    save(intervalSizes, file=paste(outstring,"intervalSizes.Robj",sep="_"))
} else {
    load(paste(outstring,"intervalSizes.Robj",sep="_"))
}


# same, but not zero-ing out A ancestry
if (!file.exists(paste(outstring,"intervalSizes_allAncs.Robj",sep="_"))) {
    intervalSizes_allAncs = lapply(positions,function(POS){
                do.call(rbind,lapply(sims.sums[[1]]$ind.ancest,function(IND){
                    get.interval.size(IND_DATA=IND,CHR=chromosome,POS=POS,ancA=TRUE,restrict.anc=FALSE)
                }))
            })
    save(intervalSizes_allAncs, file=paste(outstring,"intervalSizes_allAncs.Robj",sep="_"))
} else {
    load(paste(outstring,"intervalSizes_allAncs.Robj",sep="_"))
}

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
    layout(c(1,2),heights=c(1.5,4))
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
paramstring <- sprintf("tau=%d, s=%0.3f",params$zone_age, params$S)

# this is a (locus x deme) matrix giving mean lengths
#  for haplotype chunks of ALL ANCESTRIES
mean_per_deme = do.call(rbind, lapply(intervalSizes_allAncs,function(P){tapply(1:length(deme_ID),deme_ID,function(Z){mean(P[Z,])})}))

# this plots mean block length, regardless of ancestry
pdf(height=4, width=6.25,file=paste(outstring,"blocksAlongChromNoConditioning_nonnormalized.pdf",sep="_"))
par(mar=c(3.5,3.5,1.5,0.5),mgp=c(2.5,1,0))
matplot(100*(positions-0.5), 100*mean_per_deme, type='l', lty=1, col=rainbow(params$ndemes),
        xlab="Distance from selected locus (cM)",
        ylab="Mean block length (cM)",
        main=paste("Mean block length",paramstring) )
.spatial.legend()
dev.off()

# this gives a heatmap of the same thing
pdf(height=4, width=6.25,file=paste(outstring,"blocksAlongChromHeatmap.pdf",sep="_"))
our_image(100*mean_per_deme,
        main=paste("Mean block length (cM)",paramstring) )
dev.off()

# this plots mean block length, regardless of ancestry
#   divided by the mean value for that deme
pdf(height=4, width=6.25,file=paste(outstring,"blocksAlongChromNoConditioning.pdf",sep="_"))
par(mar=c(3.5,3.5,1.5,0.5))
plot(mean_per_deme[,1]/mean(mean_per_deme[,1],na.rm=TRUE),x=(positions-0.5),ylim=c(0.5,3),type="l",lty=1,col="white",
     main=paste("Mean block length, normalized by deme",paramstring),
     cex.main=0.8,xlab="",ylab="",yaxt="n",xaxt="n")
for(i in 1:params$ndemes){points(mean_per_deme[,i]/mean(mean_per_deme[,i],na.rm=TRUE),x=(positions-0.5),type="l",col=rainbow(params$ndemes)[i])}
axis(1,at=seq(0.48,0.52,0.01),labels=seq(-0.02,0.02,0.01),cex.axis=0.9)
axis(2,at=seq(0.5,3,0.5),cex.axis=0.9,las=2)
mtext("Distance from selected locus (M)",side=1,line=2.5)
mtext("Normalized mean block length",side=2,line=2.5)
.spatial.legend()
dev.off()

# this is a (locus x deme) matrix giving mean lengths
#   for haplotype chunks of B ancestry
mean_per_deme_AncB = do.call(rbind, lapply(intervalSizes,function(P){tapply(1:length(deme_ID),deme_ID,function(Z){all_sites = P[Z,]; mean(all_sites[which(all_sites>0)])})}))

# this plots mean block length, for only ancestry B
pdf(height=4, width=6.25,file=paste(outstring,"blocksAlongChromAncBConditioning_nonnormalized.pdf",sep="_"))
par(mar=c(3.5,3.5,1.5,0.5),mgp=c(2.5,1,0))
matplot(100*(positions-0.5), 100*mean_per_deme_AncB, type='l', lty=1, col=rainbow(params$ndemes),
        xlab="Distance from selected locus (cM)",
        ylab="Mean block length (cM)",
        main=paste("Mean block length, ancestry B",paramstring) )
.spatial.legend()
dev.off()

# and a heatmap of the same thing
pdf(height=4, width=6.25,file=paste(outstring,"blocksAlongChromHeatmapAncBConditioning.pdf",sep="_"))
our_image(100*mean_per_deme_AncB,
        main=paste("Mean block length, ancestry B (cM)",paramstring) )
dev.off()


# this plots mean block length of ancestry B,
#   divided by the mean value for that deme
pdf(height=4, width=6.25,file=paste(outstring,"blocksAlongChromAncBConditioning.pdf",sep="_"))
par(mar=c(3.5,3.5,0.5,0.5))
plot(mean_per_deme_AncB[,1]/mean(mean_per_deme_AncB[,1],na.rm=TRUE),x=(positions-0.5),ylim=c(0.5,3),type="l",lty=1,col="white",main="",cex.main=0.8,xlab="",ylab="",yaxt="n",xaxt="n")
for(i in 1:params$ndemes){points(mean_per_deme_AncB[,i]/mean(mean_per_deme_AncB[,i],na.rm=TRUE),x=(positions-0.5),type="l",col=rainbow(params$ndemes)[i])}
axis(1,at=seq(0.48,0.52,0.01),labels=seq(-0.02,0.02,0.01),cex.axis=0.9)
axis(2,at=seq(0.5,3,0.5),cex.axis=0.9,las=2)
mtext("Distance from selected locus (M)",side=1,line=2.5)
mtext("Normalized mean block length",side=2,line=2.5)
.spatial.legend()
dev.off()


####
# now on to neighboring blocks


if (!file.exists(paste(outstring,"flankingBlocks.Robj",sep="_"))) {
    # this is a list, corresponding to position along the genome,
    #  whose [[k]]th element is a vector of all lengths of blocks immediately neighboring the block surrounding the position in question.
    flanking.blocks.by.ind = lapply(positions,function(POS){
        tapply(1:length(deme_ID),deme_ID,function(DEME){do.call(c,lapply(DEME,function(IND){
                get.flanking.blocks.all(IND_DATA=sims.sums[[1]]$ind.ancest[[IND]],CHR=chromosome,POS=POS)}))})})
    save(flanking.blocks.by.ind,file=paste(outstring,"flankingBlocks.Robj",sep="_"))
} else {
    load(paste(outstring,"flankingBlocks.Robj",sep="_"))
}

if (!file.exists(paste(outstring,"flankingBlocksAncB.Robj",sep="_"))) {
flanking.blocks.by.ind.AncB = lapply(positions,function(POS){
	tapply(1:length(deme_ID),deme_ID,function(DEME){do.call(c,lapply(DEME,function(IND){
			get.flanking.blocks.AncB(IND_DATA=sims.sums[[1]]$ind.ancest[[IND]],POS=POS)}))})})
    save(flanking.blocks.by.ind.AncB,file=paste(outstring,"flankingBlocksAncB.Robj",sep="_"))
} else {
    load(paste(outstring,"flankingBlocksAncB.Robj",sep="_"))
}

        
# this is a (deme x position) matrix giving mean flanking block length
flanking.blocks.deme.matrix = do.call(rbind,lapply(1:params$ndemes,function(DEME){sapply(flanking.blocks.by.ind,function(POS){mean(POS[[DEME]],na.rm=TRUE)})}))

# this plots mean flanking block length along the genome, with lines per deme,
pdf(height=4, width=6.25,file=paste(outstring,"adjacentBlocksAlongChromNoConditioning_nonnormalized.pdf",sep="_"))
par(mar=c(3.5,3.5,1.5,0.5))
matplot(100*(positions-0.5), 100*t(flanking.blocks.deme.matrix)[,params$ndemes:1], 
        lty=1,col=rainbow(params$ndemes),type="l",
        xlab="Distance from selected locus (cM)",
        ylab="Mean flanking block length",
        main=paste("Mean flanking block length",paramstring) )
.spatial.legend()
dev.off()

# this plots mean flanking block length along the genome, with lines per deme,
# and divided by the mean length for that deme
pdf(height=4, width=6.25,file=paste(outstring,"adjacentBlocksAlongChromNoConditioning.pdf",sep="_"))
par(mar=c(3.5,3.5,0.5,0.5))
matplot(t(flanking.blocks.deme.matrix/apply(flanking.blocks.deme.matrix,1,mean,na.rm=TRUE))[,params$ndemes:1],xaxt="n",yaxt="n",xlab="",ylab="",lty=1,col=rainbow(params$ndemes),type="l",x=(positions-0.5),ylim=c(0.5,3),
        main=paste("Normalized mean flanking block length",paramstring) )
axis(1,at=seq(0.48,0.52,0.01),labels=seq(-0.02,0.02,0.01),cex.axis=0.9)
axis(2,at=seq(0.5,3,0.5),cex.axis=0.9,las=2)
mtext("Distance from selected locus (M)",side=1,line=2.5)
mtext("Normalized mean block length",side=2,line=2.5)
.spatial.legend()
dev.off()


# this is a (deme x position) matrix giving mean flanking block length
# only for indivs that are ancestry B at the position
flanking.blocks.deme.matrix.AncB = do.call(rbind,lapply(1:params$ndemes,function(DEME){sapply(flanking.blocks.by.ind.AncB,function(POS){mean(POS[[DEME]],na.rm=TRUE)})}))

# this plots mean flanking block length along the genome, with lines per deme,
# only for indivs that are ancestry B at the position
pdf(height=4, width=6.25,file=paste(outstring,"adjacentBlocksAlongChromAncBConditioning_nonnormalized.pdf",sep="_"))
par(mar=c(3.5,3.5,1.5,0.5))
matplot(100*(positions-0.5), 100*t(flanking.blocks.deme.matrix.AncB)[,params$ndemes:1], 
        lty=1,col=rainbow(params$ndemes),type="l",
        xlab="Distance from selected locus (cM)",
        ylab="Mean flanking block length",
        main=paste("Mean flanking block length",paramstring) )
.spatial.legend()
dev.off()

# this plots mean flanking block length along the genome, with lines per deme,
# only for indivs that are ancestry B at the position,
# divided by the mean value for that deme
pdf(height=4, width=6.25,file=paste(outstring,"adjacentBlocksAlongChromAncBConditioning.pdf",sep="_"))
par(mar=c(3.5,3.5,0.5,0.5))
plot(flanking.blocks.deme.matrix.AncB[1,]/mean(flanking.blocks.deme.matrix.AncB[1,],na.rm=TRUE),
     x=(positions-0.5),ylim=c(0,3),type="l",lty=1,col="white",main="",cex.main=0.8,xlab="",ylab="",yaxt="n",xaxt="n")
for(i in params$ndemes:1){points(flanking.blocks.deme.matrix.AncB[i,]/mean(flanking.blocks.deme.matrix.AncB[i,],na.rm=TRUE),x=(positions-0.5),type="l",col=rainbow(params$ndemes)[params$ndemes-i])}
axis(1,at=seq(0.48,0.52,0.01),labels=seq(-0.02,0.02,0.01),cex.axis=0.9)
axis(2,at=seq(0.5,3,0.5),cex.axis=0.9,las=2)
mtext("Distance from selected locus (M)",side=1,line=2.5)
mtext("Normalized mean block length",side=2,line=2.5)
.spatial.legend()
dev.off()

# this is the (locus x deme) matrix of ratios: (mean length)/(mean flanking length)
ratio.blocks = mean_per_deme/t(flanking.blocks.deme.matrix[params$ndemes:1,])

# this plots the ratio of (mean length)/(mean flanking length), by deme,
pdf(height=4, width=6.25,file=paste(outstring,"ratioAdjacentBlocksAlongChromNoConditioning_nonnormalized.pdf",sep="_"))
par(mar=c(3.5,3.5,0.5,0.5),mgp=c(2.5,1,0))
matplot(100*(positions-0.5), ratio.blocks, type="l",lty=1,
        col=rainbow(params$ndemes),
        xlab="Distance from selected locus (cM)",
        ylab="Ratio",
        main=paste("Mean length/mean flanking length",paramstring))
.spatial.legend()
dev.off()

# heatmap of the same thing
pdf(height=4, width=6.25,file=paste(outstring,"ratioAdjacentBlocksAlongChromHeatmapNoConditioning.pdf",sep="_"))
our_image(ratio.blocks[,params$ndemes:1],
        main=paste("Mean length/mean flanking length",paramstring))
dev.off()

# this plots the ratio of (mean length)/(mean flanking length), by deme,
# and divided by the average in that deme
pdf(height=4, width=6.25,file=paste(outstring,"ratioAdjacentBlocksAlongChromNoConditioning.pdf",sep="_"))
par(mar=c(3.5,3.5,0.5,0.5))
plot(ratio.blocks[,1]/mean(ratio.blocks[,1],na.rm=TRUE),x=(positions-0.5),ylim=c(0.5,3),type="l",lty=1,col="white",cex.main=0.8,xlab="",ylab="",yaxt="n",xaxt="n",
        main=paste("Mean length/mean flanking length",paramstring))
for(i in 1:params$ndemes){points(ratio.blocks[,i]/mean(ratio.blocks[,i],na.rm=TRUE),x=(positions-0.5),type="l",col=rainbow(params$ndemes)[i])}
axis(1,at=seq(0.48,0.52,0.01),labels=seq(-0.02,0.02,0.01),cex.axis=0.9)
axis(2,at=seq(0.5,3,0.5),cex.axis=0.9,las=2)
mtext("Distance from selected locus (M)",side=1,line=2.5)
mtext("Ratio",side=2,line=2.5)
.spatial.legend()
dev.off()

# this is a (locus x deme) matrix giving the ratio of mean lengths as above, but for individuals of ancestry B
ratio.blocks.AncB = mean_per_deme_AncB/t(flanking.blocks.deme.matrix.AncB[params$ndemes:1,])

# this plots the ratio of (mean length)/(mean flanking length), by deme,
# only for individuals that are ancestry B at that position.
pdf(height=4, width=6.25,file=paste(outstring,"ratioAdjacentBlocksAlongChromAncBConditioning_nonnormalized.pdf",sep="_"))
par(mar=c(3.5,3.5,1.5,0.5),mgp=c(2.5,1,0))
matplot(100*(positions-0.5), ratio.blocks.AncB, type="l",lty=1,
        col=rainbow(params$ndemes),
        xlab="Distance from selected locus (cM)",
        ylab="Ratio",
        main=paste("Mean length/mean flanking length, ancestry B",paramstring))
.spatial.legend()
dev.off()

# same thing, heatmap
pdf(height=4, width=6.25,file=paste(outstring,"ratioAdjacentBlocksAlongChromHeatmapAncBConditioning.pdf",sep="_"))
our_image(ratio.blocks.AncB[,params$ndemes:1],
        main=paste("Mean length/mean flanking length, ancestry B",paramstring))
dev.off()

# this plots the ratio of (mean length)/(mean flanking length), by deme,
# only for individuals that are ancestry B at that position.
# and divided by the mean value for that deme
pdf(height=4, width=6.25,file=paste(outstring,"ratioAdjacentBlocksAlongChromAncBConditioning.pdf",sep="_"))
par(mar=c(3.5,3.5,0.5,0.5))
plot(ratio.blocks.AncB[,1]/mean(ratio.blocks.AncB[,1],na.rm=TRUE),x=(positions-0.5),ylim=c(0.0,3),type="l",lty=1,col="white",cex.main=0.8,xlab="",ylab="",yaxt="n",xaxt="n",
        main=paste("Normalized mean length/mean flanking length, ancestry B",paramstring))
for(i in 1:params$ndemes){points(ratio.blocks.AncB[,i]/mean(ratio.blocks.AncB[,i],na.rm=TRUE),x=(positions-0.5),type="l",col=rainbow(params$ndemes)[i])}
axis(1,at=seq(0.48,0.52,0.01),labels=seq(-0.02,0.02,0.01),cex.axis=0.9)
axis(2,at=seq(0.5,3,0.5),cex.axis=0.9,las=2)
mtext("Distance from selected locus (M)",side=1,line=2.5)
mtext("Normalized mean block length",side=2,line=2.5)
.spatial.legend()
dev.off()



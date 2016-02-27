
source("../../sims/sim-fns.R",chdir=TRUE)
library(ReacTran)
library(Matrix)
source("../../notes/tran.1D.R",chdir=TRUE) # fix for log.A
source("../../notes/cline-fns.R",chdir=TRUE)
source("../../sims/chunks_fns.R",chdir=TRUE)


###
# Plot clines at selected and linked sites

get_simparams <- function (fname) {
    strings <- c( sigma="SIGMA", ninds="Ninds", ndemes="ndemes", s="s", tau="tau" )
    lapply( strings, function (x) { as.numeric(gsub(sprintf(".*[/_]%s([0-9.]+)_.*",x),"\\1",fname)) } )
}


positions = seq(0.2,0.8,0.01)

sel.sim.file <- "../../sims/simulation_SIGMA1_Ninds25000_ndemes50_s0.1_tau100_run2016-01-11_simsums.Robj"
load(sel.sim.file)
#CHUNKS_SEL = get.chunks.at.positions(positions,CHR=1)
sel.params <- get_simparams( sel.sim.file )

# get frequencies for selected simulation
demeID = rep(1:sel.params$ndemes, each=sel.params$ninds/sel.params$ndemes)
B_freq = get.ancestry.freqs(
        IND_DATA=sims.sums[[1]]$ind.ancest, 
        demeID=demeID, 
        VECTOR_OF_POS=positions, 
        CHR=1, 
        ndemes=sel.params$ndemes)


neutral.sim.file <- "../../sims/simulation_SIGMA1_Ninds25000_ndemes50_s0_tau100_run2016-01-15_simsums.Robj"
load(neutral.sim.file)
#CHUNKS_NEU = get.chunks.at.positions(positions,CHR=2)
neu.params <- get_simparams( neutral.sim.file )

# get frequencies for neutral simulation
demeID = rep(1:neu.params$ndemes, each=neu.params$ninds/neu.params$ndemes)
B_neu = get.ancestry.freqs(
        IND_DATA=sims.sums[[1]]$ind.ancest, 
        demeID=demeID, 
        VECTOR_OF_POS=positions, 
        CHR=2, 
        ndemes=neu.params$ndemes)


pdf(file="alleleFrequencies_sim.pdf",height=4,width=6,pointsize=10)
par(mar=c(3.5,3.5,0.5,0.5))
    matplot(B_neu/1000,xlab="",ylab="",type="l",lty=1,col="grey",lwd=0.5,x=seq(-24.5,24.6,1))
    for(i in 1:30){
        matpoints(B_freq[,c(i,62-i)]/1000,type="l",col=(rainbow(43))[32-i],lty=1,lwd=0.75,x=seq(-24.5,24.6,1))
    }
    points(B_freq[,31]/1000,col="red",lwd=2,type="l",x=seq(-24.5,24.6,1))
    mtext("Geographic position",side=1,line=2)
    mtext(expression(p[B]),side=2,line=2.5)
    legend("bottomright",legend = c("r=0","r=0.05","r=0.1","r=0.2","r=0.3","no seln."),col=c("red",rainbow(43)[c(6,11,21,31)],"darkgrey"),lwd=c(2,rep(0.75,4),0.5),cex=0.8)
dev.off()


###
# Do a "comparison to theory" equivalent.
#
# from compare-to-theory.Rmd

sel.simchunks.file <- gsub(".Robj","_chunks.Robjs",sel.sim.file)

sim.chroms <- sapply( chunks.rvals, "[", 1 )
sim.rvals <- ifelse( sim.chroms==1, sapply( chunks.rvals, "[", 2 ) - 0.5, 0.5 )
# spatial locations of the demes
xx <- seq_len(sim.params$ndemes)-(sim.params$ndemes+1)/2
theory.params <- list( sigma=trueSigma(sim.params$sigma), s=sim.params$s, tau=sim.params$tau, 
       density=sim.params$ninds/sim.params$ndemes,
       width=sim.params$ndemes )

tt <- seq(0,theory.params$tau,length.out=21)
xgrid <- setup.grid.1D(x.up=-theory.params$width/2, x.down=theory.params$width/2, N=2*theory.params$width)
fgrid <- extend_grid(xgrid)
fwds.soln <- forwards_pde(s=theory.params$s, times=tt, grid=fgrid, sigma=theory.params$sigma )

gcounts <- genotype_counts(chunks[[1]])
# sim frequencies
sfreqs <- sweep( tcrossprod(gcounts,dip.to.hap), 1, 2*rowSums(gcounts), "/" )
# theory frequencies
tfreqs <- cbind( A=fwds.soln[nrow(fwds.soln),-1], B=1-fwds.soln[nrow(fwds.soln),-1] )

pdf(file="alleleFrequencies_sim_comparison.pdf", width=6, height=4, pointsize=10)
    matplot(xx, sfreqs, type='l', lty=1, main=sprintf("Selected allele frequency, t=%d",sim.params$tau), ylab="frequency", xlab='geographic position')
    matlines(fgrid$x.mid, tfreqs, lty=2, lwd=2)
    legend("topright", legend=outer(colnames(sfreqs),c("sim","theory"),paste,sep=" : "), lty=rep(1:2,each=2), col=rep(1:ncol(sfreqs),2))
dev.off()


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
# positions = seq(0.48,0.52,0.001)
# positions = c(0,0.01,0.02,0.05,0.1,0.2,0.3)
#cols = c(2,3,6,11,21,31)

sel.sim.file <- "../../../sims/simulation_SIGMA1_Ninds25000_ndemes50_s0.1_tau100_run2016-01-11_simsums.Robj"
sel.params <- get_simparams( sel.sim.file )

neutral.sim.file <- "../../../sims/simulation_SIGMA1_Ninds25000_ndemes50_s0_tau100_run2016-01-15_simsums.Robj"
neu.params <- get_simparams( neutral.sim.file )

long.sim.file <- "../../../sims/simulation_SIGMA3_Ninds250000_ndemes500_s0.01_dir/tau100_dir/results_runid_552675_simsums.Robj"
long.params <- get_simparams(long.sim.file)

if (!file.exists("plot_frequency_clines_revision.Robj")) {
    load(sel.sim.file)
    #CHUNKS_SEL = get.chunks.at.positions(positions,CHR=1)

    # get frequencies for selected simulation
    demeID = rep(1:sel.params$ndemes, each=sel.params$ninds/sel.params$ndemes)
    B_freq = get.ancestry.freqs(
            IND_DATA=sims.sums[[1]]$ind.ancest, 
            demeID=demeID, 
            VECTOR_OF_POS=0.5-positions, 
            CHR=1, 
            ndemes=sel.params$ndemes)


    load(neutral.sim.file)
    #CHUNKS_NEU = get.chunks.at.positions(positions,CHR=2)

    # get frequencies for neutral simulation
    demeID = rep(1:neu.params$ndemes, each=neu.params$ninds/neu.params$ndemes)
    B_neu = get.ancestry.freqs(
            IND_DATA=sims.sums[[1]]$ind.ancest, 
            demeID=demeID, 
            VECTOR_OF_POS=seq(0.2,0.8,0.01), 
            CHR=1, 
            ndemes=neu.params$ndemes)

    save(B_freq, B_neu, file="plot_frequency_clines_revision.Robj")
} else {
    load("plot_frequency_clines_revision.Robj")
}

if (!file.exists("plot_frequency_clines_long_revision.Robj")) {

    load(long.sim.file)
    # get frequencies for simulation on wider range
    demeID = rep(1:long.params$ndemes, each=long.params$ninds/long.params$ndemes)
    B_long = get.ancestry.freqs(
            IND_DATA=sims.sums[[1]]$ind.ancest, 
            demeID=demeID, 
            VECTOR_OF_POS=0.5-positions, 
            CHR=1, 
            ndemes=long.params$ndemes)
    B_long_neu = get.ancestry.freqs(
            IND_DATA=sims.sums[[1]]$ind.ancest, 
            demeID=demeID, 
            VECTOR_OF_POS=0.5-positions, 
            CHR=2, 
            ndemes=long.params$ndemes)



    save(B_long, B_long_neu, file="plot_frequency_clines_long_revision.Robj")
} else {
    load("plot_frequency_clines_long_revision.Robj")
}


pdf(file="alleleFrequencies_sim_s0.01_tau100_closely_linked.pdf",height=4,width=6,pointsize=10)
par(mar=c(3.5,3.5,0.5,0.5))
    matplot(1-B_neu/1000,xlab="",ylab="",type="l",lty=1,col="grey",lwd=0.5,x=seq(-24.5,24.6,1),ylim=c(0,1))
    for(i in 1:30){
        matpoints(1-B_freq[,c(i,62-i)]/1000,type="l",col=(rainbow(43))[32-i],lty=1,lwd=0.75,x=seq(-24.5,24.6,1))
    }
    points(1-B_freq[,31]/1000,col="red",lwd=2,type="l",x=seq(-24.5,24.6,1))
    mtext("Geographic position",side=1,line=2)
    mtext(expression(p[B]),side=2,line=2.5)
    legend("bottomright",legend = c("r=0","r=0.005","r=0.01","r=0.02","r=0.03","no seln."),col=c("red",rainbow(43)[c(6,11,21,31)],"darkgrey"),lwd=c(2,rep(0.75,4),0.5),cex=0.8)
dev.off()


###
# Do a "comparison to theory" equivalent.
#
# from compare-to-theory.Rmd

# spatial locations of the demes
xx <- seq_len(sel.params$ndemes)-(sel.params$ndemes+1)/2
theory.params <- list( 
              sigma=trueSigma(sel.params$sigma), 
              s=sel.params$s, 
              tau=sel.params$tau, 
              density=sel.params$ninds/sel.params$ndemes,
              width=sel.params$ndemes )

tt <- seq(0,theory.params$tau,length.out=21)
xgrid <- setup.grid.1D(x.up=-theory.params$width/2, x.down=theory.params$width/2, N=2*theory.params$width)
fgrid <- extend_grid(xgrid)
fwds.soln <- forwards_pde(s=theory.params$s, times=tt, grid=fgrid, sigma=theory.params$sigma )
fwds.neutral.soln <- forwards_pde(s=0, times=tt, grid=fgrid, sigma=theory.params$sigma )

# compare.k <- 1+(0:(floor(length(positions)/5)-1))*5 
# compare.k <- compare.k[ positions[compare.k]>=0.5 ]
compare.k <- 31+c(0,1,2,3,5,8,13,20,30) 
#compare.positions <- positions[ compare.k ]
#compare.rvals <- compare.positions - 0.5
compare.rvals <- positions
linked.solns <- lapply( compare.rvals, function (r) {
        forwards_backwards_pde(s=theory.params$s, times=tt, grid=xgrid, r=r, 
                           sigma=theory.params$sigma, fwds.soln=fwds.soln, log.p=TRUE)
    } )
linked.clines <- lapply( linked.solns, cline_from_soln, grid=xgrid )
tfreqs <- sapply( linked.clines, function (lnc) { lnc[nrow(lnc),-1] } )

linked.neutral.solns <- lapply( compare.rvals, function (r) {
        forwards_backwards_pde(s=0, times=tt, grid=xgrid, r=r, 
                           sigma=theory.params$sigma, fwds.soln=fwds.neutral.soln, log.p=TRUE)
    } )
linked.neutral.clines <- lapply( linked.neutral.solns, cline_from_soln, grid=xgrid )
neutral.tfreqs <- sapply( linked.neutral.clines, function (lnc) { lnc[nrow(lnc),-1] } )

pdf(file="/home/alisa/alleleFrequencies_sim_comparison_revision.pdf", width=6, height=8, pointsize=10)
layout((1:2))
par(mar=c(3.5,3.5,0.5,0.5))

#COLORS=c("red",sapply(1:30,function(i){(rainbow(43))[32-i]}),sapply(30:1,function(i){(rainbow(43))[32-i]}))
COLORS=c("red",rainbow(43)[c(2,3,6,11,21,31)])
    # Neutral comparison
    matplot( xx, B_neu/(2*neu.params$ninds/neu.params$ndemes),
            xlab="",ylab="",type="l",lty=1,col="grey",lwd=0.5 )
    matlines(xgrid$x.mid, 1-neutral.tfreqs, lty=2, lwd=2, col='black' )

        mtext("Geographic position",side=1,line=2)
        mtext(expression(p[B]),side=2,line=2.5)

    # selected comparision
    r.cols <- rainbow(ncol(B_freq)+10)
    matplot( xx, B_freq/(2*sel.params$ninds/sel.params$ndemes),
            xlab="",ylab="",type="l",lty=1,lwd=0.75,
            col=COLORS )
    matlines(xgrid$x.mid, 1-tfreqs, lty=2, lwd=2, 
            col=COLORS)

        mtext("Geographic position",side=1,line=2)
        mtext(expression(p[B]),side=2,line=2.5)

    legend("bottomright",legend = c("r=0.01","r=0.05","r=0.1","r=0.2","r=0.3","no seln."),col=c(COLORS[c(2,4,5,6,7)],"darkgrey"),lwd=c(2,rep(0.75,4),0.5),cex=0.8)

dev.off()

###
## longer region

xx = (1:long.params$ndemes)-0.5-long.params$ndemes/2

pdf(file="alleleFrequencies_sim_SIGMA3_Ninds250000_ndemes500_s0.01_tau100.pdf",height=4,width=6,pointsize=10)
par(mar=c(3.5,3.5,0.5,0.5))
    matplot(xx, B_long_neu[,1:30]/1000, xlab="", ylab="", type="l", lty=1, col="grey", lwd=0.5, ylim=c(0,1))
    for(i in 1:30){
        matlines(xx, B_long[,i]/1000, col=(rainbow(43))[32-i], lty=1, lwd=0.75 )
    }
    lines(xx, B_long[,31]/1000, col="red", lwd=2 )
    mtext("Geographic position",side=1,line=2)
    mtext(expression(p[B]),side=2,line=2.5)
    legend("bottomright",legend = c("r=0","r=0.005","r=0.01","r=0.02","r=0.03","no seln."),col=c("red",rainbow(43)[c(6,11,21,31)],"darkgrey"),lwd=c(2,rep(0.75,4),0.5),cex=0.8)
dev.off()


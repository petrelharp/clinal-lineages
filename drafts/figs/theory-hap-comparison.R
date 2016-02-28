
### compare theory and simulation haplotype stats
### from compare-to-theory.Rmd

source("../../sims/sim-fns.R",chdir=TRUE)
library(ReacTran)
library(Matrix)
source("../../notes/tran.1D.R",chdir=TRUE) # fix for log.A
source("../../notes/cline-fns.R",chdir=TRUE)
source("../../sims/chunks_fns.R",chdir=TRUE)

for (simchunks.file in c( 
                         "../../sims/simulation_SIGMA3_Ninds5000_ndemes100_s0.05_dir/tau80_dir/results_runid_281321_simsums_chunks.Robj",
                         "../../sims/simulation_SIGMA3_Ninds5000_ndemes100_s0.05_dir/tau320_dir/results_runid_102364_simsums_chunks.Robj",
                         "../../sims/simulation_SIGMA3_Ninds5000_ndemes100_s0.05_dir/tau1280_dir/results_runid_379763_simsums_chunks.Robj" )) {

    load(simchunks.file)  # provides 'chunks' and 'chunks.rvals'

    get_simparams <- function (fname) {
        strings <- c( sigma="SIGMA", ninds="Ninds", ndemes="ndemes", s="s", tau="tau" )
        lapply( strings, function (x) { as.numeric(gsub(sprintf(".*[/_]%s([0-9.]+)_.*",x),"\\1",fname)) } )
    }
    sim.params <- get_simparams( simchunks.file )
    # recombination distances
    # sim.rvals <- c(0,0.5-0.01,0.5)
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

    # selected clines

    sel.gcounts <- genotype_counts(chunks[[1]])
    # sim frequencies
    sel.sfreqs <- sweep( tcrossprod(sel.gcounts,dip.to.hap), 1, 2*rowSums(sel.gcounts), "/" )
    # theory frequencies
    sel.tfreqs <- cbind( A=fwds.soln[nrow(fwds.soln),-1], B=1-fwds.soln[nrow(fwds.soln),-1] )

    # linked clines

    linked.solns <- lapply( sim.rvals, function (r) {
            forwards_backwards_pde(s=theory.params$s, times=tt, grid=xgrid, r=r, 
                               sigma=theory.params$sigma, fwds.soln=fwds.soln, log.p=TRUE)
        } )

    # could use these to matplot all at once
    # # simulation
    # all.gcounts <- lapply(chunks[-1], conditional_genotype_counts, chunks[[1]])
    # all.cfreqs <- sapply(all.gcounts,conditional_freqs)
    # dim(all.cfreqs) <- c( dim(all.cfreqs)[1]/2, 2, dim(all.cfreqs)[2] )

    # # theory
    # all.tfreqs <- sapply( linked.solns, function (x) {x[nrow(x),-1]} )
    # dim(all.tfreqs) <- c( length(xgrid$x.mid), 2, dim(all.tfreqs)[2] )

    # plot it!

    pdf(file=sprintf("cond_freqs_comparison_tau_%d.pdf",sim.params$tau), width=6.5, height=9, pointsize=10)
    layout(matrix(seq_along(sim.rvals),ncol=2,byrow=TRUE))

    matplot(xx, sweep(sel.gcounts,1,rowSums(sel.gcounts),"/"), lty=1, type='l', 
            main=sprintf("Selected genotype frequencies, t=%d",sim.params$tau), 
            ylab="frequency", xlab='geographic position')
    matlines(fgrid$x.mid, cbind(AA=sel.tfreqs[,"A"]^2,AB=2*sel.tfreqs[,"A"]*sel.tfreqs[,"B"],BB=sel.tfreqs[,"B"]^2), lty=4, lwd=2)
    legend("topright", lty=rep(1:2,each=3), col=rep(1:3,2), 
           legend=outer(c("AA","AB","BB"),c("simulation","theory"),paste,sep=": ") )

    for (kr in seq_along(sim.rvals)[-1]) {
        gcounts <- conditional_genotype_counts(chunks[[kr]], chunks[[1]])
        # frequencies from simulation
        cfreqs <- conditional_freqs(gcounts)
        # frequencies from theory
        tfreqs <- linked.solns[[kr]][nrow(linked.solns[[kr]]),-1]
        dim(tfreqs) <- c( length(xgrid$x.mid), 2 )
        matplot(xx, cfreqs, type='l', lty=1, 
                main=sprintf("A ancestry by selected background, t=%d, r=%0.3f",sim.params$tau, sim.rvals[kr]), 
                ylab="frequency", xlab='geographic position')
        matlines(xgrid$x.mid, tfreqs, lty=2, lwd=2)
        legend("topright", legend=outer(colnames(cfreqs),c("sim","theory"),paste,sep=" : "), 
               lty=rep(1:2,each=2), col=rep(1:ncol(cfreqs),2))
    }

    dev.off()
}

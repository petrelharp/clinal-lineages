# Get fine-scale haplotype length predictions.

source("../../sims/sim-fns.R",chdir=TRUE)
library(ReacTran)
library(Matrix)
source("../../notes/tran.1D.R",chdir=TRUE) # fix for log.A
source("../../notes/cline-fns.R",chdir=TRUE)
source("../../sims/chunks_fns.R",chdir=TRUE)

theory.params <- list( 
              sigma=1L,
              s=0.01,
              tau=1000L,
              density=100L,
              ndemes=50L )
theory.params$ninds <- theory.params$ndemes * theory.params$density

filebase <- sprintf("haplotypes_SIGMA%d_density%d_ndemes%d_s%0.3f_tau%d",
                 theory.params$sigma, theory.params$density, theory.params$ndemes, theory.params$s, theory.params$tau )
robj.file <- file.path("cache",paste0(filebase,".Robj"))

if (!file.exists) {
    tt <- seq(0,theory.params$tau,length.out=21)
    xgrid <- setup.grid.1D(x.up=-theory.params$ndemes/2, x.down=theory.params$ndemes/2, N=50)
    fgrid <- extend_grid(xgrid)
    fwds.soln <- forwards_pde(s=theory.params$s, times=tt, grid=fgrid, sigma=theory.params$sigma )

    # do it on a fine scale: for tau=1000, looks to be over +/- 5 or 6 cM
    rgrid <- setup.grid.1D(x.up=-0.05, x.down=0.05, N=40)
    hap.soln <- forwards_backwards_haplotypes(s=theory.params$s, times=tt, xgrid=xgrid, rgrid=rgrid,
                               sigma=theory.params$sigma, 
                               fwds.grid=fgrid, fwds.soln=fwds.soln )

    dir.create("cache",showWarnings=FALSE)
    save(theory.params, filebase, tt, xgrid, fgrid, fwds.soln, rgrid, hap.soln, file=robj.file)
} else {
    load(robj.file)
}

####

hap.pdf <- hap_cdf_to_pdf( hap.soln )


hap.len.mat <- sapply( rgrid$x.mid, function (rr) { 
               haplens <- mean_hap_len( hap.pdf, loc=rr )
               haplens[length(tt),1+seq_along(xgrid$x.mid)]
           } )

# mean length against space
matplot(rev(xgrid$x.mid), hap.len.mat, type='l', lty=1,
        main="mean B haplotype length", 
        col=rainbow(length(rgrid$x.mid)),
        xlab="spatial position", ylab="mean B haplotype length (incl. zeros)")
k.legend <- floor(seq(1,ncol(hap.len.mat),length.out=8))
legend("topleft", lty=1, col=rainbow(length(rgrid$x.mid))[k.legend],
        legend=sprintf("r=%0.0f cM",100*rgrid$x.mid[k.legend]) )

# mean length along the genome
matplot(rgrid$x.mid*100, t(hap.len.mat), type='l', lty=1,
        main="mean B haplotype length", 
        col=rainbow(length(xgrid$x.mid)),
        xlab="position along genome (cM)", 
        ylab="mean B haplotype length (incl. zeros)")
k.legend <- floor(seq(1,nrow(hap.len.mat),length.out=8))
legend("topleft", lty=1, col=rainbow(length(xgrid$x.mid))[k.legend],
        legend=sprintf("x=%0.2f",xgrid$x.mid[k.legend]) )

# mean length along the genome, normalized
matplot(rgrid$x.mid*100, t(sweep(hap.len.mat,1,rowMeans(hap.len.mat),"/")), 
        type='l', lty=1,
        main="mean B haplotype length, normalized by location", 
        col=rainbow(length(xgrid$x.mid)),
        xlab="position along genome (cM)", 
        ylab="relative mean B haplotype length (incl. zeros)")
k.legend <- floor(seq(1,nrow(hap.len.mat),length.out=8))
legend("topleft", lty=1, col=rainbow(length(xgrid$x.mid))[k.legend],
        legend=sprintf("x=%0.2f",xgrid$x.mid[k.legend]) )

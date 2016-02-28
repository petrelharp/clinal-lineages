# Get fine-scale haplotype length predictions.

source("../../sims/sim-fns.R",chdir=TRUE)
library(ReacTran)
library(Matrix)
source("../../notes/tran.1D.R",chdir=TRUE) # fix for log.A
source("../../notes/cline-fns.R",chdir=TRUE)
source("../../sims/chunks_fns.R",chdir=TRUE)

theory.params <- list( 
              sigma=1,
              s=0.01,
              tau=1000,
              density=100,
              width=50 )

tt <- seq(0,theory.params$tau,length.out=21)
xgrid <- setup.grid.1D(x.up=-theory.params$width/2, x.down=theory.params$width/2, N=50)
fgrid <- extend_grid(xgrid)
fwds.soln <- forwards_pde(s=theory.params$s, times=tt, grid=fgrid, sigma=theory.params$sigma )

# do it on a fine scale: for tau=1000, looks to be over +/- 5 or 6 cM
rgrid <- setup.grid.1D(x.up=-0.05, x.down=0.05, N=40)
hap.soln <- forwards_backwards_haplotypes(s=theory.params$s, times=tt, xgrid=xgrid, rgrid=rgrid,
                           sigma=theory.params$sigma, 
                           fwds.grid=fgrid, fwds.soln=fwds.soln )

dir.create("cache",showWarnings=FALSE)
save(list=ls(), file=sprintf("cache/haplotypes_SIGMA%d_density%d_ndemes%d_s%0.3f_tau%d.Robj", 
                             theory.params$sigma, theory.params$density, theory.params$width, theory.params$s, theory.params$tau ) )

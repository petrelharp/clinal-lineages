library(ReacTran)
library(Matrix)
source("../tran.1D.R",chdir=TRUE) # fix for log.A
source("../cline-fns.R",chdir=TRUE)
s <- 0.1
tt <- seq(0,400,length.out=21)
xgrid <- setup.grid.1D(x.up=-20.5, x.down=20.5, N=41)
rr <- c(0.05,0.2)
rcols <- rainbow(length(rr)*1.2)[seq_along(rr)]
fwds.soln <- forwards_pde(s=s,times=tt,grid=extend_grid(xgrid))
pfun <- cline_fun(fwds.soln)
linked.solns <- lapply( rr, function (r) {
            forwards_backwards_pde(s=s,times=tt,grid=xgrid,r=r,log.p=TRUE,fwds.soln=fwds.soln)
        } )
linked.clines <- lapply( linked.solns, cline_from_soln, grid=xgrid )

# check forwards solution is symmetric
xx <- fwds.soln[,-1]
yy <- 1-xx[,ncol(xx):1]
stopifnot( all( abs(as.vector(xx)-as.vector(yy)) < 1e-12 ) )

# and reverse?
lclines <- sapply( linked.clines, cline_interp, t=200 )
stopifnot( all( abs(lclines-(1-lclines[nrow(lclines):1,])) < 1e-12 ) )


cat("Everything checks out!\n")

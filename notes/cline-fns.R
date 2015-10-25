pcline <- function (x,r,p=recombfn(x,2*r,...),recombfn,...) {
    (1/2)*( 1 - sign(x)*p )
}

tanh_cline <- function (x,s,sigma=1,...) {
    (1/2)*(1+tanh(-2*x*sqrt(s)/sigma))
}

tanh_drift <- function (x,s,sigma=1) {
    # log-derivative of tanh_cline
    -(2*sqrt(s)/sigma) * (1-tanh(-2*x*sqrt(s)/sigma)^2) / (1+tanh(-2*x*sqrt(s)/sigma))
}

tanh_het <- function (x,s,sigma=1) {
    # prob of being in a het = (1-s)(1-p)/(1-s(1-p))
    pp <- tanh_cline(x,s,sigma=sigma)
    (1-s)*(1-pp)/(1-s*(1-pp))
}

zeta.reflected <- function (x,r) {
    2*pnorm(abs(x)) - 1 + 2 * exp( r*abs(x) + r^2/2 + pnorm(abs(x)+r,lower.tail=FALSE,log.p=TRUE) )
}

sim_diffusion <- function (drift,tau,nwalks,ntimes,xinit=0,...) {
    # want total time to be tau, so:
    dt <- tau/ntimes
    # put the noise in rws
    rws <- matrix( c( rep_len(xinit,nwalks), # initial conditions
                 rnorm(nwalks*(ntimes-1)) ), nrow=ntimes, ncol=nwalks, byrow=TRUE )
    # and then integrate over it
    for (k in 2:ntimes) {
        rws[k,] <- rws[k-1,]+sqrt(dt)*rws[k,]+dt*drift(rws[k-1,],...)
    }
    attr(rws,"dt") <- dt
    return(rws)
}

weighted_time <- function (rws, weight.fn, use.t=TRUE, ...) {
    colSums( weight.fn(rws[use.t,],...) ) * attr(rws,"dt")
}

clineplot <- function (soln,x=thegrid$x.mid,times=seq(1,nrow(soln),length.out=ntimes),ntimes=20,ylab="probability",ylim=c(0,1),legend=TRUE,...) {
    # plot output of deSolve as clines in a nice way
    graphics::matplot( x, t(soln[times,-1]), type='l', col=rainbow(length(times)), lty=1, 
            xlab="space", ylab=ylab, ylim=ylim, ...)
    if (legend) legend("bottomleft", lty=1, col=rainbow(length(times)), legend=paste("t=",round(soln[times,1],digits=2)), cex=0.5)
    return( invisible( t(soln[times,-1]) ) )
}

cline_interp <- function (T,soln) {
    # return the interpolated cline for a particular time
    # from the output of deSolve
    which.t <- findInterval(T,soln[,1],rightmost.closed=TRUE)
    lu <- soln[c(which.t,which.t+1),1]
    alpha <- (T-lu[1])/diff(lu)
    return( (1-alpha)*soln[which.t,-1] + alpha*soln[which.t+1,-1]  )
}

solve.pde <- function ( u, v=u, r, times, grid, ... ) {
    # Solve the pde for a lineage 
    #  where 'u' gives the drift and 'v' the rate of recombination.
    vfn <- function(x,t){v(x=x,t=t,...)}
    u.vec <- u(x=grid$x.mid,t=t,...)
    v.vec <- v(x=grid$x.mid,t=t,...)
    yinit <- c( rep(1.0,grid$N), rep(0.0,grid$N) )
    r <- 0.1
    pde.fn <- function (t,y,parms,...) {
        yA <- y[1:grid$N]
        yB <- y[grid$N+(1:grid$N)]
        tran.A <- tran.1D(C=yA, A=u.vec, D=1/2, dx=grid)$dC
        tran.B <- tran.1D(C=yB, A=1-u.vec, D=1/2, dx=grid)$dC
        # if (any(is.na(tran.A)|is.na(tran.B))) { browser() }
        list( c( 
                tran.A + r * (1-v.vec) * (yB-yA),
                tran.B + r * v.vec * (yA-yB)
                ) )
    }
    soln <- ode.1D( y=yinit, times=times, func=pde.fn, nspec=2 ) # note FIRST COLUMN IS TIME
    attr(soln,"u") <- function(x,t) {u(x=x,...)}
    attr(soln,"v") <- vfn
    return(soln)
}
cline.from.soln <- function (soln, grid, clinefn=attr(soln,"u") ) {
    # convert the output from solve.pde to the cline
    # by averaging over which allele was sampled
    A.inds <- 1+(1:grid$N)
    B.inds <- 1+grid$N+(1:grid$N)
    sel.cline <- outer( grid$x.mid, soln[nrow(soln),1], clinefn )
    pde.cline <- cbind( pde.soln[,1], pde.soln[,A.inds]*sel.cline[col(pde.soln[,A.inds])] + pde.soln[,B.inds]*(1-sel.cline)[col(pde.soln[,B.inds])] )
    attributes(pde.cline) <- c( attributes(pde.cline), attributes(pde.soln)[setdiff(names(attributes(pde.soln)),c("dim","dimnames"))] )
    attr(pde.cline,"nspec")<-1
    return(pde.cline)
}


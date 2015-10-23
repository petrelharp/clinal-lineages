pcline <- function (x,r,p=recombfn(x,2*r,...),recombfn,...) {
    (1/2)*( 1 - sign(x)*p )
}
tanh_cline <- function (x,s,sigma=1) {
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
cline.interp <- function (T,soln) {
    # return the interpolated cline for a particular time
    # from the output of deSolve
    which.t <- findInterval(T,soln[,1],rightmost.closed=TRUE)
    lu <- soln[c(which.t,which.t+1),1]
    alpha <- (T-lu[1])/diff(lu)
    return( (1-alpha)*soln[which.t,-1] + alpha*soln[which.t+1,-1]  )
}


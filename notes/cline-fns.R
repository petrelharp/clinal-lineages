source("unifprod.R")  # for the C part

pcline <- function (x,r,p=recombfn(x,2*r,...),recombfn,...) {
    (1/2)*( 1 - sign(x)*p )
}

tanh_cline <- function (x,s,sigma=1,log=FALSE,...) {
    if (log) {
        - (4*x*sqrt(s)/sigma) - log1p(exp(-4*x*sqrt(s)/sigma))
    } else {
        (1/2)*(1+tanh(-2*x*sqrt(s)/sigma))
    }
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

clineplot <- function (soln, thegrid, x=thegrid$x.mid,
               times=seq(1,nrow(soln),length.out=ntimes), ntimes=20,
               ylab="probability", ylim=c(0,1), legend=TRUE, lty=1,
               col=rainbow(length(times)),
           ...) {
    # plot output of deSolve as clines in a nice way
    graphics::matplot( x, t(soln[times,-1]), type='l', col=col, lty=lty, 
            xlab="space", ylab=ylab, ylim=ylim, ...)
    if (legend) legend("bottomleft", lty=lty, col=col, legend=paste("t=",round(soln[times,1],digits=2)), cex=0.5)
    return( invisible( t(soln[times,-1]) ) )
}

add_grid_axis <- function (thegrid,side=2,x=thegrid$x.mid,...) {
    axis( side=side, at=(pretty(x)-min(x))/diff(range(x)), labels=pretty(x))
}

cline_interp <- function (t,soln) {
    # return the interpolated cline for a particular time
    # from the output of deSolve
    sapply(t, function (TT) {
            which.t <- findInterval(TT,soln[,1],rightmost.closed=TRUE)
            lu <- soln[c(which.t,which.t+1),1]
            alpha <- (TT-lu[1])/diff(lu)
            return( (1-alpha)*soln[which.t,-1] + alpha*soln[which.t+1,-1]  )
        } )
}

cline_fun <- function (soln) {
    # return a function that interpolates the solution to (t,x)
    function (t,x) { 
            fp <- function (x,t) { approx(attr(soln,"grid")$x.mid,cline_interp(t,soln=soln),xout=x)$y }
            ans <- numeric(length(x))
            for (tval in unique(t)) {
                dothese <- (rep_len(t,length(x))==tval)
                ans[dothese] <- fp(x=x[dothese],t=tval)
            }
            return(ans)
    }
}

forwards_pde <- function (s,times,grid,sigma=1,
                      yinit=ifelse( grid$x.mid==0, 0.5, ifelse( grid$x.mid>0, 0, 1 )) 
              ) {
    # Solve the pde for frequencies, forwards in time:
    #   by default, yinit=1 to the *left* of zero;
    #   so finds the probability the selected site initiates on the left
    fwds.pde <- function (t,y,parms,...) {
        tran <- tran.1D(C=y, D=sigma^2/2, dx=grid)$dC
        list( tran + s * y * (1-y) * (2*y-1) )
    }
    soln <- ode.1D( y=yinit, times=times, func=fwds.pde, nspec=1 ) # note FIRST COLUMN IS TIME
    attr(soln,"s") <- s
    attr(soln,"grid") <- grid
    return(soln)
}

extend_grid <- function (grid) {
    # extend a given grid to have one more point on the outside
    dx.up <- grid$dx[1]
    dx.down <- grid$dx[grid$N]
    newgrid <- list(
                x.up = grid$x.up-dx.up,
                x.down = grid$x.down+dx.down,
                x.mid = c( grid$x.mid[1]-dx.up, grid$x.mid, grid$x.mid[grid$N]+dx.down ),
                x.int = c( grid$x.int[1]-dx.up, grid$x.int, grid$x.int[grid$N+1]+dx.down ),
                dx = rep_len(grid$dx,length.out=grid$N)[c(1,1:grid$N,grid$N)],
                dx.aux = c( grid$dx.aux[1], dx.up, grid$dx.aux[2:grid$N], dx.down, grid$dx.aux[grid$N+1] ),
                N =grid$N+2
            )
    class(newgrid) <- class(grid)
    return(newgrid)
}

forwards_backwards_pde <- function (s, times, grid, r, sigma=1,
        fwds.grid=extend_grid(grid),  
        fwds.soln=forwards_pde(s=s,times=times,grid=fwds.grid, # needs to be on a *wider* grid
            yinit=ifelse(fwds.grid$x.mid==0,(yinit[1:fwds.grid$N]+yinit[fwds.grid$N+(1:fwds.grid$N)])/2,
                         ifelse(fwds.grid$x.mid<0,yinit[1:fwds.grid$N],yinit[fwds.grid$N+(1:fwds.grid$N)])) 
        ),
        yinit=c( rep(1.0,grid$N), rep(0.0,grid$N) ), 
        log.p=TRUE, eps=1e-16, ... ) {
    # solve the pde for a lineage run backwards in the forwards-time profile
    #   ... note that we solve the PDE in the forwards-time direction in *both* cases.
    # HOWEVER: note that this becomes singular near t=0,
    #   so we push p away from the boundary by 'eps'.
    # By default
    #  - when the cline forms, A is on the left and B is on the right
    #  - fwds.soln is the frequency of the type A selected locus,
    #          i.e., the prob the selected lineage begins on the *left*
    #       * this must be provided on a somewhat *wider* grid.
    #  - yA, the first N numbers, is the probability that a locus linked to a type A selected locus
    #          originally was linked to a type A
    #          i.e., linked to a 'left' allele, and originates from the left side
    #  - yB, the second N numbers, is the probability that a locus linked to a type B selected locus
    #          originally was linked to a type A
    #          i.e., linked to a 'right' allele, but originates from the left side

    if ( (min(attr(fwds.soln,"grid")$x.mid)>grid$x.up) || (max(attr(fwds.soln,"grid")$x.mid)<grid$x.down) ) {
        stop("fwds.soln must be on a wider grid than 'grid' (use extend_grid( ))")
    }

    ufun <- function (x,t) { approx(attr(fwds.soln,"grid")$x.mid,cline_interp(t,soln=fwds.soln),xout=x)$y }
    # ufun is not vectorized in t, so...
    u <- function (x,t) {
        ans <- numeric(length(x))
        for (tval in unique(t)) {
            dothese <- (rep_len(t,length(x))==tval)
            ans[dothese] <- ufun(x[dothese],t=tval)
        }
        return(ans)
    }
    rev.pde.fn <- function (t,y,parms,...) {
        yA <- y[1:grid$N]
        yB <- y[grid$N+(1:grid$N)]
        p <- list( 
                  mid=u(x=grid$x.mid,t=t),
                  int=u(x=grid$x.int,t=t)
              )
        if (log.p) {
            log.p.vec <- log(pmax(p$int,min(eps,min(p$int[p$int>0]))))
            log.1mp.vec <- log(pmax(1-p$int,min(eps,min((1-p$int)[p$int<1]))))
            tran.A <- tran.1D(C=yA, A=log.p.vec, D=sigma^2/2, dx=grid, flux.up=0, flux.down=0, log.A=TRUE)$dC
            tran.B <- tran.1D(C=yB, A=log.1mp.vec, D=sigma^2/2, dx=grid, flux.up=0, flux.down=0, log.A=TRUE)$dC
        } else {
            tran.A <- tran.1D(C=yA, A=pmax(eps,p), D=sigma^2/2, dx=grid, flux.up=0, flux.down=0)$dC
            tran.B <- tran.1D(C=yB, A=pmax(eps,1-p), D=sigma^2/2, dx=grid, flux.up=0, flux.down=0)$dC
        }
        # if (any(is.na(tran.A)|is.na(tran.B))) { browser() }
        list( c( 
                tran.A + r * (1-p$mid) * (yB-yA),
                tran.B + r * p$mid * (yA-yB)
                ) )
    }
    rev.soln <- ode.1D( y=yinit, times=times, func=rev.pde.fn, nspec=2, tcrit=max(fwds.soln[,1]), ... ) # note FIRST COLUMN IS TIME
    attr(rev.soln,"r") <- r
    attr(rev.soln,"s") <- s
    attr(rev.soln,"sigma") <- sigma
    attr(rev.soln,"u") <- u
    attr(rev.soln,"grid") <- grid
    return(rev.soln)
}

solve_pde <- function ( u, v=u, r, times, grid, log.u=FALSE, um1, sigma=1,
                yinit=c( rep(0.0,grid$N), rep(1.0,grid$N) ), ... ) {
    # Solve the pde for a lineage 
    #  moving between backgrounds
    #  where 'u' gives the drift and 'v' the rate of recombination.
    #  and if log.u is TRUE then u is actually log(proportion)
    # Note that um1 should be provided if u is in log scale.
    more.args <- c(list(sigma=sigma),list(...))
    vfn <- function(x,t){do.call(v,c(list(x=x,t=t),more.args))}
    u.vec <- u(x=grid$x.int,t=t,...)
    if (missing(um1)) {
        if (log.u) {
            warning("Computing 1-u from u even though u is in log scale.")
            um1 <- function(...){-expm1(u(...))}
        } else {
            um1 <- function(...) u(...)
        }
    }
    um1.vec <- um1(x=grid$x.int,t=t,...) # 1-u.vec, but maybe log(1-u)
    v.vec <- v(x=grid$x.mid,t=t,...)
    if (log.u) { v.vec <- exp(v.vec) }
    pde.fn <- function (t,y,parms,...) {
        yA <- y[1:grid$N]
        yB <- y[grid$N+(1:grid$N)]
        tran.A <- tran.1D(C=yA, A=u.vec, D=sigma^2/2, dx=grid, flux.up=0, flux.down=0, log.A=log.u)$dC
        tran.B <- tran.1D(C=yB, A=um1.vec, D=sigma^2/2, dx=grid, flux.up=0, flux.down=0, log.A=log.u)$dC
        # if (any(is.na(tran.A)|is.na(tran.B))) { browser() }
        list( c( 
                tran.A + r * (1-v.vec) * (yB-yA),
                tran.B + r * v.vec * (yA-yB)
                ) )
    }
    soln <- ode.1D( y=yinit, times=times, func=pde.fn, nspec=2 ) # note FIRST COLUMN IS TIME
    attr(soln,"u") <- if (log.u) { 
            function(x,t) {exp(do.call(u,c(list(x=x),more.args)))}
        } else {
            function(x,t) {do.call(u,c(list(x=x),more.args))}
        }
    attr(soln,"v") <- vfn
    attr(soln,"r") <- r
    attr(soln,"sigma") <- sigma
    attr(soln,"args") <- list(...)
    attr(soln,"grid") <- grid
    return(soln)
}

cline_from_soln <- function (soln, grid, clinefn=attr(soln,"u") ) {
    # convert the output from solve_pde to the cline
    # by averaging over which allele was sampled
    # Here:
    #  - clinefn returns the local frequency of As
    #  - the first column is the time,
    #  - the next N numbers gives the value linked to As
    #  - and the next N numbers gives the value linked to Bs
    A.inds <- 1+(1:grid$N)
    B.inds <- 1+grid$N+(1:grid$N)
    sel.cline <- outer( grid$x.mid, soln[,1], clinefn )
    pde.cline <- cbind( soln[,1], soln[,A.inds]*sel.cline[col(soln[,A.inds])] + soln[,B.inds]*(1-sel.cline)[col(soln[,B.inds])] )
    attributes(pde.cline) <- c( attributes(pde.cline), attributes(soln)[setdiff(names(attributes(soln)),c("dim","dimnames"))] )
    attr(pde.cline,"nspec")<-1
    return(pde.cline)
}


####
## Linkage-related things
####

recomb_generator <- function (yA,yB,p,zeroind,nabvals,eps) {
    # the discretization of the map y(a,b) -> \int_a^b (h(a,u)*h(u,b)-h(a,b)) du
    # yA and yA is the upper triangular part of square matrices
    #   of dimensions nabvals x nabvals
    #   eps is the grid spacing between locations
    #   and zeroind is the index corresponding to the marker at zero
    eps * ( 
      p * unifprod_ut(y=yA,ncx=nabvals)
      + (1-p) * biprod_ut(yA=yA,yB=yB,zeroind=zeroind,ncx=nabvals)
    )
}

forwards_backwards_haplotypes <- function (s, times, xgrid, rgrid, sigma=1,
        fwds.grid=extend_grid(xgrid),  
        fwds.soln=forwards_pde(s=s,times=times,grid=fwds.grid,
            yinit=ifelse(fwds.grid$x.mid==0,(yinit[1:fwds.grid$N]+yinit[fwds.grid$N+(1:fwds.grid$N)])/2,
                         ifelse(fwds.grid$x.mid<0,yinit[1:fwds.grid$N],yinit[fwds.grid$N+(1:fwds.grid$N)])) 
        ),
        yinit=rep( c( rep(1.0,xgrid$N), rep(0.0,xgrid$N) ), rgrid$N*(rgrid$N+1)/2 ),
        eps=1e-16, ... ) {
    ###
    # solves the system of PDE for the probability that a haplotype across genomic region [a,b]
    # sampled at time t and location x and linked to selected allele z
    # derives entirely from side A at time 0.
    ###
    # s = selection coefficient
    # times = times to solve at
    # xgrid = grid of spatial locations
    # rgrid = grid of recombination values
    ###
    # note: doing parallel lapply below is actually a good bit slower.
    ###
    # Result is ordered by:
    #  - spatial location
    #  - linked background
    #  - left/right endpoints (in upper triangular order)
    ###
    # so can be indexed e.g. by (where hap.soln is the output)
    # hap.coords <- expand.grid( x=xgrid$x.mid, AB=c("A","B"), rk=1:nrow(attr(hap.soln,"r")) )
    # hap.coords$p <- hap.soln[nrow(hap.soln),-1]
    # hap.coords$right <- attr(hap.soln,"r")[hap.coords$rk,"right"]
    # hap.coords$left <- attr(hap.soln,"r")[hap.coords$rk,"left"]
    # hap.coords$len <- hap.coords$right-hap.coords$left
    ####
    # precompute some things
    rgrid$eps <- unique(diff(rgrid$x.mid))[1]
    rgrid$zeroind <- which(rgrid$x.mid>0)[1]

    # function to find a,b coordinates from the upper triangular coordinates
    ab.ord <- matrix( 0, nrow=rgrid$N, ncol=rgrid$N )
    avals <- rgrid$x.mid[row(ab.ord)[upper.tri(ab.ord,diag=TRUE)]]
    bvals <- rgrid$x.mid[col(ab.ord)[upper.tri(ab.ord,diag=TRUE)]]
    ab.coord <- function (ab) { cbind( a=avals[ab], b=bvals[ab] ) }
    # distance from the closest end of the segment [a,b] to zero
    rvals <- ifelse( avals*bvals<0, 0, pmin(abs(avals),abs(bvals)) )

    # we'll save the solution as an array,
    #  with dimensions (x,ab)
    # Here ab is upper-triangular ordered
    soln.dims <- c( nrow=2*xgrid$N, ncol=rgrid$N*(rgrid$N+1)/2 )

    # the local frequency of the selected allele
    ufun <- function (x,t) { approx(attr(fwds.soln,"grid")$x.mid, cline_interp(t,soln=fwds.soln), xout=x)$y }
    # not vectorized in t, so...
    u <- function (x,t) {
        ans <- numeric(length(x))
        for (tval in unique(t)) {
            dothese <- (rep_len(t,length(x))==tval)
            ans[dothese] <- ufun(x[dothese],t=tval)
        }
        return(ans)
    }

    # the gradient function
    rev.pde.fn <- function (t,y,parms,...) {
        # need to apply tran.1D to columns
        # and recomb_generator to rows
        dim(y) <- soln.dims
        p <- list( 
                  mid=u(x=xgrid$x.mid,t=t),
                  int=u(x=xgrid$x.int,t=t)
              )
        log.p <- log( pmax(p$int, min(eps,min(p$int[p$int>0]))) )
        log.1mp <- log( pmax(1-p$int, min(eps,min((1-p$int)[p$int<1]))) )
        diffusion <- do.call( cbind, lapply( 1:ncol(y), function (kcol) {
                r <- rvals[kcol]
                yA <- y[1:xgrid$N,kcol]
                yB <- y[xgrid$N+(1:xgrid$N),kcol]
                tran.A <- tran.1D(C=yA, A=log.p, D=sigma^2/2, dx=xgrid, flux.up=0, flux.down=0, log.A=TRUE)$dC
                tran.B <- tran.1D(C=yB, A=log.1mp, D=sigma^2/2, dx=xgrid, flux.up=0, flux.down=0, log.A=TRUE)$dC
                c( 
                        tran.A + r * (1-p$mid) * (yB-yA),
                        tran.B + r * p$mid * (yA-yB)
                        )
             } ) )
        recombination <- do.call( rbind, lapply( 1:nrow(y), function (krow) {
                   if (krow<=xgrid$N) {
                       recomb_generator( yA=y[krow,], yB=y[krow+xgrid$N,], p=p$mid[krow], zeroind=rgrid$zeroind, nabvals=rgrid$N, eps=rgrid$eps )
                   } else {
                       recomb_generator( yA=y[krow,], yB=y[krow-xgrid$N,], p=1-p$mid[krow-xgrid$N], zeroind=rgrid$zeroind, nabvals=rgrid$N, eps=rgrid$eps )
                   }
             } ) )
        return( list( diffusion + recombination ) )
    }
    hap.soln <- ode.1D( y=yinit, times=times, func=rev.pde.fn, nspec=2*soln.dims[2], tcrit=max(fwds.soln[,1]) )
    attr(hap.soln,"r") <- cbind(left=avals,right=bvals)
    attr(hap.soln,"s") <- s
    attr(hap.soln,"sigma") <- sigma
    attr(hap.soln,"soln.dims") <- soln.dims
    attr(hap.soln,"grid") <- xgrid
    attr(hap.soln,"u") <- u
    return(hap.soln)
}

AD_from_soln <- function (soln,xgrid) {
    # returns a time x space matrix of ancestral disequilibrium values
    #  ** for only the first two loci in the solution **
    # i.e.,
    #   p u_A(r_0,r_0+r) + (1-p) u_B(r_0,r_0+r) - (p u_A(r_0) + (1-p) u_A(r_0))(p u_A(r_0+r) + (1-p) u_A(r_0+r))
    rcombs <- attr(soln,"r")
    this.rr <- sort( unique( as.vector(rcombs) ) )
    # must interpret the coordinates of the solution:
    #  factor saying if this coordinate is linked to A or to B
    ABfac <- rep( rep(c("A","B"),each=xgrid$N), nrow(rcombs) ) 
    #  factor saying which loci the coordinate covers
    rfac <- rep( paste( match(rcombs[,1],this.rr), match(rcombs[,2],this.rr), sep="" ), each=xgrid$N*2 )
    pp <- t( outer( xgrid$x.mid, soln[,1], attr(soln,"u") ) )
    return
    ( 
           pp * soln[,1+which(ABfac=="A"&rfac=="12")] + (1-pp) * soln[,1+which(ABfac=="B"&rfac=="12")] 
           - ( pp * soln[,1+which(ABfac=="A"&rfac=="11")] + (1-pp) * soln[,1+which(ABfac=="B"&rfac=="11")] )
               * ( pp * soln[,1+which(ABfac=="A"&rfac=="22")] + (1-pp) * soln[,1+which(ABfac=="B"&rfac=="22")] )
       )

}


get_hap_probs <- function ( hap.soln, left, right ) {
    # Return a deSolve-class object giving haplotype probabilites for the given (left,right) endpoints,
    # across time, space, and linked to A and B.
    nx <- attr(hap.soln,"soln.dims")[1]/2 # number of sites in spatial grid
    # find closest interval in rgrid (hack!)
    k.ab <- which.min( (attr(hap.soln,"r")[,"left"]-left)^2 + (attr(hap.soln,"r")[,"right"]-right)^2 )
    # columns of hap.soln are time, then ordered by space, linked allele, rows of r
    kk <- which( 1+(((1:ncol(hap.soln))-2)%/%(2*nx)) == k.ab )
    dummy.soln <- hap.soln[,c(1,kk)]
    class(dummy.soln) <- class(hap.soln)
    attr(dummy.soln,"r") <- attr(hap.soln,"r")[k.ab,]
    attr(dummy.soln,"nspec") <- 2
    attr(dummy.soln,"dimens") <- attr(hap.soln,"dimens")
    return(dummy.soln)
}

hap_cdf_to_pdf <- function ( hap.soln ) {
    # we compute probability that (a,b) is within the maximal segment
    # so need to convert to probability that (a,b) *is* the maximal segment:
    # prob(a,>b) = prob(<a,>b) - prob(<a-1,>b) )
    # prob(a,b) = prob(a,>b) - prob(a,>b+1)
    #             prob(<a,>b) - prob(<a-1,>b) ) - prob(<a,>b+1) + prob(<a-1,>b+1) )
    # except at boundaries: so, turn NAs below into zeros.
    r.coords <- attr(hap.soln,"r")
    rvals <- sort(unique(as.vector(r.coords)))
    segs.k <- cbind( left=match(r.coords[,"left"],rvals), right=match(r.coords[,"right"],rvals) )
    .hash <- function (a,b) { a*2*length(rvals)+b }  # factor of two avoids wraparound
    segs.neighbors <- cbind(
            left = match( .hash(segs.k[,1]-1,segs.k[,2]), .hash(segs.k[,1],segs.k[,2]) ),    # (a-1,b)
            right = match( .hash(segs.k[,1],segs.k[,2]+1), .hash(segs.k[,1],segs.k[,2]) ),   # (a,b+1)
            both = match( .hash(segs.k[,1]-1,segs.k[,2]+1), .hash(segs.k[,1],segs.k[,2]) ) ) # (a-1,b+1)
    # here's a nrow(r.coords) x nrow(r.coords) matrix to do the conversion
    cdf.to.pdf <- sparseMatrix( i = c((1:nrow(r.coords)),segs.neighbors[!is.na(segs.neighbors)]), 
                                j = rep((1:nrow(r.coords)),4)[c(rep(TRUE,nrow(r.coords)),!is.na(segs.neighbors))], 
                                x = rep(c(1,-1,-1,1),c(nrow(r.coords),colSums(!is.na(segs.neighbors)))),
                                dims = c(nrow(r.coords),nrow(r.coords)) )
    nxab <- attr(hap.soln,"soln.dims")[1]  # number of space x AB combinations
    for (xk in 1:nxab) {
        kk <- 1 + which( ( (0:(ncol(hap.soln)-2)) %% nxab )+1 == xk )
        hap.soln[,kk] <- as.matrix( hap.soln[,kk] %*% cdf.to.pdf )
    }
    return(hap.soln)
}

mean_hap_len <- function ( hap.pdf, loc ) {
    # Return a deSolve-class object giving mean haplotype lengths about location loc
    # across time, space, and linked to A and B.
    # MUST first convert to probabilities (not 1-cumulative probs)
    # i.e. do hap.pdf <- hap_cdf_to_pdf(hap.soln)
    nx <- attr(hap.pdf,"soln.dims")[1]/2 # number of sites in spatial grid
    r.coords <- attr(hap.pdf,"r")
    haplens <- r.coords[,"right"] - r.coords[,"left"]
    # first column of hap.pdf is time; after that columns are ordered by space, linked allele, rows of r
    # these are rows of r that loc falls in
    k.in <- which( ( loc >= r.coords[,"left"] ) & ( loc <= r.coords[,"right"] ) )
    k.ab <- k.in[1]
    kk <- which( 1+(((1:ncol(hap.pdf))-2)%/%(2*nx)) == k.ab )
    if (length(kk)==0) { warning("Location outside the range simulated."); return(NULL) }
    dummy.soln <- hap.pdf[,c(1,kk)]
    dummy.soln[,-1] <- dummy.soln[,-1] * haplens[k.in[1]]
    for (k.ab in k.in[-1]) {
        kk <- which( 1+(((1:ncol(hap.pdf))-2)%/%(2*nx)) == k.ab )
        dummy.soln[,-1] <- dummy.soln[,-1] + hap.pdf[,kk] * haplens[k.ab]
    }
    class(dummy.soln) <- class(hap.pdf)
    attr(dummy.soln,"nspec") <- 2
    attr(dummy.soln,"dimens") <- attr(hap.pdf,"dimens")
    return(dummy.soln)
}

## test unifprod

source("../unifprod.R",chdir=TRUE)


###
# Test correctness

    unifprod_R <- function (x,subtract=TRUE) {
        z <- matrix(0,nrow=NROW(x),ncol=NCOL(x))
        stopifnot( NCOL(x) == NROW(x) )
        for (i in seq(1,length.out=nrow(z)-1)) {
            for (k in seq(i+1,ncol(z))) {
                if (subtract) {
                    z[i,k] <- sum(x[i,seq(i,k-1)]*x[seq(i+1,k),k]-x[i,k])
                } else {
                    z[i,k] <- sum(x[i,seq(i,k-1)]*x[seq(i+1,k),k])
                }
            }
        }
        return(z)
    }

    unifprod_ut_R <- function (y,ncx,subtract=TRUE) {
        # same as unifprod_R, but where y = x[upper.tri(x,diag=TRUE)]
        z <- numeric(length(y))
        index <- 0    # where are we in the input vector
        for (j in seq(1,ncx)) {    # column
            for (i in seq(1,j)) {  # row
                index <- index + 1
                pindex <- (j+1)*(j+2)/2    # index of first partner at [j+1,j+1] in product
                outindex <- index + j     # where we are in the output (also equal to j*(j+1)/2 + i)
                if (subtract) { 
                    z[index] <- z[index] - (j-i)*y[index] 
                }
                if ( i<ncx && j<ncx ) {
                    for (k in (j+1):ncx) {  # column of product
                        # A[i,j] contributes to C[i,k] multiplied by A[j+1,k]
                        #   and if A[u,k] is at index m then A[u,k+1] is at position m+k
                        # outindex <- (k-1)*k/2 + i
                        z[outindex] <- z[outindex] + y[index]*y[pindex]
                        pindex <- pindex+k
                        outindex <- outindex+k
                    }
                }
            }
        }
        return(z)
    }

    biprod_R <- function (xA,xB,zeroind,subtract=TRUE) {
        z <- matrix(0,nrow=NROW(x),ncol=NCOL(x))
        stopifnot( NCOL(x) == NROW(x) )
        for (i in seq(1,length.out=nrow(z)-1)) {
            for (k in seq(i+1,ncol(z))) {
                j1 <- seq(i,k-1)
                j2 <- j1+1
                if (i>=zeroind) {  # both to the right of zero
                    xy <- xA[i,j1]*xB[j2,k]
                } else if (k<=zeroind) {  # both to the left of zero
                    xy <- xB[i,j1]*xA[j2,k]
                } else {  # i <= 0 <= k
                    xy <- ifelse(j1>=zeroind, xA[i,j1]*xB[j2,k], xB[i,j1]*xA[j2,k] )
                }
                if (subtract) {
                    z[i,k] <- sum(xy-xA[i,k])
                } else {
                    z[i,k] <- sum(xy)
                }
            }
        }
        return(z)
    }


    biprod_ut_R <- function (yA,yB,zeroind,ncx,subtract=TRUE) {
        # discretization of  y -> z where
        #   z[a,b] = \int_a^0 ( yA(a,x) yB(x,b) - yA(a,b) dx + \int_0^b ( yB(a,x) yA(x,b) - yA(a,b) dx 
        # to discretize this we sum over recombinations occuring between marker (k-1) and marker k
        #   and zeroind is the index of the marker at zero
        #   so whichever interval has zeroind in it should be yA
        z <- numeric(length(yA))
        index <- 0    # where are we in the input vector
        for (j in seq(1,ncx)) {    # column
            for (i in seq(1,j)) {  # row
                index <- index + 1
                pindex <- (j+1)*(j+2)/2    # index of first partner at [j+1,j+1] in product
                outindex <- index + j     # where we are in the output (also equal to j*(j+1)/2 + i)
                if (subtract) { 
                    # cat("right: i=", i, ", j=", j, ", index=", index, "subtr=", (j-i)*yA[index], "\n")
                    z[index] <- z[index] - (j-i)*yA[index] 
                }
                if ( i<ncx && j<ncx ) {
                    for (k in (j+1):ncx) {  # column of product
                        # A[i,j] contributes to C[i,k] multiplied by A[j+1,k]
                        #   and if A[u,k] is at index m then A[u,k+1] is at position m+k
                        # outindex <- (k-1)*k/2 + i
                        if (j<zeroind) { # right-hand interval contains zero
                            # cat("right: i=", i, ", j=", j, ", k=", k, "\n")
                            z[outindex] <- z[outindex] + yB[index]*yA[pindex]
                        } else {
                            # cat(" left: i=", i, ", j=", j, ", k=", k, "\n")
                            z[outindex] <- z[outindex] + yA[index]*yB[pindex]
                        }
                        pindex <- pindex+k
                        outindex <- outindex+k
                    }
                }
            }
        }
        return(z)
    }


    # must run `R CMD SHLIB unifprod.c` first
    dyn.load("unifprod.so")
    unifprod2 <- function (x) {
        .Call("unifprod2",as.numeric(x),nrow(x),ncol(x))
    }
    unifprod_ut2 <- function (x,ncx) {
        .Call("unifprod_ut2",as.numeric(x),ncx)
    }

    x <- matrix( 1:9, nrow=3 )
    xut <- x[upper.tri(x,diag=TRUE)]
    stopifnot( all.equal( unifprod_ut_R(xut,3,subtract=FALSE), c(0, 1*5, 0, 1*8+4*9, 5*9, 0) ) ) 
    stopifnot( all.equal( unifprod_ut_R(xut,3), c(0, 1*5-4, 0, 1*8+4*9-2*7, 5*9-8, 0) ) ) 

    y <- 9+x
    yut <- 9+xut
    stopifnot( all.equal( 
                         biprod_R(xA=x,xB=y,zeroind=2,subtract=FALSE), 
                         matrix( c(0, 0, 0, 10*5, 0, 0, 10*8 + 4*18, 5*18, 0), nrow=3 ) 
                     ) )
    stopifnot( all.equal( 
                         biprod_R(xA=x,xB=y,zeroind=2),
                         matrix( c(0, 0, 0, 10*5 - 4, 0, 0, 10*8 + 4*18 - 2*7, 5*18 - 8, 0), nrow=3 ) 
                     ) )
    stopifnot( all.equal( 
                         biprod_ut_R(yA=xut,yB=yut,zeroind=2,ncx=3,subtract=FALSE), 
                         c(0, 10*5, 0, 10*8 + 4*18, 5*18, 0) 
                     ) )
    stopifnot( all.equal( 
                         biprod_ut_R(yA=xut,yB=yut,zeroind=2,ncx=3),
                         c(0, 10*5-4, 0, 10*8+4*18-2*7, 5*18-8, 0) 
                     ) ) 
    stopifnot( all.equal( 
                         biprod_ut_R(yA=xut,yB=yut,zeroind=2,ncx=3),
                         biprod_ut(yA=xut,yB=yut,zeroind=2,ncx=3),
                     ) )


    for (x in list( matrix(1:9,nrow=3), matrix(rnorm(1e4),nrow=100) ) ) {
        y <- matrix(rnorm(length(x)),nrow=nrow(x))
        zind <- floor(nrow(x)*0.75)
        stopifnot( 
            all.equal( unifprod(x) , unifprod2(x) )
            && all.equal( unifprod(x) , unifprod_R(x) )
            && all.equal( unifprod2(x) , unifprod_R(x) )
            && all.equal( unifprod_R(x,subtract=FALSE)[upper.tri(x,diag=TRUE)], unifprod_ut_R(x[upper.tri(x,diag=TRUE)],ncol(x),subtract=FALSE) )
            && all.equal( unifprod_R(x)[upper.tri(x,diag=TRUE)], unifprod_ut_R(x[upper.tri(x,diag=TRUE)],ncol(x)) )
            && all.equal( 
                    unifprod_ut_R(x[upper.tri(x,diag=TRUE)],ncol(x)), 
                    unifprod_ut(x[upper.tri(x,diag=TRUE)],ncol(x)) 
                )
            && all.equal( 
                    unifprod_ut_R(x[upper.tri(x,diag=TRUE)],ncol(x)), 
                    unifprod_ut2(x[upper.tri(x,diag=TRUE)],ncol(x)) 
                )
            && all.equal(
                    biprod_ut_R(yA=x[upper.tri(x,diag=TRUE)],yB=y[upper.tri(x,diag=TRUE)],zeroind=zind,ncx=ncol(x)),
                    biprod_ut(yA=x[upper.tri(x,diag=TRUE)],yB=y[upper.tri(x,diag=TRUE)],zeroind=zind,ncx=ncol(x))
                )
        )
    }

    x <- matrix(rnorm(1e4),nrow=100)
    xut <- x[upper.tri(x,diag=TRUE)]


###
# Test speed, comparing different versions

    library(microbenchmark)
    microbenchmark(unifprod_R(x),times=100)
    microbenchmark(unifprod(x),times=1000)
    microbenchmark(unifprod2(x),times=1000)
    microbenchmark(unifprod_ut_R(xut,ncol(x)),times=100)
    microbenchmark(unifprod_ut(xut,ncol(x)),times=1000)
    microbenchmark(unifprod_ut2(xut,ncol(x)),times=1000)


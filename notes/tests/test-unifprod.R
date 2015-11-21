## test unifprod

source("unifprod.R")

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

    # must run `R CMD SHLIB unifprod.c` first
    dyn.load("unifprod.so")
    unifprod2 <- function (x) {
        .Call("unifprod2",as.numeric(x),nrow(x),ncol(x))
    }
    unifprod_ut2 <- function (x,ncx) {
        .Call("unifprod_ut2",as.numeric(x),ncx)
    }

    for (x in list( matrix(1:9,nrow=3), matrix(rnorm(1e4),nrow=100) ) ) {
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
        )
    }

    x <- matrix(rnorm(1e4),nrow=100)
    xut <- x[upper.tri(x,diag=TRUE)]

    library(microbenchmark)
    microbenchmark(unifprod_R(x),times=100)
    microbenchmark(unifprod(x),times=1000)
    microbenchmark(unifprod2(x),times=1000)
    microbenchmark(unifprod_ut_R(xut,ncol(x)),times=100)
    microbenchmark(unifprod_ut(xut,ncol(x)),times=1000)
    microbenchmark(unifprod_ut2(xut,ncol(x)),times=1000)


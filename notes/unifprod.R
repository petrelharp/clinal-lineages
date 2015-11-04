library(Rcpp)

# Discretization of u(a,b) -> \int_a^b (u(a,t) u(t,b) - u(a,b)) dt
#  so output is y[i,k] = sum( x[i,i:(k-1)]*x[(i+1):k,k] )
#
# The key step is to do something like:
#     1 | x y z
#   a 2 | . . u
#     3 | . . v
#       +------
#     b   1 2 3
# --> should  return x*u + y*v
# the inner product of a partial row with a partial column

unifprod.cpp <- 'NumericMatrix unifprod(NumericMatrix x) {
    int nrx = x.nrow(), ncx = x.ncol();
    double sum, xik;
    NumericMatrix z(nrx,nrx);

    if (nrx != ncx) stop("Matrix must be square.");

    if (nrx > 0 && ncx > 0) {
        for (int i = 0; i < nrx; i++)
        for (int k = 0; k < ncx; k++) {
            sum = 0.0;
            xik = x(i,k);
            for (int j = i; j < k; j++) {
                sum += x(i , j) * x(j + 1 , k) - x(i, k);
            }
            z(i , k) = sum;
        }
    } else { /* zero-extent operations should return zeroes */
        for(int i = 0; i < nrx*ncx; i++) z[i] = 0;
    }
    return z;
}'

cppFunction(unifprod.cpp)


if (FALSE) { ## testing
    unifprod_R <- function (x) {
        z <- matrix(0,nrow=NROW(x),ncol=NCOL(x))
        stopifnot( NCOL(x) == NROW(x) )
        for (i in seq(1,length.out=nrow(z)-1)) {
            for (k in seq(i+1,ncol(z))) {
                z[i,k] <- sum(x[i,seq(i,k-1)]*x[seq(i+1,k),k]-x[i,k])
            }
        }
        return(z)
    }

    unifprod_ut_R <- function (x,ncx) {
        # same, but for upper triangular ordered x
        z <- matrix(0,nrow=ncx,ncol=ncx)
        index <- 0  # where are we in the vector
        for (j in seq(1,ncx)) {    # column
            for (i in seq(1,j)) {  # row
                index <- index + 1
                pindex <- index + 1      # index of partner in product
                # z[i,j] <- z[i,j] - (col(z)-row(z))[i,j]*x[index]
                if ( i<ncx && j<ncx ) {
                    for (k in j:(ncx-1)) {  # column of product - 1
                        # A[i,j] contributes to C[i,k] multiplied by A[j+1,k]
                        #   and the next number to the right is at position +j
                        pindex <- pindex+k
                        z[i,k+1] <- z[i,k+1] + x[index]*x[pindex]
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

    for (x in list( matrix(1:9,nrow=3), matrix(rnorm(1e4),nrow=100) ) ) {
        stopifnot( 
            all.equal( unifprod(x) , unifprod2(x) )
            && all.equal( unifprod(x) , unifprod_R(x) )
            && all.equal( unifprod2(x) , unifprod_R(x) )
        )
    }
    library(microbenchmark)
    microbenchmark(unifprod_R(x),times=100)
    microbenchmark(unifprod(x),times=1000)
    microbenchmark(unifprod2(x),times=1000)

    all.equal( unifprod_R(x)[upper.tri(x,diag=TRUE)], unifprod_ut_R(x[upper.tri(x,diag=TRUE)],ncol(x)) )

}

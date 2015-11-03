library(Rcpp)

unifprod.cpp <- 'NumericMatrix unifprod(NumericMatrix x) {
    int nrx = x.nrow(), ncx = x.ncol();
    double sum;
    NumericMatrix z(nrx,nrx);

    if (nrx != ncx) stop("Matrix must be square.");

    if (nrx > 0 && ncx > 0) {
        for (int i = 0; i < nrx; i++)
        for (int k = 0; k < ncx; k++) {
            sum = 0.0;
            if (i<k)
            for (int j = i; j < k; j++) {
                sum += x(i , j) * x(j + 1 , k);
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
                z[i,k] <- sum(x[i,seq(i,k-1)]*x[seq(i+1,k),k])
            }
        }
        return(z)
    }

    # must run `R CMD SHLIB unifprod2.c` first
    dyn.load("unifprod2.so")
    unifprod2 <- function (x) {
        .Call("unifprod2",x,nrow(x),ncol(x))
    }

    x <- matrix(rnorm(1e4),nrow=100)
    all( unifprod(x) == unifprod2(x) )
    microbenchmark(unifprod_R(x),times=100)
    microbenchmark(unifprod(x),times=1000)
    microbenchmark(unifprod2(x),times=1000)
}

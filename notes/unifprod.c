#include <R.h>
#include <Rinternals.h>

/* modified from ./src/main/array.c in R-3.2.2 source */

/*
 dyn.load("unifprod.so")
 unifprod <- function (x,y) { .C("unifprod",x=as.numeric(x),nrx=NROW(x),ncx=NCOL(x),y=as.numeric(y),nry=NROW(y),ncy=NCOL(y),z=matrix(0.0,nrow=NROW(x),ncol=NCOL(y)))[[7]] }
 unifprod(diag(3)*1.0,diag(3)*1.0)
 */

SEXP unifprod2(SEXP x, SEXP Rnrx, SEXP Rncx)
{
    /* compute the matrix product (A,B) -> C where
     *   C[i,j] = A[i]*B[i+1] + A[i+1]*B[i+2] + ... + A[j-1]*B[j]
     * and so if A is (m x n) and B is (n x o) then C is (m x o) and
     *   C[i,i] = 0
     */
    double sum;
    int nrx = INTEGER(Rnrx)[0];
    int ncx = INTEGER(Rncx)[0];
    double *X = REAL(x);
    SEXP z;
    PROTECT(z = allocMatrix(REALSXP, nrx, ncx));

    if (nrx > 0 && ncx > 0) {
        for (int i = 0; i < nrx; i++)
        for (int k = 0; k < ncx; k++) {
            sum = 0.0;
            if (i<k)
            for (int j = i; j < k; j++) {
                sum += X[i + j * nrx] * X[j + 1 + k * nrx];
            }
            REAL(z)[i + k * nrx] = (double) sum;
        }
    } else { /* zero-extent operations should return zeroes */
        for(int i = 0; i < nrx*ncx; i++) REAL(z)[i] = 0.0;
    }
    UNPROTECT(1);
    return z;
}

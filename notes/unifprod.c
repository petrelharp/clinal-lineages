#include <R.h>
#include <Rinternals.h>

/* modified from ./src/main/array.c in R-3.2.2 source */

SEXP unifprod2(SEXP x, SEXP Rnrx, SEXP Rncx)
{
    /* compute the matrix transformation A -> C where
     *   C[i,j] = A[i,i]*A[i+1,j] + A[i,i+1]*A[i+2,j] + ... + A[i,j-1]*A[j,j] - (j-i) A[i,j]
     *   C[i,i] = 0
     */
    double sum, xik;
    int nrx = INTEGER(Rnrx)[0];
    int ncx = INTEGER(Rncx)[0];
    double *X = REAL(x);
    SEXP z;
    PROTECT(z = allocMatrix(REALSXP, nrx, ncx));

    if (nrx > 0 && ncx > 0) {
        for (int i = 0; i < nrx; i++)
        for (int k = 0; k < ncx; k++) {
            sum = 0.0;
            xik = X[i + k * nrx];
            for (int j = i; j < k; j++) {
                sum += X[i + j * nrx] * X[j + 1 + k * nrx] - xik;
            }
            REAL(z)[i + k * nrx] = (double) sum;
            REAL(z)[i + k * nrx] = (double) sum;
        }
    } else { /* zero-extent operations should return zeroes */
        for(int i = 0; i < nrx*ncx; i++) REAL(z)[i] = 0.0;
    }
    UNPROTECT(1);
    return z;
}

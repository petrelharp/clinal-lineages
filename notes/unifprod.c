#include <R.h>
#include <Rinternals.h>

/* modified from ./src/main/array.c in R-3.2.2 source */

/*
 dyn.load("unifprod.so")
 unifprod <- function (x,y) { .C("unifprod",x=as.numeric(x),nrx=NROW(x),ncx=NCOL(x),y=as.numeric(y),nry=NROW(y),ncy=NCOL(y),z=matrix(0.0,nrow=NROW(x),ncol=NCOL(y)))[[7]] }
 unifprod(diag(3)*1.0,diag(3)*1.0)
 */

void unifprod(double *x, int *nrx, int *ncx,
                    double *y, int *nry, int *ncy, double *z)
{
    /* compute the matrix product (A,B) -> C where
     *   C[i,j] = A[i]*B[i+1] + A[i+1]*B[i+2] + ... + A[j-1]*B[j]
     * and so if A is (m x n) and B is (n x o) then C is (m x o) and
     *   C[i,i] = 0
     */
    double sum;
    R_xlen_t NRX = *nrx, NRY = *nry;
    R_xlen_t NCX = *ncx, NCY = *ncy;

    if (NRX > 0 && NCX > 0 && NRY > 0 && NCY > 0) {
        for (int i = 0; i < NRX; i++)
        for (int k = 0; k < NCY; k++) {
            sum = 0.0;
            for (int j = i; j < NCX-1; j++) {
                sum += x[i + j * NRX] * y[j + 1 + k * NRY];
            }
            z[i + k * NRX] = (double) sum;
        }
    } else { /* zero-extent operations should return zeroes */
        for(R_xlen_t i = 0; i < NRX*NCY; i++) z[i] = 0;
    }
}


// // "ORIGINAL"
// static void unifprod(double *x, int nrx, int ncx,
//                     double *y, int nry, int ncy, double *z)
// {
//     /* compute the matrix product (A,B) -> C where
//      *   C[i,j] = A[i]*B[i+1] + A[i+1]*B[i+2] + ... + A[j-1]*B[j]
//      * and so if A is (m x n) and B is (n x o) then C is (m x o) and
//      *   C[i,i] = 0
//      */
//     double sum;
//     R_xlen_t NRX = nrx, NRY = nry;
// 
//     if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
//         for (int i = 0; i < nrx; i++)
//         for (int k = 0; k < ncy; k++) {
//             sum = 0.0;
//             for (int j = i; j < ncx-1; j++)
//             sum += x[i + j * NRX] * y[j + 1 + k * NRY];
//             z[i + k * NRX] = (double) sum;
//         }
//     } else { /* zero-extent operations should return zeroes */
//         for(R_xlen_t i = 0; i < NRX*ncy; i++) z[i] = 0;
//     }
// }

// SEXP unifprod(double *x, int nrx, int ncx,
//                     double *y, int nry, int ncy)
// {
//     /* compute the matrix product (A,B) -> C where
//      *   C[i,j] = A[i]*B[i+1] + A[i+1]*B[i+2] + ... + A[j-1]*B[j]
//      * and so if A is (m x n) and B is (n x o) then C is (m x o) and
//      *   C[i,i] = 0
//      */
//     double sum;
//     R_xlen_t NRX = asInteger(nrx), NRY = asInteger(nry);
//     R_xlen_t NCX = asInteger(ncx), NCY = asInteger(ncy);
//     SEXP z = PROTECT(allocVector(REALSXP, NRX*NCY));
// 
//     if (NRX > 0 && NCX > 0 && NRY > 0 && NCY > 0) {
//         for (int i = 0; i < NRX; i++)
//         for (int k = 0; k < NCY; k++) {
//             sum = 0.0;
//             for (int j = i; j < NCX-1; j++)
//             sum += x[i + j * NRX] * y[j + 1 + k * NRY];
//             z[i + k * NRX] = (double) sum;
//         }
//     } else { /* zero-extent operations should return zeroes */
//         for(R_xlen_t i = 0; i < NRX*NCY; i++) z[i] = 0;
//     }
//     return z;
// }

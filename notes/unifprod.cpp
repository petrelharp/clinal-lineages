NumericMatrix unifprod(NumericMatrix x) {
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
}


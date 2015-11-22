NumericVector biprod_ut(
        NumericVector yA, 
        NumericVector yB,
        int zeroind,
        int ncx) {
    NumericVector z(yA.size());
    int index = 0;
    int pindex, outindex;
    for (int i=0; i<yA.size(); i++) {
        z[i] = 0.0;
    }
    for (int j=1; j<=ncx; j++) {
        for (int i=1; i<=j; i++) {
            pindex = (j+1)*(j+2)/2 - 1;
            outindex = index + j;
            z[outindex-j] -= (j-i)*yA[index];
            // printf("i=%i, j=%i, index=%i, yi=%f \\n",i,j,index,yA[index]);
            if ( (i<ncx) && (j<ncx) ) {
                for (int k=j+1; k<=ncx; k++) {
                    // printf("i=%i, j=%i, k=%i, index=%i, pindex=%i, outindex=%i yy=%f \\n",i,j,k,index+1,pindex+1,outindex+1,yA[index]*yA[pindex]);
                    if ( j < zeroind ) {
                        z[outindex] += yB[index]*yA[pindex];
                    } else {
                        z[outindex] += yA[index]*yB[pindex];
                    }
                    pindex += k;
                    outindex += k;
                }
            }
            index++;
        }
    }
    return z;
}



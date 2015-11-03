dyn.load("unifprod.so")
unifprod <- function (x,y) { 
    .C("unifprod",
            x=as.numeric(x),
            nrx=NROW(x),
            ncx=NCOL(x),
            y=as.numeric(y),
            nry=NROW(y),
            ncy=NCOL(y),
            z=matrix(0.0,nrow=NROW(x),ncol=NCOL(y))
        )[[7]] 
}
unifprod(diag(3)*1.0,diag(3)*1.0)


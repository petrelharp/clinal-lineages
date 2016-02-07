pdf('cline_centers.pdf', width=6,height=4, pointsize=10)

for (x in list.files(".","simul.*_chunks.Robj")) {

    load(x)

    B_freq = do.call(cbind,lapply(chunks,function(C){sapply(C,function(X){length(which(X>0))})}))


    matplot(B_freq/1000,xlab="geo.pos",ylab="p_B", type='l', lty=1, main=x, cex.main=0.5)
    matlines(1-B_freq/1000, lty=2)
    abline(v=25.5)
}

dev.off()

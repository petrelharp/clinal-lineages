load("~/Downloads/number_of_ancestors_tau1000.Robj")

number_of_founders_per_ind_neutral = do.call(cbind,lapply(number_of_ancestors[[1]],function(X){sapply(X,function(Z){length(unique(Z))/length(Z)})
}))
number_of_founders_per_ind_unlinked = do.call(cbind,lapply(number_of_ancestors[[2]],function(X){sapply(X,function(Z){length(unique(Z))/length(Z)})
}))
number_of_founders_per_ind_linked = do.call(cbind,lapply(number_of_ancestors[[3]],function(X){sapply(X,function(Z){length(unique(Z))/length(Z)})
}))

r2y=colorRampPalette(c("red","yellow"))

pdf(file="fmaily_size_tau1000.pdf_revision.pdf",width=6,height=5)
par(mar=c(3,3.5,0.5,0.5))
matplot(1/number_of_founders_per_ind_neutral[,seq(20,61,2)],type="l",lty=1,col="grey",xaxt="n",yaxt="n",ylab="",xlab="",ylim=c(0,30))
matpoints(1/number_of_founders_per_ind_unlinked[,seq(20,61,2)],type="l",lty=1,col="yellow",xaxt="n",yaxt="n",ylab="",xlab="")
points(1/number_of_founders_per_ind_unlinked[,41],type="l",col="yellow",lwd=5)

matpoints(1/number_of_founders_per_ind_linked[,seq(20,61,2)],type="l",lty=1,col=c(rev(r2y(15)[3:12]),"red",r2y(15)[3:12]),lwd=1.5,xaxt="n",yaxt="n",ylab="",xlab="")
points(1/number_of_founders_per_ind_linked[,41],type="l",col="red",lwd=5)


legend("bottomright",legend = c("neutral","unlinked","linked"),col=c("grey","yellow","red"),lty=1,horiz=TRUE)
axis(1,at=c(seq(0.5,20.5,5),25.5,seq(30.5,50.5,5)),labels=c(seq(-25,-5,5),0,seq(5,25,5)),cex.axis=0.8)
mtext("Distance from HZ center",side=1,line=2)
axis(2,at=seq(0,30,5),cex.axis=0.8,las=2)
mtext("mean number of descendents per haplotype",side=2,line=2.5)
dev.off()


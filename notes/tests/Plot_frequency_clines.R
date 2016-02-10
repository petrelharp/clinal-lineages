demeID=do.call(c,lapply(1:ndemes,function(X){rep(X,500)}))
focal_sites = lapply(qtl,function(Z){Z$pos})
positions = seq(0.2,0.8,0.01)

load("sims/simulation_SIGMA1_Ninds25000_ndemes50_s0.1_tau100_run2016-01-11_simsums.Robj")
#CHUNKS_SEL = get.chunks.at.positions(positions,CHR=1)
B_freq = get.ancestry.freqs(positions,CHR=1)


load("sims/simulation_SIGMA1_Ninds25000_ndemes50_s0_tau100_run2016-01-15_simsums.Robj")
#CHUNKS_NEU = get.chunks.at.positions(positions,CHR=2)
B_neu = get.ancestry.freqs(positions,CHR=2)


pdf(file="drafts/figs/alleleFrequencies_sim.pdf",height=4,width=6)
par(mar=c(3.5,3.5,0.5,0.5))
matplot(B_neu/1000,xlab="",ylab="",type="l",lty=1,col="grey",lwd=0.5,x=seq(-24.5,24.6,1))
for(i in 1:30){
	matpoints(B_freq[,c(i,62-i)]/1000,type="l",col=(rainbow(43))[32-i],lty=1,lwd=0.75,x=seq(-24.5,24.6,1))
}
points(B_freq[,31]/1000,col="red",lwd=2,type="l",x=seq(-24.5,24.6,1))
mtext("Geographic position",side=1,line=2)
mtext(expression(p[B]),side=2,line=2.5)

legend("bottomright",legend = c("r=0","r=0.05","r=0.1","r=0.2","r=0.3","no seln."),col=c("red",rainbow(43)[c(6,11,21,31)],"darkgrey"),lwd=c(2,rep(0.75,4),0.5),cex=0.8)

dev.off()


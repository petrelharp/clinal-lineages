load("simulation_SIGMA1_Ninds25000_ndemes50_s0.1_tau1000_chr1_start0.48_stop0.52_by0.001_intervalSizes_allAncs.Robj")
load("simulation_SIGMA1_Ninds25000_ndemes50_s0.1_tau1000_chr1_start0.48_stop0.52_by0.001_intervalSizes.Robj")

Ndemes=50
deme_ID = do.call(c,lapply(1:Ndemes,function(X){rep(X,Ninds/Ndemes)}))

positions=seq(0.48,0.52,0.001)

pdf(height=4, width=6.25,file="blocksAlongChromNoConditioning.pdf")
par(mar=c(3.5,3.5,0.5,0.5))
mean_per_deme = do.call(rbind, lapply(intervalSizes_allAncs,function(P){tapply(1:25000,deme_ID,function(Z){mean(P[Z,])})}))
plot(mean_per_deme[,1]/mean(mean_per_deme[,1],na.rm=T),x=positions,ylim=c(0.5,3),type="l",lty=1,col=c(rainbow(25),rev(rainbow(25))),main="",cex.main=0.8,xlab="",ylab="",yaxt="n",xaxt="n")
for(i in 1:50){points(mean_per_deme[,i]/mean(mean_per_deme[,i],na.rm=T),x=positions,type="l",col=rainbow(50)[i])}
axis(1,at=seq(0.48,0.52,0.01),labels=seq(-0.02,0.02,0.01),cex.axis=0.9)
axis(2,at=seq(0.5,3,0.5),cex.axis=0.9,las=2)
mtext("Distance from selected locus (M)",side=1,line=2.5)
mtext("Normalized mean block length",side=2,line=2.5)
legend("topright",title="distance from center:",legend = rev(c(seq(-25,-5,5),seq(5,25,5))),cex=0.7,col=rainbow(50)[seq(1,50,5)],lty=1,xjust=1)
dev.off()

pdf(height=4, width=6.25,file="blocksAlongChromAncBConditioning.pdf")
par(mar=c(3.5,3.5,0.5,0.5))

mean_per_deme_AncB = do.call(rbind, lapply(intervalSizes,function(P){tapply(1:25000,deme_ID,function(Z){all_sites = P[Z,]; mean(all_sites[which(all_sites>0)])})}))
plot(mean_per_deme_AncB[,1]/mean(mean_per_deme_AncB[,1],na.rm=T),x=positions,ylim=c(0.5,3),type="l",lty=1,col=c(rainbow(25),rev(rainbow(25))),main="",cex.main=0.8,xlab="",ylab="",yaxt="n",xaxt="n")
for(i in 1:50){points(mean_per_deme_AncB[,i]/mean(mean_per_deme_AncB[,i],na.rm=T),x=positions,type="l",col=rainbow(50)[i])}
axis(1,at=seq(0.48,0.52,0.01),labels=seq(-0.02,0.02,0.01),cex.axis=0.9)
axis(2,at=seq(0.5,3,0.5),cex.axis=0.9,las=2)
mtext("Distance from selected locus (M)",side=1,line=2.5)
mtext("Normalized mean block length",side=2,line=2.5)
legend("topright",title="distance from center:",legend = rev(c(seq(-25,-5,5),seq(5,25,5))),cex=0.7,col=rainbow(50)[seq(1,50,5)],lty=1,xjust=1)
dev.off()


load("drafts/figs/simulation_SIGMA1_Ninds25000_ndemes50_s0_tau1000_chr1_start0.48_stop0.52_by0.001_intervalSizes_allAncs.Robj")
load("drafts/figs/simulation_SIGMA1_Ninds25000_ndemes50_s0_tau1000_chr1_start0.48_stop0.52_by0.001_intervalSizes.Robj")

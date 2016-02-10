load("simulation_SIGMA1_Ninds25000_ndemes50_s0.1_tau1000_chr1_start0.48_stop0.52_by0.001_intervalSizes_allAncs.Robj")
load("simulation_SIGMA1_Ninds25000_ndemes50_s0.1_tau1000_chr1_start0.48_stop0.52_by0.001_intervalSizes.Robj")

Ndemes=50
deme_ID = do.call(c,lapply(1:Ndemes,function(X){rep(X,Ninds/Ndemes)}))

positions=seq(0.48,0.52,0.001)

pdf(height=4, width=6.25,file="blocksAlongChromNoConditioning.pdf")
par(mar=c(3.5,3.5,0.5,0.5))
mean_per_deme = do.call(rbind, lapply(intervalSizes_allAncs,function(P){tapply(1:25000,deme_ID,function(Z){mean(P[Z,])})}))
plot(mean_per_deme[,1]/mean(mean_per_deme[,1],na.rm=T),x=positions,ylim=c(0.5,3),type="l",lty=1,col="white",main="",cex.main=0.8,xlab="",ylab="",yaxt="n",xaxt="n")
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
plot(mean_per_deme_AncB[,1]/mean(mean_per_deme_AncB[,1],na.rm=T),x=positions,ylim=c(0.5,3),type="l",lty=1,col="white",main="",cex.main=0.8,xlab="",ylab="",yaxt="n",xaxt="n")
for(i in 1:50){points(mean_per_deme_AncB[,i]/mean(mean_per_deme_AncB[,i],na.rm=T),x=positions,type="l",col=rainbow(50)[i])}
axis(1,at=seq(0.48,0.52,0.01),labels=seq(-0.02,0.02,0.01),cex.axis=0.9)
axis(2,at=seq(0.5,3,0.5),cex.axis=0.9,las=2)
mtext("Distance from selected locus (M)",side=1,line=2.5)
mtext("Normalized mean block length",side=2,line=2.5)
legend("topright",title="distance from center:",legend = rev(c(seq(-25,-5,5),seq(5,25,5))),cex=0.7,col=rainbow(50)[seq(1,50,5)],lty=1,xjust=1)
dev.off()

pdf(height=4, width=6.25,file="blocksAlongChromHeatmap.pdf")

layout(c(1,2),heights=c(1,4))
par(mar=c(0.5,3.5,2,0.5))
	BREAKS = seq(0.004,0.013,0.0005)
	COLORS = heat.colors(length(BREAKS)-1)
	plot(0,0,col="white",xlim=range(BREAKS),ylim=c(0,1),yaxt="n",ylab="",xlab="",bty="n",xaxt="n")
	polygon(c(BREAKS[1],BREAKS[1],BREAKS[length(BREAKS)],BREAKS[length(BREAKS)]),y=c(0,1,1,0))
	for(i in 1:(length(BREAKS)-1)){
		polygon(c(BREAKS[i],BREAKS[i],BREAKS[i+1],BREAKS[i+1]),y=c(0,1,1,0),col=COLORS[i],border=COLORS[i])	
	}	

	axis(line=,side=3,(BREAKS[-1]+BREAKS[1:(length(BREAKS)-1)])/2,labels=(BREAKS[-1]+BREAKS[1:(length(BREAKS)-1)])/2*100,cex.axis=0.7,lwd=0,lwd.ticks=1)

	par(mar=c(3.5,3.5,0,0.5))
	image(mean_per_deme[,rev(1:50)],x=positions,y=-25:25,main="",yaxt="n",xaxt="n",ylab="",xlab="",breaks=BREAKS,col=COLORS)
	axis(1,at=seq(0.48,0.52,0.01),labels=seq(-0.02,0.02,0.01),cex.axis=0.9)
	axis(2,at=seq(-24.5,0,5),labels=seq(-25,-5,5),cex.axis=0.9,las=2)
	axis(2,at=seq(0,24.5,5),labels=seq(5,25,5),cex.axis=0.9,las=2)
	mtext("Distance from selected locus (M)",side=1,line=2.5)
	mtext("Distance from HZ center",side=2,line=2.5)

dev.off()

pdf(height=4, width=6.25,file="blocksAlongChromHeatmapAncBConditioning.pdf")

layout(c(1,2),heights=c(1,4))
par(mar=c(0.5,3.5,2,0.5))

	BREAKS = seq(0.0015,0.0145,0.0005)
	COLORS = heat.colors(length(BREAKS)-1)

	plot(0,0,col="white",xlim=range(BREAKS),ylim=c(0,1),yaxt="n",ylab="",xlab="",bty="n",xaxt="n")
	polygon(c(BREAKS[1],BREAKS[1],BREAKS[length(BREAKS)],BREAKS[length(BREAKS)]),y=c(0,1,1,0))
	for(i in 1:(length(BREAKS)-1)){
		polygon(c(BREAKS[i],BREAKS[i],BREAKS[i+1],BREAKS[i+1]),y=c(0,1,1,0),col=COLORS[i],border=COLORS[i])	
	}	

	axis(line=,side=3,(BREAKS[-1]+BREAKS[1:(length(BREAKS)-1)])/2,labels=(BREAKS[-1]+BREAKS[1:(length(BREAKS)-1)])/2*100,cex.axis=0.7,lwd=0,lwd.ticks=1)

	par(mar=c(3.5,3.5,0,0.5))
	image(mean_per_deme_AncB[,50:1],x=positions,y=-25:25,main="",yaxt="n",xaxt="n",ylab="",xlab="",breaks=BREAKS,col=heat.colors(length(BREAKS)-1))
	axis(1,at=seq(0.48,0.52,0.01),labels=seq(-0.02,0.02,0.01),cex.axis=0.9)
	axis(2,at=seq(-24.5,0,5),labels=seq(-25,-5,5),cex.axis=0.9,las=2)
	axis(2,at=seq(0,24.5,5),labels=seq(5,25,5),cex.axis=0.9,las=2)
	mtext("Distance from selected locus (M)",side=1,line=2.5)
	mtext("Distance from HZ center",side=2,line=2.5)
dev.off()


load("simulation_SIGMA1_Ninds25000_ndemes50_s0.01_tau1000_chr1_start0.48_stop0.52_by0.001_flanking.blocks.by.ind.Robj")
load("simulation_SIGMA1_Ninds25000_ndemes50_s0.01_tau1000_chr1_start0.48_stop0.52_by0.001_flanking.blocks.by.ind.AncB.Robj")

pdf(height=4, width=6.25,file="adjacentBlocksAlongChromNoConditioning.pdf")
par(mar=c(3.5,3.5,0.5,0.5))
flanking.blocks.deme.matrix = do.call(rbind,lapply(1:50,function(DEME){sapply(flanking.blocks.by.ind,function(POS){mean(POS[[DEME]])})}))
matplot(t(flanking.blocks.deme.matrix/apply(flanking.blocks.deme.matrix,1,mean))[,50:1],xaxt="n",yaxt="n",xlab="",ylab="",lty=1,col=rainbow(50),type="l",x=positions,ylim=c(0.5,3))
axis(1,at=seq(0.48,0.52,0.01),labels=seq(-0.02,0.02,0.01),cex.axis=0.9)
axis(2,at=seq(0.5,3,0.5),cex.axis=0.9,las=2)
mtext("Distance from selected locus (M)",side=1,line=2.5)
mtext("Normalized mean block length",side=2,line=2.5)
legend("topright",title="distance from center:",legend = rev(c(seq(-25,-5,5),seq(5,25,5))),cex=0.7,col=rainbow(50)[seq(1,50,5)],lty=1,xjust=1)
dev.off()

pdf(height=4, width=6.25,file="adjacentBlocksAlongChromAncBConditioning.pdf")
par(mar=c(3.5,3.5,0.5,0.5))
flanking.blocks.deme.matrix.AncB = do.call(rbind,lapply(1:50,function(DEME){sapply(flanking.blocks.by.ind.AncB,function(POS){mean(POS[[DEME]])})}))

plot(flanking.blocks.deme.matrix.AncB[1,]/mean(flanking.blocks.deme.matrix.AncB[1,],na.rm=T),x=positions,ylim=c(0.5,3),type="l",lty=1,col="white",main="",cex.main=0.8,xlab="",ylab="",yaxt="n",xaxt="n")
for(i in 50:1){points(flanking.blocks.deme.matrix.AncB[i,]/mean(flanking.blocks.deme.matrix.AncB[i,],na.rm=T),x=positions,type="l",col=rainbow(50)[50-i])}
axis(1,at=seq(0.48,0.52,0.01),labels=seq(-0.02,0.02,0.01),cex.axis=0.9)
axis(2,at=seq(0.5,3,0.5),cex.axis=0.9,las=2)

axis(1,at=seq(0.48,0.52,0.01),labels=seq(-0.02,0.02,0.01),cex.axis=0.9)
axis(2,at=seq(0.5,3,0.5),cex.axis=0.9,las=2)
mtext("Distance from selected locus (M)",side=1,line=2.5)
mtext("Normalized mean block length",side=2,line=2.5)
legend("topright",title="distance from center:",legend = rev(c(seq(-25,-5,5),seq(5,25,5))),cex=0.7,col=rainbow(50)[seq(1,50,5)],lty=1,xjust=1)
dev.off()

pdf(height=4, width=6.25,file="ratioAdjacentBlocksAlongChromNoConditioning.pdf")
par(mar=c(3.5,3.5,0.5,0.5))
ratio.blocks = mean_per_deme/t(flanking.blocks.deme.matrix[50:1,])
plot(ratio.blocks[,1]/mean(ratio.blocks[,1],na.rm=T),x=positions,ylim=c(0.5,3),type="l",lty=1,col="white",main="",cex.main=0.8,xlab="",ylab="",yaxt="n",xaxt="n")
for(i in 1:50){points(ratio.blocks[,i]/mean(ratio.blocks[,i],na.rm=T),x=positions,type="l",col=rainbow(50)[i])}
axis(1,at=seq(0.48,0.52,0.01),labels=seq(-0.02,0.02,0.01),cex.axis=0.9)
axis(2,at=seq(0.5,3,0.5),cex.axis=0.9,las=2)

mtext("Distance from selected locus (M)",side=1,line=2.5)
mtext("Normalized mean block length",side=2,line=2.5)
legend("topright",title="distance from center:",legend = rev(c(seq(-25,-5,5),seq(5,25,5))),cex=0.7,col=rainbow(50)[seq(1,50,5)],lty=1,xjust=1)
dev.off()



pdf(height=4, width=6.25,file="ratioAdjacentBlocksAlongChromHeatmapNoConditioning.pdf")

layout(c(1,2),heights=c(1,4))
par(mar=c(0.5,3.5,2,0.5))

	BREAKS = seq(1.4,5.5,0.25)
	COLORS = heat.colors(length(BREAKS)-1)

	plot(0,0,col="white",xlim=range(BREAKS),ylim=c(0,1),yaxt="n",ylab="",xlab="",bty="n",xaxt="n")
	polygon(c(BREAKS[1],BREAKS[1],BREAKS[length(BREAKS)],BREAKS[length(BREAKS)]),y=c(0,1,1,0))
	for(i in 1:(length(BREAKS)-1)){
		polygon(c(BREAKS[i],BREAKS[i],BREAKS[i+1],BREAKS[i+1]),y=c(0,1,1,0),col=COLORS[i],border=COLORS[i])	
	}	

	axis(line=,side=3,(BREAKS[-1]+BREAKS[1:(length(BREAKS)-1)])/2,labels=(BREAKS[-1]+BREAKS[1:(length(BREAKS)-1)])/2,cex.axis=0.7,lwd=0,lwd.ticks=1)

	par(mar=c(3.5,3.5,0,0.5))
	image(ratio.blocks[,50:1],x=positions,y=-25:25,main="",yaxt="n",xaxt="n",ylab="",xlab="",breaks=BREAKS,col=heat.colors(length(BREAKS)-1))
	axis(1,at=seq(0.48,0.52,0.01),labels=seq(-0.02,0.02,0.01),cex.axis=0.9)
	axis(2,at=seq(-24.5,0,5),labels=seq(-25,-5,5),cex.axis=0.9,las=2)
	axis(2,at=seq(0,24.5,5),labels=seq(5,25,5),cex.axis=0.9,las=2)
	mtext("Distance from selected locus (M)",side=1,line=2.5)
	mtext("Distance from HZ center",side=2,line=2.5)
dev.off()




pdf(height=4, width=6.25,file="ratioAdjacentBlocksAlongChromAncBConditioning.pdf")
par(mar=c(3.5,3.5,0.5,0.5))
ratio.blocks.AncB = mean_per_deme_AncB/t(flanking.blocks.deme.matrix.AncB[50:1,])
plot(ratio.blocks.AncB[,1]/mean(ratio.blocks.AncB[,1],na.rm=T),x=positions,ylim=c(0.5,3),type="l",lty=1,col="white",main="",cex.main=0.8,xlab="",ylab="",yaxt="n",xaxt="n")
for(i in 1:50){points(ratio.blocks.AncB[,i]/mean(ratio.blocks.AncB[,i],na.rm=T),x=positions,type="l",col=rainbow(50)[i])}
axis(1,at=seq(0.48,0.52,0.01),labels=seq(-0.02,0.02,0.01),cex.axis=0.9)
axis(2,at=seq(0.5,3,0.5),cex.axis=0.9,las=2)

mtext("Distance from selected locus (M)",side=1,line=2.5)
mtext("Normalized mean block length",side=2,line=2.5)
legend("topright",title="distance from center:",legend = rev(c(seq(-25,-5,5),seq(5,25,5))),cex=0.7,col=rainbow(50)[seq(1,50,5)],lty=1,xjust=1)
dev.off()

pdf(height=4, width=6.25,file="ratioAdjacentBlocksAlongChromHeatmapAncBConditioning.pdf")

layout(c(1,2),heights=c(1,4))
par(mar=c(0.5,3.5,2,0.5))

	BREAKS = seq(0.25,8.5,0.25)
	COLORS = heat.colors(length(BREAKS)-1)

	plot(0,0,col="white",xlim=range(BREAKS),ylim=c(0,1),yaxt="n",ylab="",xlab="",bty="n",xaxt="n")
	polygon(c(BREAKS[1],BREAKS[1],BREAKS[length(BREAKS)],BREAKS[length(BREAKS)]),y=c(0,1,1,0))
	for(i in 1:(length(BREAKS)-1)){
		polygon(c(BREAKS[i],BREAKS[i],BREAKS[i+1],BREAKS[i+1]),y=c(0,1,1,0),col=COLORS[i],border=COLORS[i])	
	}	

	axis(line=,side=3,(BREAKS[-1]+BREAKS[1:(length(BREAKS)-1)])/2,labels=(BREAKS[-1]+BREAKS[1:(length(BREAKS)-1)])/2,cex.axis=0.7,lwd=0,lwd.ticks=1)

	par(mar=c(3.5,3.5,0,0.5))
	image(ratio.blocks.AncB[,50:1],x=positions,y=-25:25,main="",yaxt="n",xaxt="n",ylab="",xlab="",breaks=BREAKS,col=heat.colors(length(BREAKS)-1))
	axis(1,at=seq(0.48,0.52,0.01),labels=seq(-0.02,0.02,0.01),cex.axis=0.9)
	axis(2,at=seq(-24.5,0,5),labels=seq(-25,-5,5),cex.axis=0.9,las=2)
	axis(2,at=seq(0,24.5,5),labels=seq(5,25,5),cex.axis=0.9,las=2)
	mtext("Distance from selected locus (M)",side=1,line=2.5)
	mtext("Distance from HZ center",side=2,line=2.5)
dev.off()


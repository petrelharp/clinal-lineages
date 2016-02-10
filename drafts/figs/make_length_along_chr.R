image.scale <- function(z, zlim, col = heat.colors(12),
breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
 if(!missing(breaks)){
  if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
 }
 if(missing(breaks) & !missing(zlim)){
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
 }
 if(missing(breaks) & missing(zlim)){
  zlim <- range(z, na.rm=TRUE)
  zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
  zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
 }
 poly <- vector(mode="list", length(col))
 for(i in seq(poly)){
  poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
 }
 xaxt <- ifelse(horiz, "s", "n")
 yaxt <- ifelse(horiz, "n", "s")
 if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
 if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
 if(missing(xlim)) xlim=XLIM
 if(missing(ylim)) ylim=YLIM
 plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
 for(i in seq(poly)){
  if(horiz){
   polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
  }
  if(!horiz){
   polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
  }
 }
}




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


load("drafts/figs/simulation_SIGMA1_Ninds25000_ndemes50_s0_tau1000_chr1_start0.48_stop0.52_by0.001_intervalSizes_allAncs.Robj")
load("drafts/figs/simulation_SIGMA1_Ninds25000_ndemes50_s0_tau1000_chr1_start0.48_stop0.52_by0.001_intervalSizes.Robj")

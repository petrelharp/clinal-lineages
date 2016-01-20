#Want to produce a plot of the haplotypes
setwd("~/Documents/Hybrid_Zones/clinal-lineages/")

plot.chroms = function(DATA=sims.sums,INDS=random_inds,CHR="chr1"){
	
		plot(0:1,0:1,col="white",bty="n",xaxt="n",yaxt="n",xlab="",ylab="")

	for(i in 1:length(random_inds)){
		data=DATA[[1]][["ind.ancest"]][[INDS[i]]][[CHR]]$X1
		rect(data$starts,(i-1)/50+0.001,data$stops,i/50-0.001,col=ifelse(data$sp2,"orange","red"),border=FALSE)
		#segments(0,i/50,1,i/50,col="orange",lwd=5,lend="butt")
		#if(length(which(data$sp2))>0){
		#segments(data$starts[data$sp2],i/50,data$stops[data$sp2],i/50,col="red",lwd=5,lend="butt")	}
	} 	
}

ndemes=50
random_inds = sapply(1:ndemes,function(X){ind_indices=1:500+(X-1)*500;return(sample(ind_indices,size=1))})


#TAU=100

chroms_sel= get(load("sims/simulation_SIGMA1_Ninds25000_ndemes50_s0.1_tau100_run2016-01-11_simsums.Robj"))
chroms_neut = get(load("sims/simulation_SIGMA1_Ninds25000_ndemes50_s0_tau100_run2016-01-15_simsums.Robj"))
	#output is:  a list: [[zone_age]][[c("ind.ancest","ancest.prop")]]
		#"ind.ancestry" is a list of length ninds: [[ind][[chr]]][[X1/X2]][[data.frame where col1 = (chunk)starts, col2 = (chunk)stops, col3 = sp2(T/F)]]
		#"ancest.prop" is a ninds*2 matrix giving the ancestry proportion of each individual on each of their two copies of the chromosome.

pdf(file=c("sims/plot_chromosomes_tau100.pdf"))
par(mfrow=c(1,2),mar=c(0,0,0,0),oma=c(4,4,2,1))

#compare s=0 simulation to s=0.1 simulation
plot.chroms(DATA = chroms_neut,INDS=random_inds,CHR="chr1")
	#x-axis
	axis(side=1,at=c(0.01,seq(0.2,0.8,0.2),0.98),labels=seq(0,100,20),tick=F,cex.axis=0.75,line=-1.5)
	mtext("Position (cM)",side=1,cex=0.75,line=0.5)
	#y-axis
	axis(side=2,at=seq(1-4/ndemes,0,-1/ndemes*5),labels=seq(5,ndemes,5),las=1,tick=F,cex.axis=0.8,line=-1)
	mtext("deme",side=3,at=-0.08,cex=0.8,line=-0.5)
	mtext("s=0",line=0.5,cex=1.2,font=2)


plot.chroms(DATA=chroms_sel,INDS=random_inds,CHR="chr1")
	axis(side=1,at=c(0.01,seq(0.2,0.8,0.2),0.98),labels=seq(0,100,20),tick=F,cex.axis=0.75,line=-1.5)
	mtext("Position (cM)",side=1,cex=0.75,line=0.5)
	mtext("s=0.1",line=0.5,cex=1.2,font=2)
	arrows(0.5, -0.05, 0.5, 0, xpd = TRUE,length=0.075,lwd=1.5)

#compare selected chrom to non-selected chrom in s=0.1 simulation
plot.chroms(DATA = chroms_sel,INDS=random_inds,CHR="chr2")
	#x-axis
	axis(side=1,at=c(0.01,seq(0.2,0.8,0.2),0.98),labels=seq(0,100,20),tick=F,cex.axis=0.75,line=-1.5)
	mtext("Position (cM)",side=1,cex=0.75,line=0.5)
	#y-axis
	axis(side=2,at=seq(1-4/ndemes,0,-1/ndemes*5),labels=seq(5,ndemes,5),las=1,tick=F,cex.axis=0.8,line=-1)
	mtext("deme",side=3,at=-0.08,cex=0.8,line=-0.5)
	mtext("neutral chromosome",line=0.5,cex=1.2,font=2)
	

plot.chroms(DATA=chroms_sel,INDS=random_inds,CHR="chr1")
	axis(side=1,at=c(0.01,seq(0.2,0.8,0.2),0.98),labels=seq(0,100,20),tick=F,cex.axis=0.75,line=-1.5)
	mtext("Position (cM)",side=1,cex=0.75,line=0.5)
	mtext("s=0.1",line=0.5,cex=1.2,font=2)
	arrows(0.5, -0.05, 0.5, 0, xpd = TRUE,length=0.075,lwd=1.5)

dev.off()




chroms_sel = get(load("~/Documents/Hybrid_Zones/clinal-lineages/sims/simulation_SIGMA1_Ninds25000_ndemes50_s0.1_tau1000_simsums.Robj"))
chroms_neut = get(load("~/Documents/Hybrid_Zones/clinal-lineages/sims/simulation_SIGMA1_Ninds25000_ndemes50_s0_tau1000_simsums.Robj"))

pdf(file=c("sims/plot_chromosomes_tau1000.pdf"))
par(mfrow=c(1,2),mar=c(0,0,0,0),oma=c(4,4,2,1))

#compare s=0 simulation to s=0.1 simulation
plot.chroms(DATA = chroms_neut,INDS=random_inds,CHR="chr1")
	#x-axis
	axis(side=1,at=c(0.01,seq(0.2,0.8,0.2),0.98),labels=seq(0,100,20),tick=F,cex.axis=0.75,line=-1.5)
	mtext("Position (cM)",side=1,cex=0.75,line=0.5)
	#y-axis
	axis(side=2,at=seq(1-4/ndemes,0,-1/ndemes*5),labels=seq(5,ndemes,5),las=1,tick=F,cex.axis=0.8,line=-1)
	mtext("deme",side=3,at=-0.08,cex=0.8,line=-0.5)
	mtext("s=0",line=0.5,cex=1.2,font=2)
	

plot.chroms(DATA=chroms_sel,INDS=random_inds,CHR="chr1")
	axis(side=1,at=c(0.01,seq(0.2,0.8,0.2),0.98),labels=seq(0,100,20),tick=F,cex.axis=0.75,line=-1.5)
	mtext("Position (cM)",side=1,cex=0.75,line=0.5)
	mtext("s=0.1",line=0.5,cex=1.2,font=2)
	arrows(0.5, -0.05, 0.5, 0, xpd = TRUE,length=0.075,lwd=1.5)
	
	
	
#compare selected chrom to non-selected chrom in s=0.1 simulation
plot.chroms(DATA = chroms_sel,INDS=random_inds,CHR="chr2")
	#x-axis
	axis(side=1,at=c(0.01,seq(0.2,0.8,0.2),0.98),labels=seq(0,100,20),tick=F,cex.axis=0.75,line=-1.5)
	mtext("Position (cM)",side=1,cex=0.75,line=0.5)
	#y-axis
	axis(side=2,at=seq(1-4/ndemes,0,-1/ndemes*5),labels=seq(5,ndemes,5),las=1,tick=F,cex.axis=0.8,line=-1)
	mtext("deme",side=3,at=-0.08,cex=0.8,line=-0.5)
	mtext("neutral chromosome",line=0.5,cex=1.2,font=2)
	

plot.chroms(DATA=chroms_sel,INDS=random_inds,CHR="chr1")
	axis(side=1,at=c(0.01,seq(0.2,0.8,0.2),0.98),labels=seq(0,100,20),tick=F,cex.axis=0.75,line=-1.5)
	mtext("Position (cM)",side=1,cex=0.75,line=0.5)
	mtext("s=0.1",line=0.5,cex=1.2,font=2)

	arrows(0.5, -0.05, 0.5, 0, xpd = TRUE,length=0.075,lwd=1.5)
dev.off()
	

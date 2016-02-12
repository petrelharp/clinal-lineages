load("simulation_SIGMA1_Ninds25000_ndemes50_s0.01_tau1000_run2016-01-27_simsums.Robj")
chromosome = "chr1"
plot_start = 0.48
plot_stop=0.52
increments = 0.001
Ngen = 1000
sel=0
SIGMA=1
Ninds=25000
Ndemes = 50

positions = seq(plot_start,plot_stop,increments)

outstring = sprintf("simulation_SIGMA%d_Ninds%d_ndemes50_s%g_tau%d_%s_start%g_stop%g_by%g",SIGMA,Ninds,sel,Ngen,chromosome,plot_start,plot_stop,increments)

deme_ID = do.call(c,lapply(1:Ndemes,function(X){rep(X,Ninds/Ndemes)}))

#1. For each site, add up scores of all overlapping haplotypes
#2. Ways to score:
#	1. length
#	2. score = $\log$ of ECDF of length (i.e. log of prob that haplotype is longer than this one)

#lengths:
#lengths given a particular ancestry
get.interval.size = function(IND_DATA=sims.sums[[1]]$ind.ancest[[1]],CHR=chromosome,POS=0.5,ancB=TRUE){
	chunk = sapply(IND_DATA[[CHR]],function(X){diff(as.numeric(X[which(X$starts<POS & X$stops>POS),1:2]))})
	identity = sapply(IND_DATA[[CHR]],function(X){X[which(X$starts<POS & X$stops>POS),3]})
	replace(chunk,which(identity==ancB),0)
	
}

intervalSizes = lapply(positions,function(POS){
				do.call(rbind,lapply(sims.sums[[1]]$ind.ancest,function(IND){
				get.interval.size(IND_DATA=IND,CHR=chromosome,POS=POS,ancB=TRUE)	
				}))
				})
				
#lengths of any ancestry:
get.interval.size = function(IND_DATA=sims.sums[[1]]$ind.ancest[[1]],CHR=chromsoome,POS=0.5,ancB=TRUE){
	chunk = sapply(IND_DATA[[CHR]],function(X){diff(as.numeric(X[which(X$starts<POS & X$stops>POS),1:2]))})	
}

intervalSizes_allAncs = lapply(positions,function(POS){
				do.call(rbind,lapply(sims.sums[[1]]$ind.ancest,function(IND){
				get.interval.size(IND_DATA=IND,CHR=chromosome,POS=POS,ancB=TRUE)	
				}))
				})

#[[POS]][ind,chr]

#COMPARE:
pdf(file=paste(outstring,"statistic_length.pdf",sep="_"),height = 3,width=5)
#mean length across ALL indivduals in the zone (including non-recombinant chroms) vs chromosome position

plot(positions,sapply(intervalSizes_allAncs,mean),type="l",xlab="physical position",ylab="mean length",main="mean length, whole zone")

#mean length across ALL individuals of ancestryB at given locus:
points(positions,sapply(intervalSizes,function(X){mean(X[which(X>0)])}),type="l",col="red")

legend("topright",legend=c("ALL","AncB"),col=c("black","red"),lty=1)

#PLOT CURVE FOR EACH DEME:

#mean length across ALL individuals in deme (i.e. ancestry doesn't matter)
mean_per_deme = do.call(rbind, lapply(intervalSizes_allAncs,function(P){tapply(1:25000,deme_ID,function(Z){mean(P[Z,])})}))
plot(mean_per_deme[,1]/mean(mean_per_deme[,1],na.rm=T),ylim=c(0.5,2),type="l",lty=1,col=c(rainbow(25),rev(rainbow(25))),main="Any ancestry, normalized by deme (col reflects dist from center)",cex.main=0.8)
for(i in 1:50){points(mean_per_deme[,i]/mean(mean_per_deme[,i],na.rm=T),type="l",col=c(rainbow(25),rev(rainbow(25)))[i])}

#mean length across ALL individuals of ancestryB at given locus:
mean_per_deme_AncB = do.call(rbind, lapply(intervalSizes,function(P){tapply(1:25000,deme_ID,function(Z){all_sites = P[Z,]; mean(all_sites[which(all_sites>0)])})}))
plot(mean_per_deme_AncB[,1]/mean(mean_per_deme_AncB[,1],na.rm=T),ylim=c(0.5,2),type="l",lty=1,col=c(rainbow(25),rev(rainbow(25))),main="Any ancestry, normalized by deme (col reflects deme pos)",cex.main=0.8)
for(i in 1:50){points(mean_per_deme_AncB[,i]/mean(mean_per_deme_AncB[,i],na.rm=T),type="l",col=rainbow(50)[i])}


	image(mean_per_deme,x=positions,y=-25:25,main="all inds")
	image(mean_per_deme_AncB,x=positions,y=-25:25,main="AncB")


dev.off()

pdf(file=paste(outstring,"statistic_ecdf.pdf",sep="_"),height = 3,width=5)

#	2. score = $\log$ of ECDF of length (i.e. log of prob that haplotype is longer than this one)
#relative to ALL individuals in ALL populations:
everyone_ecdf = ecdf(do.call(rbind, intervalSizes_allAncs))
mean_Pr_hap_longer = sapply(intervalSizes_allAncs,function(POS){mean(everyone_ecdf(POS))})
plot(mean_Pr_hap_longer,type="l",ylim=c(0.2,0.6))

#relative to ANCESTRYB individuals in ALL populations:
AncB_everyone_ecdf = ecdf(do.call(c, lapply(intervalSizes,function(X){X[which(X>0)]})))
mean_Pr_hap_longer_AncB = sapply(intervalSizes,function(POS){mean(AncB_everyone_ecdf(POS))})
points(mean_Pr_hap_longer_AncB,type="l",col="red")


#Plot for each deme:
#Everyone:
	everyone_deme_ecdf = lapply(1:50,function(DEME){ecdf(do.call(c,lapply(intervalSizes_allAncs,function(POS){POS[which(deme_ID==DEME),]})))})
	
	matrix_of_mean_probs_anyone_by_deme = do.call(rbind,lapply(1:Ndemes,function(DEME){
		ECDF = everyone_deme_ecdf[[DEME]];sapply(1:length(positions),function(POS){
			mean(ECDF(intervalSizes_allAncs[[POS]][which(deme_ID==DEME),]))
			})}))
	
	matplot(x=positions,log(t(matrix_of_mean_probs_anyone_by_deme)),type="l",lty=1,col=rainbow(50),main="by deme, all inds")	
		
#condition on being AncB:
	#Get rid of 0 ancestry guys in AncB:
	intervalSizes_AncBinds = lapply(intervalSizes,function(POS){lapply(1:50,function(DEME){bydeme=POS[which(deme_ID==DEME),];bydeme[which(bydeme>0)]})})
	
	everyone_deme_ecdf_AncB = lapply(1:50,function(DEME){ecdf(do.call(c,lapply(intervalSizes_AncBinds,function(POS){POS[[DEME]]})))})

	matrix_of_mean_probs_AncB_by_deme = do.call(rbind,lapply(1:50,function(DEME){
		ECDF = everyone_deme_ecdf_AncB[[DEME]];sapply(1:length(positions),function(POS){
			mean(ECDF(intervalSizes_AncBinds[[POS]][[DEME]]))
			})}))

	matplot(x=positions,log(t(matrix_of_mean_probs_AncB_by_deme)),type="l",lty=1,col=rainbow(50),main="by deme | AncB")	

	image(log(t(matrix_of_mean_probs_anyone_by_deme)),x=positions,y=-25:25,main="all_inds")
	image(log(t(matrix_of_mean_probs_AncB_by_deme)),x=positions,y=-25:25,main="AncB")

dev.off()

#5. Look also at length of *next* block over?

#For ALL inds:

pdf(file=paste(outstring,"statistic_length_of_adjacent_blocks.pdf",sep="_"),height = 3,width=5)
get.flanking.blocks.all = function(IND_DATA=sims.sums[[1]]$ind.ancest[[1]],CHR=chromosome,POS=0.5,ancB=TRUE){
	focal_chunks = do.call(rbind,lapply(IND_DATA[[CHR]],function(X){
		FOCUS = which(X$starts<POS & X$stops>POS);
		if(max(FOCUS-1,1)!=min(FOCUS+1,nrow(X))){return(X[c(max(FOCUS-1,1),min(FOCUS+1,nrow(X))),])}else return(X[max(FOCUS-1,1),])}
		))
	focal_chunks$stops-focal_chunks$starts
	}
	
flanking.blocks.by.ind = lapply(positions,function(POS){
	tapply(1:25000,deme_ID,function(DEME){do.call(c,lapply(DEME,function(IND){
			get.flanking.blocks.all(IND_DATA=sims.sums[[1]]$ind.ancest[[IND]],POS=POS)}))})})
	
plot(positions,sapply(flanking.blocks.by.ind,function(POS){mean(do.call(c,POS))}),type="l",xlab="physical position",ylab="mean length")

#By DEME:

flanking.blocks.deme.matrix = do.call(rbind,lapply(1:50,function(DEME){sapply(flanking.blocks.by.ind,function(POS){mean(POS[[DEME]])})}))

matplot(t(flanking.blocks.deme.matrix/apply(flanking.blocks.deme.matrix,1,mean)),x=positions,type="l",col=rainbow(50),lty=1,main="mean length, all inds in deme",cex.main=0.7)	


#For AncB inds:
get.flanking.blocks.AncB = function(IND_DATA=sims.sums[[1]]$ind.ancest[[1]],CHR=chromosome,POS=0.5,ancB=TRUE){
	focal_chunks = do.call(rbind,lapply(IND_DATA[[CHR]],function(X){
		FOCUS = which(X$starts<=POS & X$stops>=POS);
		if(X[FOCUS,]$sp2){
		if(max(FOCUS-1,1)!=min(FOCUS+1,nrow(X))){return(X[c(max(FOCUS-1,1),min(FOCUS+1,nrow(X))),])}else return(X[max(FOCUS-1,1),])}else return()}
		))
	focal_chunks$stops-focal_chunks$starts
	}

flanking.blocks.by.ind.AncB = lapply(positions,function(POS){
	tapply(1:25000,deme_ID,function(DEME){do.call(c,lapply(DEME,function(IND){
			get.flanking.blocks.AncB(IND_DATA=sims.sums[[1]]$ind.ancest[[IND]],POS=POS)}))})})

plot(positions,sapply(flanking.blocks.by.ind.AncB,function(POS){mean(do.call(c,POS))}),type="l",xlab="physical position",ylab="mean length",main="mean length | AncB",cex.main=0.7)

flanking.blocks.deme.matrix.AncB = do.call(rbind,lapply(1:50,function(DEME){sapply(flanking.blocks.by.ind.AncB,function(POS){mean(POS[[DEME]])})}))

#matplot(t(flanking.blocks.deme.matrix.AncB/apply(flanking.blocks.deme.matrix.AncB,1,mean)),x=positions,type="l",col=rainbow(50),lty=1)	
plot(xlim=c(0,20),ylim=c(0,2),0,0,col="white",ylab="mean length",xlab="physical position",main="flanking block length| AncB, normalized deme",cex.main=0.7)
for(i in 1:50){points(flanking.blocks.deme.matrix.AncB[i,]/mean(flanking.blocks.deme.matrix.AncB[i,],na.rm=T),type="l",col=rainbow(50)[i])}



plot(xlim=c(0,20),ylim=c(0,2),0,0,col="white",ylab="mean length",xlab="physical position",main="flanking block length| AncB, normalized deme",cex.main=0.7)
for(i in 1:50){points(flanking.blocks.deme.matrix.AncB[i,]/mean(flanking.blocks.deme.matrix.AncB[i,],na.rm=T),type="l",col=rainbow(50)[i])}

matplot((mean_per_deme)/t(flanking.blocks.deme.matrix),type="l",col=rainbow(50),lty=1,main="mean_focal/mean_flanking length ALL" ,cex.main=0.7)

matplot((mean_per_deme_AncB)/t(flanking.blocks.deme.matrix.AncB[50:1,]),type="l",col=rainbow(50),lty=1,main="mean_focal/mean_flanking length AncB" ,cex.main=0.7)

image(t(flanking.blocks.deme.matrix),x=positions,y=-25:25,ylab="deme",xlab="phys.pos",main="heatmap ALL (log)")
image(t(flanking.blocks.deme.matrix.AncB),x=positions,y=-25:25,ylab="deme",xlab="phys.pos",main="heatmap AncB (log)")



image(((mean_per_deme)/t(flanking.blocks.deme.matrix[50:1,])),x=positions,y=-25:25,ylab="deme",xlab="phys.pos",main="focal/flanking heatmap ALL (log)")
image(((mean_per_deme_AncB)/t(flanking.blocks.deme.matrix.AncB[50:1,])),x=positions,y=-25:25,ylab="deme",xlab="phys.pos",main="focal/flanking heatmap AncB (log)")

dev.off()
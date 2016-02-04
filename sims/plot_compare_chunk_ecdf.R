load("~/Documents/Hybrid_Zones/clinal-lineages/sims/simulation_SIGMA1_Ninds25000_ndemes50_s0.1_tau1000_simsums_chunks.Robj")
#this is a list of [selected,near,unlinked chunk distribution].The object name is "chunks"

load("~/Documents/Hybrid_Zones/clinal-lineages/sims/simulation_SIGMA1_Ninds25000_ndemes50_s0_tau1000_chunks_by_deme.Robj")
#this is only around the selected locus. The object name is "chunks_by_deme"


transparent_rainbow = adjustcolor(rainbow(ndemes),0.7)


empty_deme = rep(0,2*deme_size)

chunks_ecdf = lapply(chunks,function(X){
	lapply(X,function(D){
		segregating =which(D<1 & D>0)
		if(length(segregating)>0){
		return(ecdf(D[segregating]))}else return(ecdf(D))	
	})})


noSelection_ecdf = lapply(chunks_by_deme,function(D){segregating =which(D<1 & D>0)
		if(length(segregating)>0){
		return(ecdf(D[segregating]))}else return(ecdf(D))})

#vanilla ecdfs	
pdf(file="simulation_SIGMA1_Ninds25000_ndemes50_s0_tau1000_ecdf_comparisons.pdf")	
par(mfrow=c(2,2),mar=c(2,2,2,2))
for(type in names(chunks)){
	plot(1, type="n",ylim=c(0,1),xlim=c(0,0.2),main=type,ylab="1-ecdf(x)",xlab="x=chunk length")	
	#legend("topright",legend=seq(-24.5,24.5,5),col=transparent_rainbow[seq(1,50,5)],lty=1)
	for(i in rev(1:ndemes)){points(1-chunks_ecdf[[type]][[i]](seq(0,1,0.005))~seq(0,1,0.005),ylim=c(0,1),type="l",col=transparent_rainbow[i])}
}

	plot(1, type="n",ylim=c(0,1),xlim=c(0,0.2),main=type,ylab="1-ecdf(x)",xlab="x=chunk length")	
	legend("topright",legend=seq(-24.5,24.5,5),col=transparent_rainbow[seq(1,50,5)],lty=1)
	for(i in rev(1:ndemes)){points(1-noSelection_ecdf[[i]](seq(0,1,0.005))~seq(0,1,0.005),ylim=c(0,1),type="l",col=transparent_rainbow[i])}

dev.off()
#ecdf diffs:
pdf(file="simulation_SIGMA1_Ninds25000_ndemes50_s0_tau1000_ecdf_difference_from_neutral.pdf")	

par(mfrow=c(2,2),mar=c(2,2,2,2))
for(type in names(chunks)){
	plot(1, type="n",ylim=c(-1,1),xlim=c(0,0.2),main=type,ylab="1-ecdf(x)",xlab="x=chunk length")	
	#legend("topright",legend=seq(-24.5,24.5,5),col=transparent_rainbow[seq(1,50,5)],lty=1)
	for(i in rev(1:ndemes)){points((chunks_ecdf[[type]][[i]](seq(0,1,0.005))-noSelection_ecdf[[i]](seq(0,1,0.005)))~seq(0,1,0.005),ylim=c(0,1),type="l",col=transparent_rainbow[i])}
}

	plot(1, type="n",ylim=c(0,1),xlim=c(0,0.2),main="deme position legend",ylab="1-ecdf(x)",xlab="x=chunk length")	
	legend("topright",legend=seq(-24.5,24.5,5),col=transparent_rainbow[seq(1,50,5)],lty=1,bty="n")

dev.off()

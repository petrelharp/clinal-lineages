start_coordinate=0.49
end_coordinate=0.51

load("simulation_SIGMA1_Ninds25000_ndemes50_s0.01_tau1000_run2016-01-27_rawoutput.Robj")	#output is: a list with hierarchy: [[zone_age]][[c("inds","sp.ids","deme","pars","QTL")]]
		#where:
			#"inds" is a list of length ninds, such that [[ind]][[chr]][[X1/X2]][[data.frame where col1 = coordinates of rec., col2 = ancestor ID]]
			#"sp.ids" is a character vector of length ninds, giving the ancestry given an ancestor ID
			#"deme" is an integer vector of length ninds, giving the deme number of an ancestor, given their ancestor ID
			#"pars" is an integer that gives the number of generations of the simulation
			#"QTL" returns the QTL vector used in the simulation


#Let's get the number of ancestors represented in chunks in each deme:

get.interval.size = function(IND_DATA=sims.sums[[1]]$ind.ancest[[1]],CHR=1,POS=0.5,ancA = TRUE){
	#return the interval containing focal site for a individual	
	chunk = sapply(IND_DATA[[CHR]],function(X){diff(as.numeric(X[which(X$starts<POS & X$stops>POS),1:2]))})	
	identity = sapply(IND_DATA[[CHR]],function(X){X[which(X$starts<POS & X$stops>POS),3]})
	replace(chunk,which(identity==ancA),0)
}


get.deme.chunks = function(IND_DATA=sims.sums[[1]]$ind.ancest, DEME = 1, CHR=1,POS=0.5,ancA = TRUE){
	#get distribution of chunks within a deme
	INDS = IND_DATA[which(sims.s0.1[[1]]$deme==DEME)]
	intervals = do.call(rbind,lapply(INDS,get.interval.size,CHR=CHR,POS=POS,ancA=ancA)) 	
	return(intervals)
}




#returns table (ancestor at locus, and whether it is AncB (T/F))
ancestry.at.position = function(IND=1,DATA=sims.s0.1,POS = 0.5,CHROM="chr1",SP="B"){
	distance_from_focal_pos = lapply(DATA[[1]][["inds"]][[IND]][[CHROM]],function(X){d = POS-X$pos;return(X[max(which(d>0)),])})
	sp_id = lapply(distance_from_focal_pos,function(X){data.frame(ancest=X$ancest,spB=DATA[[1]][["sp.ids"]][ceiling(X$ancest/2)]==SP)})
}

ninds=25000

selected_chromsome_ancestor_bydeme = lapply(seq(start_coordinate,end_coordinate,0.00025),function(P){tapply(1:ninds,sims.s0.1[[1]][["deme"]],function(X){lapply(X,ancestry.at.position,DATA=sims.s0.1,POS=P,CHROM="chr1",SP="B")})})
unlinked_chromsome_ancestor_bydeme = lapply(seq(start_coordinate,end_coordinate,0.00025),function(P){tapply(1:ninds,sims.s0.1[[1]][["deme"]],function(X){lapply(X,ancestry.at.position,DATA=sims.s0.1,POS=P,CHROM="chr2",SP="B")})})

selected_locus_AncB_ancestors = lapply(selected_chromsome_ancestor_bydeme,function(POS){
	lapply(POS,function(X){table_of_inds=do.call(rbind,lapply(X,function(Z){
			do.call(rbind,Z)})); return(table_of_inds$ancest[which(table_of_inds$spB)])})})
unlinked_locus_AncB_ancestors = lapply(unlinked_chromsome_ancestor_bydeme,function(POS){lapply(POS,function(X){
	table_of_inds=do.call(rbind,lapply(X,function(Z){
		do.call(rbind,Z)})); return(table_of_inds$ancest[which(table_of_inds$spB)])})})



load("simulation_SIGMA1_Ninds25000_ndemes50_s0_tau1000_rawoutput.Robj")
neutral_chromsome_ancestor_bydeme = lapply(seq(start_coordinate,end_coordinate,0.00025),function(P){tapply(1:ninds,sims.s0.1[[1]][["deme"]],function(X){lapply(X,ancestry.at.position,DATA=sims.s0.1,POS=P,CHROM="chr1",SP="B")})})

neutral_locus_AncB_ancestors = lapply(neutral_chromsome_ancestor_bydeme,function(POS){
	lapply(POS,function(X){table_of_inds=do.call(rbind,lapply(X,function(Z){
			do.call(rbind,Z)})); return(table_of_inds$ancest[which(table_of_inds$spB)])})})



number_of_ancestors = list(neutral_locus_AncB_ancestors,unlinked_locus_AncB_ancestors,selected_locus_AncB_ancestors)

save(number_of_ancestors,file="number_of_ancestors_tau1000_revision.Robj")
















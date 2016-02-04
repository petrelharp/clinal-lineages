*_rawoutput.Robj
		#output is: a list with hierarchy: [[zone_age]][[c("inds","sp.inds","deme","pars","QTL")]]
		#where:
			#"inds" is a list of length ninds, such that [[ind]][[chr]][[X1/X2]][[data.frame where col1 = coordinates of rec., col2 = ancestor ID]]
			#"sp.inds" is a character vector of length ninds, giving the ancestry given an ancestor ID
			#"deme" is an integer vector of length ninds, giving the deme number of an ancestor, given their ancestor ID
			#"pars" is an integer that gives the number of generatsion of the simulation
			#"QTL" returns the QTL vector used in the simulation

*_simsums.Robj
		#output is:  a list: [[zone_age]][[c("ind.ancestry","ancest.prop")]]
			#"ind.ancestry" is a list of length ninds: [[ind][[chr]]][[X1/X2]][[data.frame where col1 = (chunk)starts, col2 = (chunk)stops, col3 = sp2(T/F)]]
			#"ancest.prop" is a ninds*2 matrix giving the ancestry proportion of each individual on each of their two copies of the chromosome.

*_simsums_chunks.Robj
		chunks = list(selected=testB,distant=testB_far,unlinked=testB_unlinked)
		[[locus]][[deme]][[matrix of chunk sizes]]
		where the matrix is deme_size by ploidy. Each matrix entry gives the chunk size of Ancestry B surrounding locus (0 means the locus is of Ancestry A)


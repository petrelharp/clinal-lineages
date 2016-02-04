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

*_simsums_chunks.Robj : has a single object, "chunks", 
    - which is a list with components:
        - selected : at the selected site
        - distant : at recombination distance 0.01
        - unlinked : at the center of the other chromosome
    - and each of these is a list with one component per deme, each of which is
		- a matrix of chunk sizes, of dimensions deme_size by ploidy, with columns named e.g. "X1" and "X2" for the chromosomes.
        - Each matrix entry gives the chunk size of Ancestry B surrounding the locus (0 means the locus is of Ancestry A)


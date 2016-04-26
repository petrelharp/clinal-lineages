##
# Functions for extracting information from a _simsums_chunks object.
##

genotype_counts <- function (ch) { 
    gcounts <- t(sapply( ch, function (x) { tabulate( 1+rowSums( x>0 ), nbins=3 ) } )) 
    colnames(gcounts) <- c("AA","AB","BB")
    return(gcounts)
}

conditional_genotype_counts <- function (x,y) { 
    # returns joint genotype counts across space from two chunks objects
    #  a (ndemes) x 4 x 4 array (4 phased genotypes)
    deme <- rep(seq_along(x),sapply(x,NROW))
    x.alleles <- ifelse(do.call(rbind,x)==0,"A","B")
    y.alleles <- ifelse(do.call(rbind,y)==0,"A","B")
    x.genotypes <- factor( paste(x.alleles[,1],x.alleles[,2],sep=""), levels=c("AA","AB","BA","BB") )
    y.genotypes <- factor( paste(y.alleles[,1],y.alleles[,2],sep=""), levels=c("AA","AB","BA","BB") )
    return( as.array( table( deme, x=x.genotypes, y=y.genotypes ) ) )
}

dip.to.hap <- rbind( A=c(AA=2,AB=1,BB=0), B=c(AA=0,AB=1,BB=2) )

bidip.to.bihap <- rbind(
            "A-A" = c( 2, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0 ),
            "B-A" = c( 0, 1, 1, 2, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0 ),
            "A-B" = c( 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 2, 1, 1, 0 ),
            "B-B" = c( 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 2 ) )
colnames( bidip.to.bihap ) <- as.vector( outer( c("AA","AB","BA","BB"), c("AA","AB","BA","BB"), paste, sep="-" ) )

conditional_haplotype_counts <- function (gcounts) {
    # uses output of conditional_genotype_counts to get conditional haplotype counts
    # a (ndemes) x 2 x 2 array
    dim(gcounts) <- c(dim(gcounts)[1],prod(dim(gcounts)[2:3]))
    out <- tcrossprod( gcounts, bidip.to.bihap )
    dim(out) <- c(dim(gcounts)[1],2,2)
    dimnames(out) <- list( dimnames(gcounts)[[1]], c("A","B"), c("A","B") )
    return(out)
}

conditional_freqs <- function (gcounts) {
    # uses output of conditional_genotype_counts to get conditional allele frequencies
    # of A alleles at the FIRST locus, given the identity of the allele at the SECOND locus
    # (so, selected allele should be SECOND argument to conditional_genotype_counts)
    # a (ndemes) x 2 matrix
    hcounts <- conditional_haplotype_counts(gcounts)
    cbind( A = hcounts[,"A","A"]/(hcounts[,"A","A"]+hcounts[,"B","A"]),
           B = hcounts[,"A","B"]/(hcounts[,"A","B"]+hcounts[,"B","B"]) )
}


# do testing

# a chunks object is a list, corresponding to loci, each element is 
#   a list demes; each deme is a (nind) x (ploidy) matrix of chunk lengths of type B

test.chunk <- list(
        "r0.0" = list(
                "x1" = matrix( c(
                        1, 2,  # BB
                        0, 0,  # AA
                        0, 0), # AA
                    ncol=2, byrow=TRUE ),
                "x2" = matrix( c(
                        1, 0,  # BA
                        1, 3,  # BB
                        0, 0), # AA
                    ncol=2, byrow=TRUE ) ),
        "r0.1" = list(
                "x1" = matrix( c(
                        0, 2,  # AB
                        1, 1,  # BB
                        0, 0), # AA
                    ncol=2, byrow=TRUE ),
                "x2" = matrix( c(
                        1, 2,  # BB
                        1, 3,  # BB
                        0, 1), # AB
                    ncol=2, byrow=TRUE ) )
        )

# test genotype_counts
test.gcounts <- lapply( test.chunk, genotype_counts )
stopifnot( all.equal( test.gcounts,
        list(
            "r0.0" = rbind( # AA AB BB
                        x1 = c(2, 0, 1),
                        x2 = c(1, 1, 1) ),
            "r0.1" = rbind( # AA AB BB
                        x1 = c(1, 1, 1),
                        x2 = c(0, 1, 2) ) ),
        check.attributes=FALSE ) )

# test conditional_genotype_counts
#  here "AB-BB" has haplotype AB on one chromosome and BB on the other
test.cgcounts <- conditional_genotype_counts( test.chunk[[1]], test.chunk[[2]] )
true.genotypes <- list( x1=c( "BB-AB", "AA-BB", "AA-AA" ),
                         x2=c( "BA-BB", "BB-BB", "AA-AB" ) )
true.genocounts <- do.call( table, c( list( rep(c("x1","x2"),each=3), 
                            factor(sapply(strsplit(unlist(true.genotypes),"-"),"[",1), levels=c("AA","AB","BA","BB")),
                            factor(sapply(strsplit(unlist(true.genotypes),"-"),"[",2), levels=c("AA","AB","BA","BB")) ) ) )
stopifnot( all.equal( test.cgcounts, true.genocounts, check.attributes=FALSE ) )

# test conditional_haplotype_counts
test.chcounts <- conditional_haplotype_counts( test.cgcounts )
true.haplotypes <- list( x1=c("B-A","B-B","A-B","A-B","A-A","A-A"),
                         x2=c("B-B","A-B","B-B","B-B","A-A","A-B") )
true.haplocounts <- do.call( table, c( list( rep(c("x1","x2"),each=6), 
                            factor(sapply(strsplit(unlist(true.haplotypes),"-"),"[",1), levels=c("A","B")),
                            factor(sapply(strsplit(unlist(true.haplotypes),"-"),"[",2), levels=c("A","B")) ) ) )
stopifnot( all.equal( as.vector(test.chcounts), as.vector(true.haplocounts), check.attributes=FALSE ) )

# test conditional_freqs
test.cfreqs <- conditional_freqs( test.cgcounts )
true.cfreqs <- rbind( x1=c( A=2/3, B=2/3 ),
                      x2=c( A=1/1, B=2/5 ) )
stopifnot( all.equal( test.cfreqs, true.cfreqs, check.attributes=FALSE ) )

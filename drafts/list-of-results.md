
Main points
    1. tract lengths given ancestry tend to be longer/shorter conditional on linked ancestry
    2. selected chunks on the wrong side should be long (but rare)

Figures
    1. 
    2. distribution and difference of haplotype lengths
        a. conditioned on ancestry
        b. not conditioned on ancestry?
        c. for big tau?

## Outline

1.  Theory
    * lineages moving back in cline and jumping between backgrounds
        - illustrative picture
    * how fast should cline flatten 
    * how big should haplotypes be
    * talk about coalescent expectations?

2.  Simulations
    * methods
    * look at weak selection also?
    * can selection be strong enough that we see it on the other chromosome?  (at what time?)

9.  Appendix: Differential equations

## Ideas for statistics

**Sum of haplotype scores:**

1. For each site, add up scores of all overlapping haplotypes

2. Ways to score:
    1. length
    2. score = $\log$ of ECDF of length (i.e. log of prob that haplotype is longer than this one)

3. How to incorporate information from nearby locations?
    1. use fact that nearby haps of the other one should be shorter?
    2. do the above plot *conditioned* on ancestry *and* spatial location, and look for a dip with a peak in the middle

4. How to incorporate information from space?  No need?
    1. to check, compute at different spatial locations, should be the same

5. Look also at length of *next* block over?


**Density of transitions:**

Plot the density of transitions along the chromosome.

**Difference in chunk lengths:**

Look at the difference in chunk lengths between neighboring chunks.

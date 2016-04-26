These deal with environmental variation in fitness
* Haldane1948 : does 1D stable cline including dominance for step fitness gradient
* Fisher1950 : continuous gradient in selection
* hanson1966effects : does numerical solns to discretized ODE including observing minimal patch size
* endler1973population : does numerical solns for a few different varieties of externally specified fitness gradients,
    concluding that "there are many possible spatial patterns of selection and gene flow that can produce a given cline structure"
* Nagylaki1975 : writes down equations in 2D and conditions for stable solns
* Slatkin1975 and Maruyama : does numerical iteration of the discrete equations for a 1D cline, showing the expected cline is flatted by drift (but see Polechova2011)
* Felsenstein1975, genetic drift in clines : also looks at how genetic drift affects slope and center of clines but with different assumptions than Slatkin1975

These deal with hybrid zones
* Bazykin1969 : writes down the stable solution for the one-dimensional cline
* barton1979dynamics of hybrid zones: looks at how hybrid zones move using a potential; also looks at the strength of stochastic effects. seems to NOT say how fast they relax?
* barton1979gene flow past a cline : step in frequency of linked neutral allele 
     proportional to gradient in frequency and inverse to recomb rate
     adaptive introgression occurs quickly 
* barton1985 and hewitt, analysis of hybrid zones : a REVIEW of everything up 'til 85 (and some not yet published)
* barton1986 and bengtsson, the barrier to genetic exchange : looks at many linked incompatibilities
     with barrier B, neutral allele delayed by (B/sigma)^2, 
     while adaptively introgressing allele delayed only if (B/sigma)^2 > 1/s.
     also does size of step of allele frequency
* barton1986 effects of linkage : gets strength of barrier to linked neutral alleles, produces a "stepped cline"
* barton1989 and hewitt, adaptation, speciation, and hybrid zones : another review, shorter.
* barton2000 and shpak, effect of epistasis : more about LD on many linked selected loci
* barton2008effect of a barrier to gene flow : on fluctuations in allele frequency, using diffusion approximation; 
    finds is more sensitive in 2d than 1d

This deals with both:

* Slatkin1973 : does environmental cline and minimal width of a patch;
    - shows discontinuous freqs at barriers to dispersal
    - AND observes that het disadvantage leads to steeper clines,
    - so maybe environment determines location of cline; het disadvantage the shape

* may1975gene (with Endler, McMutrie) : spatial resolution of adaptation; also does het disadvantage;
    - reviews Slatkin and Endler's recent papers;
    - writes down time-evolution and does space-time scaling argument for environmental cline;
    - says that behavior of het disadvantage and environmental clines have "exactly the same scaling" properties, 
    - and "the clines differ only in numerical details".

* Barton1993 and gale, genetic analysis of hybrid zones: compares cline shapes of different mechanisms:
    - het disadvantage; an ecotone; an ecotone with dominance; stabilizing selection on a quantitative trait
    - also writes down PDE for stable pattern of linkage

* Kruuk1999 and Baird, Gale, and Barton : "a comparison of multilocus clines maintained by environmental clines or selection against hybrids"
    - "the shape of clines [produced by the two are] indistinguishable"
    - shows that pairwise LD generated primarily by migration
    - analytical prediction of environmental clines consistently too narrow
    - "Barton’s (1983, 1986) multilocus models invoked endogenous selection,
        specifically heterozygote disadvantage; simulations showed that the
        effects of epistatic endogenous selection were qualitatively similar (Barton
        and Gale 1993; Barton and Shpak 1999).  However, there is also substantial
        evidence of an important role for exogenous selection in hybrid zones (Arnold
        1997; Harrison 1990). In single-locus models, the shape of the clines generated
        by exogenous vs. endogenous viability selection is effectively indistinguish-
        able (Barton and Hewitt 1989; Barton and Gale 1993)"

* Polechova2011 and Barton, genetic drift widens the expected cline but narrows individual clines; 
     somewhat different answers for environmental step-clines and for tension zones

Other:

* cox1995 and durrett, "Hybrid zones and voter model interfaces"
    - show that the interface spreads out by diffusion in the neutral case

* durrett2007 and Zahle, "On the width of hybrid zones"
    - study the cline in an environmental step cline, showing (in low population density) it's
    - of width 1/s in d=1
    - of width sqrt( log(1/s) / s ) in d=2
    - of width 1/sqrt(s) in d>2

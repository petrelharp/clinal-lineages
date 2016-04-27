# Theory and numerical simulation

This directory contains Rmarkdown documents describing development of the theory,
implementation of numerical solutions,
and testing.

These contain functions used in the documents:

- [cline-fns.R](cline-fns.R) : functions to do most of the work of setting up and solving the relevant PDE; to see how they are used, read inline comments and the Rmarkdown documents.

- [tran.1D.R](tran.1D.R) : a slightly modified version of `ReacTran::tran.1D` that allows specifying fluxes on a log scale (thus avoiding numerical issues);
    
    * this is described in [reactran-notes.Rmd](reactran-notes.html)

- [unifprod.R](unifprod.R) : this is an implementation of the discretization of equation (5) in the text;
    or equivalently, finds multilocus probabilities of ancestral identity (so is used to compute LD, for instance); 
    see Appendix B for details.   It is used in the function `recomb_generator()` in `cline_fns.R`.


And, these are documents written while developing the code:

- [pde-clines.Rmd](pde-clines.html) : a good walk-through of what PDE we need to solve, what functions in `cline-fns.R` are used to solve them,
    and how those functions are structured.  
    Includes computation of some things we didn't get to in the paper, e.g., LD.

- [time-scales.Rmd](time-scales.html) : The same set of descriptive plots at a few time points (tens, hundreds, and thousands of generations).

- [rough-calculations.Rmd](rough-calculations.html) : Theoretical order-of-magnitude arguments, that are mostly in the paper.

- [diffusion-local-times.Rmd](diffusion-local-times.html) : Probably ignore this one. 
    Talks through some ideas for approximations to the actual dynamics of lineages,
    in the hopes of getting some better analytical approximations.
    We ended up just solving the correct equations; this document talks about how to set these up in a few places, also.



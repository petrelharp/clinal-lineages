# Beyond clines: lineages and haplotype blocks in hybrid zones

by [Alisa Sedghifar](https://github.com/asedghifar), [Yaniv Brandvain](https://brandvainlab.wordpress.com/), and [Peter Ralph](https://github.com/petrelharp)

To appear in Molecular Ecology; here's the preprint on the [bioRxiv](http://biorxiv.org/content/early/2016/03/12/043190.abstract).


This is a "working" repository.  It's not a software package.  
Our simulation code is not very efficient,
but the code for numerical solutions to the PDE is substantially more tidy and usable.
There's a number of dangling ends that we never finished off.
We're providing it so that others can see exactly what we did 
and so that code and methods can be useful for others
(but if you use our code, please acknowledge and/or cite us!).

Here's the important parts:

- [sims/](tree/master/sims/) : code to run, and plot results of, individual-based simulations.

    - [sims/README.md](tree/master/sims/README.md) : describes how to run a simulation and make some plots of the result

- [notes/](tree/master/notes/) : implementation of numerical solutions to the PDE, and documents describing development of the theory

    - [notes/README.md](tree/master/notes/README.md) : describes the relevant files in this directory

- [drafts/](tree/master/drafts/) : the writeup, in LaTeX. 

    - [drafts/hybrid_zone.tex](tree/master/drafts/hybrid_zone.tex) : the LaTeX for the document
    - [drafts/figs/](tree/master/drafts/figs/) : figures for the paper, and the code to make them.
    - [drafts/figs/README.md](tree/master/drafts/figs/README.md) : describes some of the R scripts that produce the figures
        *If you want to see how a given figure was produced, and it doesn't say here,
        try looking in the LaTeX source to find the filename, 
        then search for that string in the scripts in `figs/`.*

- [resources/](tree/master/resources/) : scripts to help out with making html reports from Rmarkdown, using [templater](https://github.com/petrelharp/templater).
    You can ignore this directory.


We also have produced a fairly large number of pdfs of simulation results under various parameter combinations.
These are found in `sims/`, in directories like `simulation_SIGMA3_Ninds6000_ndemes120_s0.01_dir/`.

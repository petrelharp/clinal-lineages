# To run a simulation 

Edit `run_sims.R` and execute it.
The parameters are at the top of the file:
```
params = list(
    SIGMA = 3,        # sigma for dispersal, in deme spacings
    ndemes = 100,     # number of demes
    deme_size = 50,   # number of diploid individuals per deme
    S = 0.05,         # selection coefficient, s
    zone_age = 10     # number of generations to run simulation for
)
```

It will create files in the format `simulation_SIGMA%s_Ninds%s_ndemes%s_s%s_dir/tau%s_dir/results_`,
including some useful plots.


# Description of resulting .Robj files

* `_rawoutput.Robj`
		#output is: a list with hierarchy: [[zone_age]][[c("inds","sp.inds","deme","pars","QTL")]]
		#where:
			#"inds" is a list of length ninds, such that [[ind]][[chr]][[X1/X2]][[data.frame where col1 = coordinates of rec., col2 = ancestor ID]]
			#"sp.inds" is a character vector of length ninds, giving the ancestry given an ancestor ID
			#"deme" is an integer vector of length ninds, giving the deme number of an ancestor, given their ancestor ID
			#"pars" is an integer that gives the number of generatsion of the simulation
			#"QTL" returns the QTL vector used in the simulation

* `_simsums.Robj`
		#output is:  a list: [[zone_age]][[c("ind.ancestry","ancest.prop")]]
			#"ind.ancestry" is a list of length ninds: [[ind][[chr]]][[X1/X2]][[data.frame where col1 = (chunk)starts, col2 = (chunk)stops, col3 = sp2(T/F)]]
			#"ancest.prop" is a ninds*2 matrix giving the ancestry proportion of each individual on each of their two copies of the chromosome.

* `_simsums_chunks.Robj` : has a single object, "chunks", 
    - which is a list with components:
        - selected : at the selected site
        - distant : at recombination distance 0.01
        - unlinked : at the center of the other chromosome
    - and each of these is a list with one component per deme, each of which is
		- a matrix of chunk sizes, of dimensions deme_size by ploidy, with columns named e.g. "X1" and "X2" for the chromosomes.
        - Each matrix entry gives the chunk size of Ancestry B surrounding the locus (0 means the locus is of Ancestry A)


# To make more plots, e.g. of haplotypes, of the results

Many plots can be produced by `haplotype-plots.R`;
this runs on the `sims.sums` object saved out by the script in the file `results_runid_%s_simsums.Robj`.
The script `make-plots.R` will load the required things and run `haplotype-plots.R`;
it works as follows:
```
get_simparams <- function (fname) {
    strings <- c( SIGMA="SIGMA", ninds="Ninds", ndemes="ndemes", S="s", tau="tau", run.id="runid_" )
    out <- lapply( strings, function (x) { as.numeric(gsub(sprintf(".*[/_]%s([0-9.]+)_.*",x),"\\1",fname)) } )
    out$zone_age <- out$tau
    out$deme_size <- out$ninds/out$ndemes
    return(out)
}

simdir <- "simulation_SIGMA3_Ninds250000_ndemes500_s0.05_dir/tau100_dir"
simsum.files <- list.files(simdir, ".*simsums.Robj", full.names=TRUE)
sapply( simsum.files, function (fname) {
        params <- get_simparams(fname)
        load(fname)
        source("haplotype-plots.R")
    } )
```

# To compile the resulting figure across time points at the same parameter values

Since each set of parameter values gets a directory,
with time points in subdirectories,
we compile all the figures of each sort into multipage pdfs to put in the parameter value directory.
This is done with `compile-pdfs.sh`: 
change to the directory corresponding to parameter values, then run `compile-pdfs.sh`, e.g.:
```
cd simulation_SIGMA3_Ninds6000_ndemes120_s0.001_dir
../compile-pdfs.sh
```
This creates the files:

- `ratioAdjacentBlocksAlongChromAncBConditioning_nonnormalized.pdf`
- `blocksAlongChromNoConditioning_nonnormalized.pdf`
- `blocksAlongChromAncBConditioning_nonnormalized.pdf`
- `adjacentBlocksAlongChromNoConditioning_nonnormalized.pdf`
- `adjacentBlocksAlongChromAncBConditioning_nonnormalized.pdf`
- `LD.pdf`
- `freqplot.pdf`
- `mean.pdf`
- `ecdf.pdf`
- `density.pdf`
- `ratioAdjacentBlocksAlongChromNoConditioning.pdf`
- `ratioAdjacentBlocksAlongChromHeatmapNoConditioning.pdf`
- `ratioAdjacentBlocksAlongChromHeatmapAncBConditioning.pdf`
- `ratioAdjacentBlocksAlongChromAncBConditioning.pdf`
- `blocksAlongChromNoConditioning.pdf`
- `blocksAlongChromHeatmap.pdf`
- `blocksAlongChromHeatmapAncBConditioning.pdf`
- `blocksAlongChromAncBConditioning.pdf`
- `adjacentBlocksAlongChromNoConditioning.pdf`
- `adjacentBlocksAlongChromAncBConditioning.pdf`

# To create the "comparison to theory" document

First `devtools::install_github("petrelharp/templater")` and then run e.g.
```
make simulation_SIGMA1_Ninds5000_ndemes100_s0.1_dir/tau500_dir/results_comparison-to-theory.html
```
or even
```
SIMDIR=simulation_SIGMA1_Ninds5000_ndemes100_s0.1_dir/
for x in $(find $SIMDIR -name "*_simsums_chunks.Robj")
do 
    y=$(echo $x | sed -e 's/_simsums_chunks.Robj/_comparison-to-theory.html/') 
    make $y
done
```

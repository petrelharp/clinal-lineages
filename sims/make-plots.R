# run haplotype-plots.R
# on all _simsums.Robj files passed in on the command line

the.files <- commandArgs(TRUE)

get_simparams <- function (fname) {
    strings <- c( SIGMA="SIGMA", ninds="Ninds", ndemes="ndemes", S="s", tau="tau", run.id="runid_" )
    out <- lapply( strings, function (x) { as.numeric(gsub(sprintf(".*[/_]%s([0-9.]+)_.*",x),"\\1",fname)) } )
    out$zone_age <- out$tau
    out$deme_size <- out$ninds/out$ndemes
    return(out)
}

for (fname in the.files) {
    if (!file.exists(fname)) {
        warn(sprintf("File %s does not exist.\n",fname))
        break
    }
    params <- get_simparams(fname)
    load(fname)
    source("haplotype-plots.R")
}

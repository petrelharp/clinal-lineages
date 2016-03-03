# run this in a directory to compile information from all *intervalSizes.Robj and *intervalSizes_allAncs.Robj files within that directory.

source("sim-fns.R")

B.files <- list.files(".",".*intervalSizes.Robj",recursive=TRUE,full.names=TRUE)
all.files <- list.files(".",".*intervalSizes_allAncs.Robj",recursive=TRUE,full.names=TRUE)

B.outfile <- "mean-haplotype-lengths_AncB.csv"
all.outfile <- "mean-haplotype-lengths_AllAnc.csv"

get_simparams <- function (fname) {
    strings <- c( SIGMA="SIGMA", ninds="Ninds", ndemes="ndemes", S="s", tau="tau", run.id="runid_", start="start", stop="stop", by="by", chromosome="chr" )
    out <- lapply( strings, function (x) { as.numeric(gsub(sprintf(".*[/_]%s([0-9.]+)_.*",x),"\\1",fname)) } )
    out$zone_age <- out$tau
    out$deme_size <- out$ninds/out$ndemes
    return(out)
}

B.lengths <- lapply( B.files, function (filename) {
            load(filename)

            # this is a (locus x deme) matrix giving mean lengths
            #   for haplotype chunks of B ancestry
            # mean_per_deme_AncB = 
            do.call(rbind, lapply(intervalSizes,function(P){tapply(1:length(deme_ID),deme_ID,function(Z){all_sites = P[Z,]; mean(all_sites[which(all_sites>0)])})}))
        } )


all.lengths <- lapply( all.files, function (filename) {
            load(filename)
            # this is a (locus x deme) matrix giving mean lengths
            #  for haplotype chunks of ALL ANCESTRIES
            # mean_per_deme = 
            do.call(rbind, lapply(intervalSizes_allAncs,function(P){tapply(1:length(deme_ID),deme_ID,function(Z){mean(P[Z,])})}))
        } )

positions = with( list2env(params), seq(plot_start,plot_stop,increments) )

outstring = with( list2env(params), sprintf("simulation_SIGMA%s_Ninds%s_ndemes%s_s%s_dir/tau%s_dir/results_runid_%s_%s_start%g_stop%g_by%g",
                                            SIGMA, ninds, ndemes, S, zone_age, run.id,  
                                            chromosome, plot_start, plot_stop, increments ) )

deme_ID = rep( 1:params$ndemes, each=params$deme_size )



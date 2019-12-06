library( tidyverse )

Rcpp::sourceCpp("tanimoto.cpp")

main <- function()
{
    X <- readRDS( "data/morgan_normal.rds" ) %>% slice( 1:1e5 )

    Rprof("prof1.out")
    v <- map_dbl( X$fingerprint, hex_jaccard, X$fingerprint[1] )
    Rprof(NULL)
    summaryRprof("prof1.out")
}

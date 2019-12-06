library( tidyverse )

Rcpp::sourceCpp("tanimoto.cpp", rebuild=TRUE)

X <- readRDS( "data/morgan_normal.rds" ) %>% slice( 1:1e5 )
fp <- map( X$fingerprint, MorganFP$new )

Rprof()
res <- map_dbl( fp, fp[[1]]$tanimoto )
Rprof(NULL)

summaryRprof()

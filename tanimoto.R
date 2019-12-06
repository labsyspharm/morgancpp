library( tidyverse )

Rcpp::sourceCpp("tanimoto.cpp", rebuild=TRUE)

X <- readRDS( "data/morgan_normal.rds" )
fps <- MorganFPS$new( X$fingerprint )

Rprof("prof1.out")
res <- map_dbl( 1:nrow(X), fps$tanimoto, 1 )
Rprof(NULL)

Rprof("prof2.out")
res2 <- fps$tanimoto_all(1)
Rprof(NULL)
      
summaryRprof("prof1.out")
summaryRprof("prof2.out")

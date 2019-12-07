library( tidyverse )

Rcpp::sourceCpp("morgan.cpp", rebuild=TRUE)

## Load the data
fps <- readRDS( "data/morgan_normal.rds" )$fingerprint

## Similarity over hex strings
tanimoto( fps[1], fps[1] )
tanimoto( fps[1], fps[2] )

## Create a fingerprint collections
m <- MorganFPS$new( fps )
m$size()

## Similarity between pairs
m$tanimoto( 1, 1 )
m$tanimoto( 1, 2 )
m$tanimoto( 1, 2e6 )
m$tanimoto( 1, -1 )

## Similarity profile of a single compound
Rprof()
res <- m$tanimoto_all( 1 )
Rprof(NULL)
summaryRprof()

## Similarity profile of an external hex string
m$tanimoto_ext( "ABC" )
res2 <- m$tanimoto_ext( fps[1] )

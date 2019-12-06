library( tidyverse )

Rcpp::sourceCpp("tanimoto.cpp")

main <- function()
{
    X <- readRDS( "data/morgan_normal.rds" ) %>% slice( 1:1e6 )
    vhex <- X$fingerprint
    vbin <- map_chr( vhex, hex2bin )

    rbin <- bin_jaccard_many( vbin[1], vbin )
}

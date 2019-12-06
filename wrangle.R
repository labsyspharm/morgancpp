library( tidyverse )

main <- function()
{
    X <- read_rds( "data/canonical_fingerprints.rds" )
    X1 <- X %>% filter( fp_name == "morgan_normal" ) %>%
        pluck( "data", 1 ) %>% filter( fp_name == "morgan_normal" ) %>%
        select( lspci_id, fingerprint )
    saveRDS( X1, file="data/morgan_normal.rds" )

    X2 <- X1 %>% slice( 1:10000 )
    saveRDS( X2, file="data/morgan_normal_10k.rds" )
}

# morgancpp: Efficient structure for storing and comparing Morgan fingerprints

## Installation

The package can be installed directly from GitHub:

``` r
if( !require(devtools) ) install.packages("devtools")
devtools::install_github("ArtemSokolov/morgancpp")
```

Once installed, the package can be loaded through the standard `library()` interface:

``` r
library( morgancpp )
```

## Example data

The package works with Morgan fingerprints specified as hexadecimal strings of length 512.
An example set of 100,000 fingerprints is included with the package. It can be loaded as follows:

``` r
fn <- system.file( "examples/example1.txt.gz", package="morgancpp" )
fps <- scan( fn, what=character() )
# [1] "020000000000000000000100...
# [2] "020004000000000000808101...
# [3] "200000000000000000000000...
```

## Computing similarity of fingerprints

The tanimoto similarity can be computed directly on hexadecimal strings:

``` r
tanimoto( fps[1], fps[1] )    # 1
tanimoto( fps[1], fps[2] )    # 0.1627907
tanimoto( fps[1], "FFF" )
# Error in tanimoto(fps[1], "FFF") : Input hex string must be of length 512
```

However, doing this for a large number of fingerprints will be slow. Instead, we can initialize the c++ data structure to store fingerprints in an effcient manner and query that structure for similarity:

``` r
m <- MorganFPS$new( fps )
# C++ object <0x559fa18f2c40> of class 'MorganFPS' <0x559fa484d570>

object.size(fps)
57467704 bytes

m$size()
[1] 25600000
```

Similarity between fingerprints in the structure can be computed by indexing:

``` r
m$tanimoto( 1, 1 )   # 1
m$tanimoto( 1, 2 )   # 0.1627907
```

The entire similarity profile against all other fingerprints in the collection can be obtained for compounds already in the collection (specified as index `i`) or new external compounds (specified as hexadecimal strings):

``` r
res1 <- m$tanimoto_all( 1 )       # Compound 1 against all 100,000 fingerprints
res2 <- m$tanimoto_ext( fps[1] )  # External compound against all 100,000 fingerprints
```

## Additional documentation

Obtaining additional information about functionality of the package can be done through the standard R interface:

``` r
# All available functions
library( help = morgancpp )

# Help for individual functions / data structures
?morgancpp::tanimoto
?morgancpp::MorganFPS

# C++ docstrings providing full signatures of each method
morgancpp::MorganFPS
```

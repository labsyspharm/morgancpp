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


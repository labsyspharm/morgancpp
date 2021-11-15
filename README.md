morgancpp: Efficient structure for storing and comparing Morgan
fingerprints
================
4/27/2020

## Installation

The package can be installed directly from GitHub:

``` r
if(!require(devtools)) install.packages("devtools")
devtools::install_github("labsyspharm/morgancpp")
```

Once installed, the package can be loaded through the standard
`library()` interface:

``` r
library(morgancpp)
```

## Example data

The package works with Morgan fingerprints specified as hexadecimal
strings of length 512. An example set of 100,000 fingerprints is
included with the package. It can be loaded as follows:

``` r
fn <- system.file("examples/example1.txt.gz", package = "morgancpp")
fps <- scan(fn, what = character())
```

    ## [1] "02000000000000000000010000100000000000000000000000000000020000000000000802000000..."
    ## [2] "02000400000000000080810100008000000000000000000000000040082000000400000800000000..."
    ## [3] "20000000000000000000000000000000000000000000000040000000000000080020000000000020..."

## Computing similarity of fingerprints

The tanimoto similarity can be computed directly on hexadecimal strings:

``` r
tanimoto(fps[1], fps[1])
```

    ## [1] 1

``` r
tanimoto(fps[1], fps[2])
```

    ## [1] 0.1627907

``` r
tanimoto(fps[1], "FFF")
```

    ## Error in tanimoto(fps[1], "FFF"): Input hex string must be of length 512

However, doing this for a large number of fingerprints will be slow.
Instead, we can initialize the c++ data structure to store fingerprints
in an efficient manner and query that structure for similarity:

``` r
m <- MorganFPS$new(fps)
object.size(fps)
```

    ## 57467704 bytes

``` r
m$size()
```

    ## [1] 25600000

Similarity between fingerprints in the structure can be computed by
indexing:

``` r
m$tanimoto(1, 1)
```

    ## [1] 1

``` r
m$tanimoto(1, 2)
```

    ## [1] 0.1627907

The entire similarity profile against all other fingerprints in the
collection can be obtained for compounds already in the collection
(specified as index `i`) or new external compounds (specified as
hexadecimal strings):

``` r
res1 <- m$tanimoto_all( 1 )       # Compound 1 against all 100,000 fingerprints
```

|  id | structural_similarity |
|----:|----------------------:|
|   1 |             1.0000000 |
|   2 |             0.1627907 |
|   3 |             0.0957447 |
|   4 |             0.0930233 |
|   5 |             0.2179487 |
|   6 |             0.0750000 |

``` r
res2 <- m$tanimoto_ext( fps[1] )  # External compound against all 100,000 fingerprints
```

|  id | structural_similarity |
|----:|----------------------:|
|   1 |             1.0000000 |
|   2 |             0.1627907 |
|   3 |             0.0957447 |
|   4 |             0.0930233 |
|   5 |             0.2179487 |
|   6 |             0.0750000 |

## Tresholded searches

Often the similarities between all NxN pairs of fingerprints above a
certain threshold are of interest. The `tanimoto_threshold()` computes
these similarities and outputs them as a data frame.

``` r
res3 <- m$tanimoto_threshold( 0.9 )
```

    ## 0 done
    ## 10000 done
    ## 20000 done
    ## 30000 done
    ## 40000 done
    ## 50000 done
    ## 60000 done
    ## 70000 done
    ## 80000 done
    ## 90000 done

| id_1 |  id_2 | structural_similarity |
|-----:|------:|----------------------:|
|    6 | 36184 |             1.0000000 |
|   51 | 78113 |             0.9305556 |
|  122 |  6768 |             0.9062500 |
|  185 | 93177 |             0.9850746 |
|  209 | 88214 |             0.9137931 |
|  229 | 29896 |             0.9166667 |

## Saving fingerprints

Fingerprints can be stored on disk in a binary format supporting
extremely fast writing and reading performance using
[Zstandard](https://github.com/facebook/zstd) compression.

``` r
tmp_file <- tempfile()

m$save_file(tmp_file)
```

    ## Wrinting 100000 fingerprints
    ## Fingerprints compressed 6264569 bytes
    ## Wrote fingerprints
    ## Names compressed 253792 bytes
    ## Wrote Names

``` r
m2 <- MorganFPS$new(tmp_file, from_file = TRUE)
```

    ## Reading 100000 fingerprints from file
    ## Fingerprint block has 6264569 bytes
    ## Fingerprints decompressed
    ## Names block has 253792 bytes
    ## Names decompressed

``` r
m2$size()
```

    ## [1] 25600000

## Additional documentation

Obtaining additional information about functionality of the package can
be done through the standard R interface:

``` r
# All available functions
library( help = morgancpp )

# Help for individual functions / data structures
?morgancpp::tanimoto
?morgancpp::MorganFPS

# C++ docstrings providing full signatures of each method
morgancpp::MorganFPS
```

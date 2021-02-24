context( "MorganFPS functionality" )

load_example1 <- function(n)
    scan( "../../inst/examples/example1.txt.gz", n=n,
         what=character(), quiet=TRUE )

test_that("Self-similarity is always 1", {
    ## Test the first 100 strings from the example
    v <- load_example1(1000)
    res <- sapply( v, function(hx) tanimoto(hx,hx) )
    expect_equal( unname(res), rep(1,1000) )
})

test_that("Hex strings can be compared directly", {
    ## Spot-check several pairwise values
    v <- load_example1(10)
    expect_equal( tanimoto(v[1], v[2]), 0.1627907 )
    expect_equal( tanimoto(v[5], v[6]), 0.08450704 )
    expect_equal( tanimoto(v[9], v[10]), 0.09677419 )
})

test_that("Hex strings have to be of length 512", {
    expect_error( tanimoto("ABC", "123"),
                 "Input hex string must be of length 512" )
})

test_that("New collections can be instatiated from hex strings", {
    v <- load_example1(1000)
    m <- MorganFPS$new(v)
    expect_equal( m$size(), 256000 )
})

test_that("Collections can be queried for pairwise similarities", {
    v <- load_example1(10)
    m <- MorganFPS$new(v)

    ## Compare to stand-alone function
    expect_equal( m$tanimoto(1,2), tanimoto(v[1], v[2]) )
    expect_equal( m$tanimoto(5,6), tanimoto(v[5], v[6]) )
    expect_equal( m$tanimoto(9,10), tanimoto(v[9], v[10]) )

    ## Index has to be in-range
    expect_error( m$tanimoto(-1,1) )
})

test_that("Collections can be queried for full similarity profiles", {
    v <- load_example1(100)
    m <- MorganFPS$new(v)
    v0 <- sapply( 1:100, function(i) m$tanimoto(1,i) )
    v1 <- m$tanimoto_all(1)
    v2 <- m$tanimoto_ext(v[1])

    expect_equal( nrow(v1), 100 )
    expect_equal( nrow(v2), 100 )
    expect_identical( v0, v1[["structural_similarity"]] )
    expect_identical( v0, v2[["structural_similarity"]] )
})

test_that("Collection indexing is 1-based", {
    v <- load_example1(1000)
    m <- MorganFPS$new(v)

    ## Pair-wise function
    expect_error( m$tanimoto(0,1), "not found" )
    expect_error( m$tanimoto(1,0), "not found" )
    expect_error( m$tanimoto(-1,1), "not found" )
    expect_error( m$tanimoto(1,-1), "not found" )
    expect_error( m$tanimoto(1001,1), "not found" )
    expect_error( m$tanimoto(1,1001), "not found" )

    expect_identical( m$tanimoto(1000,1000), 1 )

    ## Full-profile function
    expect_error( m$tanimoto_all(-1), "not found" )
    expect_error( m$tanimoto_all(0), "not found" )
    expect_error( m$tanimoto_all(1001), "not found" )

    expect_equal( nrow(m$tanimoto_all(1000)), 1000 )
})

test_that("Collections can be saved to and loaded from binary files", {
    v <- load_example1(100)
    m <- MorganFPS$new(v)
    tmp <- tempfile()
    m$save_file(tmp)
    browser()
    expect_gt(file.size(tmp), 5000)
    m2 <- MorganFPS$new(tmp, from_file = TRUE)
    expect_equal(m2$size(), 25600)
    expect_equal(m$tanimoto_all(1), m2$tanimoto_all(1))
})

test_that("Fingerprint ids are respected", {
    set.seed(42)
    v <- load_example1(100)
    vn <- sample(1:1000, 100)
    names(v) <- vn
    m <- MorganFPS$new(v)
    res <- m$tanimoto_all(vn[[1]])
    expect_equal(res[["id"]], sort(vn))
})

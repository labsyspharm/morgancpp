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
    expect_identical( v0, v1[["similarity"]] )
    expect_identical( v0, v2[["similarity"]] )
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
    expect_gt(file.size(tmp), 5000)
    m2 <- MorganFPS$new(tmp, from_file = TRUE)
    expect_equal(m2$size(), 25600)
    expect_equal(m$tanimoto_all(1), m2$tanimoto_all(1))
})

test_that("Fingerprint ids are respected", {
    set.seed(42)
    v <- load_example1(100)
    vn <- c(1e09L, 1e03L, sample(1e09L, 98))
    vn_s <- sort(vn)
    names(v) <- vn
    m <- MorganFPS$new(v)
    res <- m$tanimoto_all(vn[[1]])
    expect_equal(res[["id"]], vn_s)
    names(v) <- as.character(vn)
    m <- MorganFPS$new(v)
    res <- m$tanimoto_all(vn[[1]])
    expect_equal(res[["id"]], vn_s)
    names(v) <- as.numeric(vn)
    m <- MorganFPS$new(v)
    res <- m$tanimoto_all(vn[[1]])
    expect_equal(res[["id"]], vn_s)
    expect_error(m$tanimoto_all("66"), "Fingerprint 66 not found")
    expect_error(m$tanimoto_subset(c("66", "67"), NULL), "Fingerprint 66 not found")
})

test_that("Duplicate fingerprints are rejected", {
    v <- load_example1(5)
    vn <- c(1, 2, 3, 5, 3)
    names(v) <- vn
    expect_error(MorganFPS$new(v), "Duplicate names are not allowed")
})

test_that("Subset queries are accurate", {
    set.seed(42)
    v <- load_example1(100)
    vn <- c(1e09L, 1e03L, sample(1e09L, 98))
    vn_s <- sort(vn)
    names(v) <- vn
    vn_x <- sample(vn, 10)
    vn_y <- sample(vn, 20)
    m <- MorganFPS$new(v)
    res <- m$tanimoto_subset(vn_x, vn_y)
    expect_equal(nrow(res), 200)
    combos <- expand.grid(sort(vn_y), sort(vn_x))
    manual_similarity <- with(
        combos,
        mapply(function(x, y) m$tanimoto(x, y), Var2, Var1)
    )
    expect_equal(res$similarity, manual_similarity)
    res <- m$tanimoto_subset(vn_x, NULL)
    expect_equal(nrow(res), 1000)
    combos <- expand.grid(sort(vn), sort(vn_x))
    manual_similarity <- with(
        combos,
        mapply(function(x, y) m$tanimoto(x, y), Var2, Var1)
    )
    expect_equal(res$similarity, manual_similarity)
})


test_that("Thresholded queries work", {
  v <- load_example1(10)
  m <- MorganFPS$new(v)
  res <- m$tanimoto_threshold(0.1)
  expect_equal(nrow(res), 18)
})

test_that("Identity matching works", {
  v <- load_example1(100)
  m <- MorganMap$new(v)
  set.seed(42)
  v2 <- sample(load_example1(300), 10)
  ma <- m$find_matches(v2)
  expect_equal(nrow(ma), 4)
})

test_that("Loading fingprints from rdkit hex code works", {
  fps <- fingerprints(
    c(
      "e0ffffff000400001a00000026582eb62a09002c2a620c18309638069c0644327826440c324e3e98",
      "e0ffffff0004000037000000043c221a0c0e0c044a140ede321e3e2c06000e104c14362c040c206a3a76022a06401c1800540402021e002a0a00183e1640101a462e1208",
      "e0ffffff000400003200000002042c04182a1c125a0006ae8e2a00140618060622620e365c92307e0a023c021478041e02320e0c04282046000e10124c3208"
    ),
    "rle"
  )
  m <- MorganFPS$new(fps)
  expect_equal(m$n(), 3)
  expect_equal(nrow(m$tanimoto_all(1)), 3)

  m <- MorganMap$new(fps)
  ma <- m$find_matches(
    fingerprints(
      c(
        "e0ffffff0004000037000000043c221a0c0e0c044a140ede321e3e2c06000e104c14362c040c206a3a76022a06401c1800540402021e002a0a00183e1640101a462e1208"),
      "rle"
    )
  )
  expect_equal(nrow(ma), 1)
})

test_that("Full and RLE encoded fingerprints are identical", {
  full <- c('00000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000800000000000000000000000000000000000000000000000000000000000000000120000002000000000100000000000000000000014000000000000000104080000000020000000000000000000000501001000040001000000000000000000020109000000000000000000000000000000000000044000400000000080000000000000000004000000000000000000000000000000000000000000000000000000000000000000000000200000000000010000800000000000040000000000400000000000000200100800000000000',
            '40000000001000000000810000000000000000000000000000000000000000000000080000000000100110000080002000000000000000000000800100400000000000000000000000000000000000000000000020000000010000000400000000000200028000000000000000000020000000020100800000000000000000102800001000000000400000200000000080c01000000000000000000000800000000000000004000008040000088000000000000000080000000000002000000100000000000000000000000020000000000000000020000000100200000000000000200000000000000048000040000000000000000000200000000000000000',
            '00000000001000000000810000000000000400000000000000000000600000000000080002000000100000020000002000000000000000000000800100400000000000000000000000400020000000000020000030000000010010000400000000000200020000000000000000800000000000120100000000000000000000002000000048000000400000200000000080801000000000000000000000800000000000000000000000040000088000000000000000080000080000022000000100000000000000000000000000000000000000000800000000000200000000000000204000000000000048000040001000000000000010000000000000000000')
  rle <- c(
    'e0ffffff0008000023000000de31022902043650ae02820a0c4ab4020e1622229c1008045500061e4c989102682068567818105e',
    'e0ffffff00080000330000000252480cf1005c16062822aa1c12e100483a601e02a2461010940c022c4a304a0e0012a8882c102c06866a38c48e4014767c0428a08a',
           'e0ffffff000800003600000056480c6aa6006222383436aa1c12be205e2e00461622601e82640410b43c0438304a0e14a8c82c06862e320638910062760c6e0428226e98'
    )
  m1 <- MorganFPS$new(full)
  m2 <- MorganFPS$new(fingerprints(rle, format = "rle"))
  expect_equal(m1$tanimoto_all(1)$similarity, m2$tanimoto_all(1)$similarity)
})

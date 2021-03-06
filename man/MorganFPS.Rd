% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{MorganFPS}
\alias{MorganFPS}
\title{Morgan fingerprints collection}
\description{
Efficient structure for storing a set of Morgan fingerprints
}
\section{Fields}{

\describe{
\item{\code{new}}{Constructor. Accepts either a vector of fingerprints in hexadecimal
format or a path to a binary file of fingerprints using the argument
`from_file = TRUE`. The vector of fingerprints can optionally be named.
Names need to be coercible to integers. When querying, the indices i and j
refer to the given names.}

\item{\code{tanimoto}}{(i,j) similarity between fingerprints i and j}

\item{\code{tanimoto_all}}{(i) similarity between fingerprint i and all others}

\item{\code{tanimoto_ext}}{(s) similarity between external hexadecimal string s and all
fingerprints in the collection}

\item{\code{save_file}}{(path, compression_level) Save fingerprints to file in binary format}

\item{\code{size}}{number of bytes used to store the fingerprints}
}}


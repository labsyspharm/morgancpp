#' @useDynLib morgancpp
NULL

Rcpp::loadModule("morgan_cpp", TRUE)
Rcpp::loadModule("morgan_identity_cpp", TRUE)

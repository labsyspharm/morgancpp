#include <Rcpp.h>
#include "zstd/zstd.h"
#include <array>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

#include "utils.hpp"

FingerprintName convert_name(std::string x) {
  FingerprintName name;
  try {
    // Converting to float first to cover cases where integers are passed as, e.g. 1e+07
    name = static_cast<FingerprintName>(std::stold(x));
  } catch (const std::invalid_argument& e) {
    Rcpp::stop(
      "Fingerprint names must be passed as positive integers, numerics, or strings representing integers."
    );
  } catch (...) {
    Rcpp::stop("Unknown error converting Fingerprint names");
  }
  return name;
}

FingerprintName convert_name(Rcpp::RObject& x) {
  FingerprintName name;
  if (Rcpp::is<Rcpp::CharacterVector>(x)) {
    name = convert_name(Rcpp::as<std::string>(x));
  } else if (Rcpp::is<Rcpp::IntegerVector>(x) || Rcpp::is<Rcpp::NumericVector>(x)) {
    name = static_cast<FingerprintName>(Rcpp::as<int>(x));
  } else {
    Rcpp::stop(
      "Fingerprint indices (names) must be passed as positive integers, numerics, or strings representing integers."
    );
  }
  return name;
}

std::vector<FingerprintName> convert_name_vec(Rcpp::RObject& names) {
  std::vector<FingerprintName> unsorted_names;
  if (Rcpp::is<Rcpp::CharacterVector>(names)) {
    auto passed_names_vec = Rcpp::as<Rcpp::CharacterVector>(names);
    unsorted_names.reserve(passed_names_vec.size());
    for (auto& x: passed_names_vec)
      unsorted_names.push_back(convert_name(Rcpp::as<std::string>(x)));
  } else if (Rcpp::is<Rcpp::IntegerVector>(names) || Rcpp::is<Rcpp::NumericVector>(names)) {
    unsorted_names = Rcpp::as< std::vector<FingerprintName> >(names);
  } else {
    Rcpp::stop(
      "Fingerprint names must be passed as positive integers, numerics, or strings representing integers."
    );
  }
  return unsorted_names;
}

std::vector<size_t> sort_indices(std::vector<FingerprintName>& unsorted_names) {
  size_t n = unsorted_names.size();
  // Make vector that contains indices that sort the given vector
  std::vector<size_t> sort_vector(n);
  std::iota(sort_vector.begin(), sort_vector.end(), 0);
  std::sort(
    sort_vector.begin(), sort_vector.end(),
    [&](const size_t i, const size_t j){return unsorted_names.at(i) < unsorted_names.at(j);}
  );
  return sort_vector;
}

std::vector<FingerprintName> convert_sort_name_vec(Rcpp::RObject& names) {
  auto unsorted_names = convert_name_vec(names);
  auto sort_vector = sort_indices(unsorted_names);
  std::vector<FingerprintName> sorted_names;
  sorted_names.reserve(unsorted_names.size());
  for (auto &i: sort_vector) {
    sorted_names.push_back(unsorted_names.at(i));
  }
  // Checking for duplicate names
  auto duplicate_pair = std::adjacent_find(sorted_names.begin(), sorted_names.end());
  if (duplicate_pair != sorted_names.end()) {
    Rcpp::stop("Duplicate names are not allowed");
  }
  return sorted_names;
}


// Parse a single hexadecimal character to its integer representation.
int parse_hex_char(const char& c) {
  int v;
  if ((c >= '0') && (c <= '9')) {
    v = (c - '0');
  } else if ((c >= 'A') && (c <= 'F')) {
    v = (c - 'A' + 10);
  } else if ((c >= 'a') && (c <= 'f')) {
    v = (c - 'a' + 10);
  } else {
    ::Rf_error("Hex string may only contain characters in [0-9A-F]");
  }
  return v;
}

// Convert raw byte string to fingerprint.
Fingerprint raw2fp(const std::string& raw) {
  if (raw.length() != 256) {
    ::Rf_error("Input raw string must be of length 256");
  }
  Fingerprint fp;
  std::memcpy(fp.data(), &raw[0], sizeof(fp));
  return fp;
}

// Convert ASCII hex string to fingerprint.
Fingerprint hex2fp(const std::string& hex) {
  if (hex.length() != 512) {
    ::Rf_error("Input hex string must be of length 512");
  }
  // Convert hex to raw bytes
  std::string raw(256, '\0');
  auto hi = std::cbegin(hex);
  for (auto& ri : raw) {
    ri = parse_hex_char(*hi++) | (parse_hex_char(*hi++) << 4);
  }
  return raw2fp(raw);
}


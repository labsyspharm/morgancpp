#include <Rcpp.h>
#include <zlib.h>
#include <array>
#include <fstream>
#include <string>
#include <algorithm>

using Fingerprint = std::array<std::uint64_t, 32>;
using FingerprintName = std::int32_t;

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


// Compute Jaccard similarity of two fingerprints
double jaccard_fp(const Fingerprint& f1, const Fingerprint& f2) {
  int count_and = 0, count_or = 0;
  auto i1 = std::cbegin(f1), i2 = std::cbegin(f2);
  while (i1 != std::cend(f1)) {
    count_and += __builtin_popcountll(*i1 & *i2);
    count_or += __builtin_popcountll(*i1 | *i2);
    ++i1;
    ++i2;
  }
  return static_cast<double>(count_and) / count_or;
}


//' Tanimoto similarity between two Morgan fingerprints
//'
//' Computes Tanimoto similarity between two hexadecimal strings
//'
//' @param s1 Hexadecimal string of length 512
//' @param s2 Hexadecimal string of length 512
//' @return Jaccard similarity over the bits representing individual keys
//' @export
// [[Rcpp::export]]
double tanimoto(const std::string& s1, const std::string& s2) {
  const Fingerprint& fp1 = hex2fp(s1);
  const Fingerprint& fp2 = hex2fp(s2);
  return jaccard_fp(fp1, fp2);
}


//' @name MorganFPS
//' @title Morgan fingerprints collection
//' @description Efficient structure for storing a set of Morgan fingerprints
//' @field new Constructor. Accepts either a vector of fingerprints in hexadecimal
//'   format or a path to a binary file of fingerprints using the argument
//'   `from_file = TRUE`
//' @field tanimoto (i,j) similarity between fingerprints i and j
//' @field tanimoto_all (i) similarity between fingerprint i and all others
//' @field tanimoto_ext (s) similarity between external hexadecimal string s and all
//'    fingerprints in the collection
//' @field save_file (path) Save fingerprints to file in binary format
//' @field size number of bytes used to store the fingerprints
//' @importFrom Rcpp cpp_object_initializer
//' @export
class MorganFPS {

public:

  // Constructor accepts a named character vector of hex strings
  MorganFPS(Rcpp::CharacterVector fps_hex) {
    size_t n = fps_hex.length();
    Rcpp::RObject passed_names = fps_hex.names();
    fps.reserve(n);
    fp_names.reserve(n);
    if(passed_names.isNULL()) {
      for (FingerprintName i = 0; i < n; i++) {
        fp_names.push_back(i + 1);
      }
    } else {
      Rcpp::CharacterVector passed_names_vec = Rcpp::as<Rcpp::CharacterVector>(passed_names);
      if (Rcpp::unique(passed_names_vec).length() != n)
        Rcpp::stop("Names must be unique");
      for (Rcpp::CharacterVector::iterator i = passed_names_vec.begin(); i != passed_names_vec.end(); i++) {
        fp_names.push_back(std::stoll(std::string(*i)));
      }
      Rcpp::IntegerVector idx = Rcpp::seq_along(passed_names_vec) - 1;
      std::sort(idx.begin(), idx.end(), [&](int i, int j){return passed_names_vec[i] < passed_names_vec[j];});
      std::vector<FingerprintName> fp_names_sorted(n);
      for (auto &i : idx) {
        fp_names_sorted.push_back(fp_names[i]);
      }
      fp_names = fp_names_sorted;
      fps_hex = fps_hex[idx];
    }
    for (Rcpp::CharacterVector::iterator i = fps_hex.begin(); i != fps_hex.end(); i++) {
      fps.push_back(hex2fp(std::string(*i)));
    }
  }

  // Constructor accepts a file path to load fingerprints from binary file
  MorganFPS(const std::string& filename, const bool from_file) {
    gzFile in_stream = gzopen(filename.c_str(), "rb");
    if (!in_stream) {
      Rcpp::stop("gzopen of " + filename + " failed: " + strerror(errno));
    }
    size_t n;
    gzread(in_stream, reinterpret_cast<char*>(&n), sizeof(size_t));
    fps.resize(n);
    fp_names.resize(n);
    int bytes_read;
    bytes_read = gzread(in_stream, reinterpret_cast<char*>(fps.data()), n * sizeof(Fingerprint));
    if (bytes_read != n * sizeof(Fingerprint)) {
      gzclose(in_stream);
      Rcpp::stop("Error reading gzipped fingerprints. Bad file?");
    }
    bytes_read = gzread(in_stream, reinterpret_cast<char*>(fp_names.data()), n * sizeof(FingerprintName));
    if ((bytes_read != n * sizeof(FingerprintName))) {
      gzclose(in_stream);
      Rcpp::stop("Error reading gzipped fingerprint names. Bad file?");
    }
    gzclose(in_stream);
  }

  // Tanimoto similarity between drugs i and j
  double tanimoto(FingerprintName &i, FingerprintName &j) {
    return jaccard_fp(fp_index(i), fp_index(j));
  }

  // Tanimoto similarity of drug i to every other drug
  Rcpp::DataFrame tanimoto_all(FingerprintName &i) {
    const Fingerprint& fp_other = fp_index(i);
    Rcpp::NumericVector res(fps.size());
    for (int i = 0; i < fps.size(); i++) {
      res[i] = jaccard_fp(fps[i], fp_other);
    }
    return Rcpp::DataFrame::create(
      Rcpp::Named("id") = fp_names,
      Rcpp::Named("structural_similarity") = res
    );
  }

  // Tanimoto similarity of an external drug to every other drug
  //   in the collection
  Rcpp::DataFrame tanimoto_ext(const std::string& other) {
    const Fingerprint& fp_other = hex2fp(other);
    Rcpp::NumericVector res(fps.size());
    for (int i = 0; i < fps.size(); i++) {
      res[i] = jaccard_fp(fps[i], fp_other);
    }
    return Rcpp::DataFrame::create(
      Rcpp::Named("id") = fp_names,
      Rcpp::Named("structural_similarity") = res
    );
  }

  // Save binary fp file
  void save_file(const std::string& filename) {
    gzFile out_stream = gzopen(filename.c_str(), "wb8");
    size_t n = fps.size();
    gzwrite(out_stream, reinterpret_cast<char*>(&n), sizeof(size_t));
    gzwrite(out_stream, reinterpret_cast<char*>(fps.data()), size());
    gzwrite(out_stream, reinterpret_cast<char*>(fp_names.data()), fps.size() * sizeof(FingerprintName));
    gzclose(out_stream);
  }

  // Size of the dataset
  int size() {
    return fps.size() * sizeof(Fingerprint);
  }

private:

  Fingerprint& fp_index(FingerprintName &x) {
    auto fp_pt = std::lower_bound(fp_names.begin(), fp_names.end(), x);
    if (*fp_pt != x)
      Rcpp::stop("Fingerprint %s not found", std::to_string(x).c_str());
    return fps.at(fp_pt - fp_names.begin());
  }

  std::vector<Fingerprint> fps;
  std::vector<FingerprintName> fp_names;
};

// Expose all relevant classes through an Rcpp module
RCPP_EXPOSED_CLASS(MorganFPS)
RCPP_MODULE(morgan_cpp) {

  using namespace Rcpp;

  class_<MorganFPS>( "MorganFPS" )
    .constructor< CharacterVector >("Construct fingerprint collection from vector of fingerprints")
    .constructor< std::string, bool >("Construct fingerprint collection from binary file")
    .method("size", &MorganFPS::size, "Size of the data in bytes")
    .method("tanimoto", &MorganFPS::tanimoto,
	    "Similarity between two fingerprints in the collection")
    .method("tanimoto_all", &MorganFPS::tanimoto_all,
	    "Similarity of a fingerprint against all other fingerprints in the collection")
    .method("tanimoto_ext", &MorganFPS::tanimoto_ext,
	    "Similarity of an external fingerprints against all fingerprints in the collection")
    .method("save_file", &MorganFPS::save_file,
	    "Save fingerprints to file in binary format")
    ;
}

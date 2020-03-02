#include <Rcpp.h>
#include <zlib.h>
#include <array>
#include <fstream>
#include <string>

using Fingerprint = std::array<std::uint64_t, 32>;

// Read 1k records at a time
#define FP_READ_BUFFER sizeof(Fingerprint) * 1000
#define FP_NAME_READ_BUFFER sizeof(int) * 1000

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
// //' @field save_file (path) Save fingerprints to file in binary format
//' @field size number of bytes used to store the fingerprints
//' @export
class MorganFPS {

public:

  // Constructor accepts a named character vector of hex strings
  MorganFPS(Rcpp::CharacterVector fps_hex) {
    Rcpp::RObject passed_names = fps_hex.names();
    fps.reserve(fps_hex.length());
    fp_names.reserve(fps_hex.length());
    if(passed_names.isNULL()) {
      for (int i = 0; i < fps_hex.length(); i++) {
        fp_names.push_back(i + 1);
      }
    } else {
      Rcpp::CharacterVector passed_names_vec = Rcpp::as<Rcpp::CharacterVector>(passed_names);
      for (Rcpp::CharacterVector::iterator i = passed_names_vec.begin(); i != passed_names_vec.end(); i++) {
        fp_names.push_back(std::stoi(std::string(*i)));
      }
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
    int n;
    gzread(in_stream, reinterpret_cast<char*>(&n), sizeof(int));
    fps.reserve(n);
    fp_names.reserve(n);
    int ix = 0;
    int bytes_read;
    char buffer[FP_READ_BUFFER];
    while (1) {
      bytes_read = gzread(in_stream, buffer, FP_READ_BUFFER);
      std::copy(buffer, buffer + bytes_read, reinterpret_cast<char*>(fps.data()) + ix);
      ix = ix + bytes_read;
      if (bytes_read < FP_READ_BUFFER)
        break;
    }
    ix = 0;
    char buffer_names[FP_NAME_READ_BUFFER];
    while (1) {
      bytes_read = gzread(in_stream, buffer_names, FP_NAME_READ_BUFFER);
      std::copy(buffer_names, buffer_names + bytes_read, reinterpret_cast<char*>(fp_names.data()) + ix);
      ix = ix + bytes_read;
      if (bytes_read < FP_READ_BUFFER) {
        if (gzeof(in_stream)) break;
        else {
          gzclose(in_stream);
          Rcpp::stop("Error reading gzipped fingerprints. Bad file?");
        }
      }
    }
    gzclose(in_stream);
  }

  // Tanimoto similarity between drugs i and j
  double tanimoto(int i, int j) {
    return jaccard_fp(fps.at(i - 1), fps.at(j - 1));
  }

  // Tanimoto similarity of drug i to every other drug
  Rcpp::NumericVector tanimoto_all(int i) {
    const Fingerprint& fp_other = fps.at(i - 1);
    Rcpp::NumericVector res(fps.size());
    res.names() = fp_names;
    for (int i = 0; i < fps.size(); i++) {
      res[i] = jaccard_fp(fps[i], fp_other);
    }
    return res;
  }

  // Tanimoto similarity of an external drug to every other drug
  //   in the collection
  Rcpp::NumericVector tanimoto_ext(const std::string& other) {
    const Fingerprint& fp_other = hex2fp(other);
    Rcpp::NumericVector res(fps.size());
    res.names() = fp_names;
    for (int i = 0; i < fps.size(); i++) {
      res[i] = jaccard_fp(fps[i], fp_other);
    }
    return res;
  }

  // Save binary fp file
  void save_file(const std::string& filename) {
    gzFile out_stream = gzopen(filename.c_str(), "wb8");
    int n = fps.size();
    gzwrite(out_stream, reinterpret_cast<char*>(&n), sizeof(int));
    gzwrite(out_stream, reinterpret_cast<char*>(fp_names.data()), fps.size() * sizeof(int));
    gzwrite(out_stream, reinterpret_cast<char*>(fps.data()), size());
    gzclose(out_stream);
  }

  // Size of the dataset
  int size() {
    return fps.size() * sizeof(Fingerprint);
  }

private:

  std::vector<Fingerprint> fps;
  std::vector<int> fp_names;
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

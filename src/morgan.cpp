#include <Rcpp.h>
#include <array>

using Fingerprint = std::array<std::uint64_t, 32>;


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
//' @field new Constructor
//' @field tanimoto (i,j) similarity between fingerprints i and j
//' @field tanimoto_all (i) similarity between fingerprint i and all others
//' @field tanimoto_ext (s) similarity between external hexadecimal string s and all
//'    fingerprints in the collection
//' @field size number of bytes used to store the fingerprints
//' @export
class MorganFPS {

public:

  // Constructor accepts a character vector of hex strings
  MorganFPS(const std::vector<std::string>& vhx) {
    for (auto&& hex : vhx) {
      fps.push_back(hex2fp(hex));
    }
  }

  // Tanimoto similarity between drugs i and j
  // Adjust for 0-based indexing
  double tanimoto(int i, int j) {
    return jaccard_fp(fps.at(i-1), fps.at(j-1));
  }

  // Tanimoto similarity of drug i to every other drug
  std::vector<double> tanimoto_all(int i) {
    const Fingerprint& fp_other = fps.at(i-1);
    std::vector<double> res;
    res.reserve(fps.size());
    for (auto&& fp : fps) {
      res.push_back(jaccard_fp(fp, fp_other));
    }
    return res;
  }

  // Tanimoto similarity of an external drug to every other drug
  //   in the collection
  std::vector<double> tanimoto_ext(const std::string& other) {
    const Fingerprint& fp_other = hex2fp(other);
    std::vector<double> res;
    res.reserve(fps.size());
    for (auto&& fp : fps) {
      res.push_back(jaccard_fp(fp, fp_other));
    }
    return res;
  }

  // Size of the dataset
  int size() {
    return fps.size() * sizeof(Fingerprint);
  }

private:

  std::vector<Fingerprint> fps;

};


// Expose all relevant classes through an Rcpp module
RCPP_EXPOSED_CLASS(MorganFPS)
RCPP_MODULE(morgan_cpp) {

  using namespace Rcpp;

  class_<MorganFPS>( "MorganFPS" )
    .constructor< std::vector<std::string> >()
    .method("size", &MorganFPS::size, "Size of the data in bytes")
    .method("tanimoto", &MorganFPS::tanimoto,
	    "Similarity between two fingerprints in the collection")
    .method("tanimoto_all", &MorganFPS::tanimoto_all,
	    "Similarity of a fingerprint against all other fingerprints in the collection")
    .method("tanimoto_ext", &MorganFPS::tanimoto_ext,
	    "Similarity of an external fingerprints against all fingerprints in the collection")
    ;
}

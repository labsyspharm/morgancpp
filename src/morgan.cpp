#include <Rcpp.h>
#include "zstd/zstd.h"
#include <array>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

using Fingerprint = std::array<std::uint64_t, 32>;
using FingerprintName = std::int32_t;
using FingerprintN = std::uint64_t;

FingerprintName convert_name(std::string x) {
  unsigned long long temp;
  try {
    temp = std::stoll(x);
  } catch (const std::invalid_argument& e) {
      Rcpp::stop("Fingerprint names must be passed as strings of positive integers, e.g. c(\"1\", \"2\")");
  } catch (...) {
      Rcpp::stop("Unknown error converting Fingerprint names");
  }
  FingerprintName temp2 = (FingerprintName) temp;
  return temp2;
}

FingerprintName convert_name(Rcpp::RObject& x) {
  FingerprintName name;
  if (Rcpp::is<Rcpp::CharacterVector>(x)) {
    name = convert_name(Rcpp::as<std::string>(x));
  } else if (Rcpp::is<Rcpp::IntegerVector>(x) || Rcpp::is<Rcpp::NumericVector>(x)) {
    name = (FingerprintName) Rcpp::as<int>(x);
  } else {
    Rcpp::stop(
      "Fingerprint indices (names) must be passed as positive integers, numerics, or strings representing integers."
    );
  }
  return name;
}

size_t zstd_frame_decompress(
    std::ifstream &in_stream, size_t &compressed_size, std::vector<char> &out_buffer
) {
  std::vector<char> compressed_buffer;
  // compressed_buffer.resize(ZSTD_FRAMEHEADERSIZE_MAX);
  // std::streampos cur_pos = in_stream.tellg();
  // in_stream.read(compressed_buffer.data(), ZSTD_FRAMEHEADERSIZE_MAX);
  // in_stream.seekg(cur_pos, std::ios_base::beg);
  //
  // unsigned long long const decompressed_size = ZSTD_getFrameContentSize(
  //   compressed_buffer.data(), ZSTD_FRAMEHEADERSIZE_MAX
  // );
  compressed_buffer.resize(compressed_size);
  in_stream.read(compressed_buffer.data(), compressed_size);
  unsigned long long const decompressed_size = ZSTD_getFrameContentSize(
    compressed_buffer.data(), compressed_size
  );
  if (decompressed_size == ZSTD_CONTENTSIZE_ERROR || decompressed_size == ZSTD_CONTENTSIZE_UNKNOWN) {
    Rcpp::stop("Error finding decompressed frame size");
  }

  size_t const frame_compressed_size = ZSTD_findFrameCompressedSize(
    compressed_buffer.data(), compressed_size
  );
  if (ZSTD_isError(frame_compressed_size)) {
    Rcpp::stop("Error finding compressed frame size: %s", ZSTD_getErrorName(frame_compressed_size));
  }
  if (compressed_size != frame_compressed_size) {
    Rcpp::stop("Inconsistent reported compressed sizes: %i and %i", compressed_size, frame_compressed_size);
  }

  out_buffer.resize(decompressed_size);
  size_t const decompressed_bytes = ZSTD_decompress(
    out_buffer.data(), decompressed_size,
    compressed_buffer.data(), compressed_size
  );
  if (ZSTD_isError(decompressed_bytes)) {
    Rcpp::stop("Error decompressing: %s", ZSTD_getErrorName(decompressed_bytes));
  }
  if (decompressed_bytes != decompressed_size) {
    Rcpp::stop("Inconsistent decompressed size: Expected %i Actual %i", decompressed_size, decompressed_bytes);
  }

  return decompressed_size;
};

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
//'   `from_file = TRUE`. The vector of fingerprints can optionally be named.
//'   Names need to be coercible to integers. When querying, the indices i and j
//'   refer to the given names.
//' @field tanimoto (i,j) similarity between fingerprints i and j
//' @field tanimoto_all (i) similarity between fingerprint i and all others
//' @field tanimoto_ext (s) similarity between external hexadecimal string s and all
//'    fingerprints in the collection
//' @field save_file (path, compression_level) Save fingerprints to file in binary format
//' @field size number of bytes used to store the fingerprints
//' @importFrom Rcpp cpp_object_initializer
//' @export
class MorganFPS {

public:

  // Constructor accepts a named character vector of hex strings
  MorganFPS(const Rcpp::CharacterVector& fps_hex) {
    size_t n = fps_hex.length();
    Rcpp::RObject passed_names = fps_hex.names();
    fps.reserve(n);
    fp_names.reserve(n);
    if(passed_names.isNULL()) {
      for (int i = 0; i < n; i++) {
        fp_names.push_back((FingerprintName) i + 1);
        fps.push_back(hex2fp(Rcpp::as<std::string>(fps_hex(i))));
      }
    } else {
      std::vector<FingerprintName> unsorted_names;
      unsorted_names.reserve(n);
      if (Rcpp::is<Rcpp::CharacterVector>(passed_names)) {
        auto passed_names_vec = Rcpp::as<Rcpp::CharacterVector>(passed_names);
        for (auto& x: passed_names_vec)
          unsorted_names.push_back(convert_name(Rcpp::as<std::string>(x)));
      } else if (Rcpp::is<Rcpp::IntegerVector>(passed_names) || Rcpp::is<Rcpp::NumericVector>(passed_names)) {
        unsorted_names = Rcpp::as< std::vector<FingerprintName> >(passed_names);
      } else {
        Rcpp::stop(
          "Fingerprint indices (names) must be passed as positive integers, numerics, or strings representing integers."
        );
      }
      // Make vector that contains indices that sort the name vector
      std::vector<size_t> sort_vector(n);
      std::iota(sort_vector.begin(), sort_vector.end(), 0);
      std::sort(
        sort_vector.begin(), sort_vector.end(),
        [&](const size_t i, const size_t j){return unsorted_names.at(i) < unsorted_names.at(j);}
      );
      for (auto &i: sort_vector) {
        fp_names.push_back(unsorted_names.at(i));
        fps.push_back(hex2fp(Rcpp::as<std::string>(fps_hex(i))));
      }
    }
  }

  // Constructor accepts a file path to load fingerprints from binary file
  MorganFPS(const std::string& filename, const bool from_file) {
    std::ifstream in_stream;
    in_stream.open(filename, std::ios::in | std::ios::binary);

    // std::in_stream in_stream(filename, std::ios::in | std::ios::binary);
    // std::vector<char> file_buffer(std::istreambuf_iterator<char>(input), {});
    // std::vector<char> decompressed_buffer;

    std::vector<char> decompressed_buffer;

    char magic[] = "xORGANFPS";
    in_stream.read(magic, 9);

    if (strcmp(magic, "MORGANFPS") != 0) {
      Rcpp::stop("File is incompatible, doesn't start with 'MORGANFPS': '%s'", magic);
    }

    FingerprintN n;
    in_stream.read(reinterpret_cast<char*>(&n), sizeof(FingerprintN));
    Rcpp::Rcerr << "Reading " << n << " fingerprints from file\n";
    fps.resize(n);
    fp_names.resize(n);

    size_t size_next_block;
    in_stream.read(reinterpret_cast<char*>(&size_next_block), sizeof(size_t));
    Rcpp::Rcerr << "Fingerprint block has " << size_next_block << " bytes\n";

    size_t bytes_decompressed = zstd_frame_decompress(
      in_stream, size_next_block, decompressed_buffer
    );

    if (fps.size() * sizeof(Fingerprint) != bytes_decompressed) {
      Rcpp::stop("Size of decompressed fingerprints differs from expected size.");
    };

    std::memcpy(reinterpret_cast<char*>(fps.data()), decompressed_buffer.data(), bytes_decompressed);
    Rcpp::Rcerr << "Fingerprints decompressed\n";

    in_stream.read(reinterpret_cast<char*>(&size_next_block), sizeof(size_t));
    Rcpp::Rcerr << "Names block has " << size_next_block << " bytes\n";

    bytes_decompressed = zstd_frame_decompress(
      in_stream, size_next_block, decompressed_buffer
    );

    if (fp_names.size() * sizeof(FingerprintName) != bytes_decompressed) {
      Rcpp::stop("Size of decompressed fingerprint names differs from expected size.");
    };

    std::memcpy(reinterpret_cast<char*>(fp_names.data()), decompressed_buffer.data(), bytes_decompressed);
    Rcpp::Rcerr << "Names decompressed\n";
  }

  // Tanimoto similarity between drugs i and j
  double tanimoto(Rcpp::RObject &i, Rcpp::RObject &j) {
    return jaccard_fp(fp_index(i), fp_index(j));
  }

  // Tanimoto similarity of drug i to every other drug
  Rcpp::DataFrame tanimoto_all(Rcpp::RObject &x) {
    const Fingerprint fp_other = fp_index(x);
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
    const Fingerprint fp_other = hex2fp(other);
    Rcpp::NumericVector res(fps.size());
    for (int i = 0; i < fps.size(); i++) {
      res[i] = jaccard_fp(fps[i], fp_other);
    }
    return Rcpp::DataFrame::create(
      Rcpp::Named("id") = fp_names,
      Rcpp::Named("structural_similarity") = res
    );
  }

  void save_file(const std::string& filename) {
    save_file(filename, 3);
  }

  // Save binary fp file
  void save_file(const std::string& filename, const int& compression_level=3) {
    if (compression_level < 1 || compression_level > 22)
      Rcpp::stop("Compression level must be between 0 and 22. Default = 3");

    FingerprintN n = fps.size();
    Rcpp::Rcerr << "Wrinting " << n << " fingerprints\n";

    std::ofstream out_stream;
    out_stream.open(filename, std::ios::out | std::ios::binary | std::ios::trunc);
    size_t input_size;

    std::vector<char> out_buffer;

    // out_stream << "MORGANFPS";
    out_stream.write("MORGANFPS", 9);
    out_stream.write(reinterpret_cast<char*>(&n), sizeof(FingerprintN));

    input_size = fps.size() * sizeof(Fingerprint);
    out_buffer.resize(ZSTD_compressBound(input_size));
    const size_t fingerprints_compressed = ZSTD_compress(
      out_buffer.data(), out_buffer.size(),
      reinterpret_cast<char *>(fps.data()), input_size,
      compression_level
    );
    if (ZSTD_isError(fingerprints_compressed)) {
      Rcpp::stop("Error compressing fingerprints: %s", ZSTD_getErrorName(fingerprints_compressed));
    }

    Rcpp::Rcerr << "Fingerprints compressed " << fingerprints_compressed << " bytes\n";
    // Save number of bytes of the compressed data. Important for finding
    // second block with names for decompression
    out_stream.write(reinterpret_cast<const char*>(&fingerprints_compressed), sizeof(size_t));
    out_stream.write(out_buffer.data(), fingerprints_compressed);
    Rcpp::Rcerr << "Wrote fingerprints\n";

    input_size = fp_names.size() * sizeof(FingerprintName);
    out_buffer.resize(ZSTD_compressBound(input_size));
    const size_t names_compressed = ZSTD_compress(
      out_buffer.data(), out_buffer.size(),
      reinterpret_cast<char *>(fp_names.data()), input_size,
      compression_level
    );
    if (ZSTD_isError(names_compressed)) {
      Rcpp::stop("Error compressing fingerprint names: %s", ZSTD_getErrorName(names_compressed));
    }

    Rcpp::Rcerr << "Names compressed " << names_compressed << " bytes\n";
    // Save number of bytes of the compressed data. Important for finding
    // second block with names for decompression
    out_stream.write(reinterpret_cast<const char*>(&names_compressed), sizeof(size_t));
    out_stream.write(out_buffer.data(), names_compressed);
    Rcpp::Rcerr << "Wrote Names\n";

    out_stream.close();
  }

  // Size of the dataset in bytes
  int size() {
    return fps.size() * sizeof(Fingerprint);
  }

  // Number of elements
  size_t n() {
    return fps.size();
  }

  std::vector<Fingerprint> fps;
  std::vector<FingerprintName> fp_names;

private:

  Fingerprint& fp_index(Rcpp::RObject& x) {
    FingerprintName x_name = convert_name(x);
    auto fp_pt = std::lower_bound(fp_names.begin(), fp_names.end(), x_name);
    if (*fp_pt != x_name)
      Rcpp::stop("Fingerprint %u not found", x_name);
    return fps.at(fp_pt - fp_names.begin());
  }

};

// Expose all relevant classes through an Rcpp module
RCPP_EXPOSED_CLASS(MorganFPS)
RCPP_MODULE(morgan_cpp) {

  using namespace Rcpp;

  class_<MorganFPS>( "MorganFPS" )
    .constructor<Rcpp::CharacterVector>("Construct fingerprint collection from vector of fingerprints")
    .constructor<std::string, bool>("Construct fingerprint collection from binary file")
    .method("size", &MorganFPS::size, "Size of the data in bytes")
    .method("n", &MorganFPS::n, "Number of elements")
    .method("tanimoto", &MorganFPS::tanimoto,
	    "Similarity between two fingerprints in the collection")
    .method("tanimoto_all", &MorganFPS::tanimoto_all,
	    "Similarity of a fingerprint against all other fingerprints in the collection")
    .method("tanimoto_ext", &MorganFPS::tanimoto_ext,
	    "Similarity of an external fingerprints against all fingerprints in the collection")
    .method("save_file", (void (MorganFPS::*)(const std::string&, const int&)) (&MorganFPS::save_file),
	    "Save fingerprints to file in binary format. Compression level between 0 and 22 (default 3).")
    .method("save_file", (void (MorganFPS::*)(const std::string&)) (&MorganFPS::save_file),
	    "Save fingerprints to file in binary format Compression level between 0 and 22 (default 3).")
    .field_readonly("fingerprints", &MorganFPS::fps)
    .field_readonly("names", &MorganFPS::fp_names)
    ;
}

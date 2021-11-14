#include <Rcpp.h>
#include "zstd/zstd.h"
#include <array>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

#include "utils.hpp"

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

Fingerprint convert_fp(const Rcpp::CharacterVector& fps_hex) {
  if (fps_hex.length() != 1)
    Rcpp::stop("Requires exactly one fingerprint");
  std::string format = guess_fp_format(fps_hex);
  auto string_to_fp = select_fp_reader(format);
  return string_to_fp(Rcpp::as<std::string>(fps_hex(0)));
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
double tanimoto(const Rcpp::CharacterVector& s1, const Rcpp::CharacterVector& s2) {
  const Fingerprint& fp1 = convert_fp(s1);
  const Fingerprint& fp2 = convert_fp(s2);
  return jaccard_fp(fp1, fp2);
}

void convert_fps(
    const Rcpp::CharacterVector& fps_hex,
    std::vector<FingerprintName>& out_names,
    std::vector<Fingerprint>& out_fps
) {
  std::string format = guess_fp_format(fps_hex);
  auto string_to_fp = select_fp_reader(format);
  size_t n = fps_hex.length();
  Rcpp::RObject passed_names = fps_hex.names();
  out_fps.reserve(n);
  out_names.reserve(n);
  if(passed_names.isNULL()) {
    for (int i = 0; i < n; i++) {
      out_names.push_back((FingerprintName) i + 1);
      out_fps.push_back(string_to_fp(Rcpp::as<std::string>(fps_hex(i))));
    }
  } else {
    std::vector<FingerprintName> unsorted_names = convert_name_vec(passed_names);
    // Make vector that contains indices that sort the name vector
    std::vector<size_t> sort_vector = sort_indices(unsorted_names);
    for (auto i: sort_vector) {
      out_names.push_back(unsorted_names.at(i));
      out_fps.push_back(string_to_fp(Rcpp::as<std::string>(fps_hex(i))));
    }
    // Checking for duplicate names
    auto duplicate_pair = std::adjacent_find(out_names.begin(), out_names.end());
    if (duplicate_pair != out_names.end())
      Rcpp::stop("Duplicate names are not allowed");
  }
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
//' @field tanimoto_threshold (threshold) similarity of all NxN combinations of fingerprints above the set threshold
//' @field tanimoto_subset (i,j) similarity of a set of fingerprints against another set, or all other fingerprints in the collection when j is NULL
//' @field tanimoto_ext (s) similarity between external hexadecimal string s and all
//'    fingerprints in the collection
//' @field save_file (path, compression_level) Save fingerprints to file in binary format
//' @field size number of bytes used to store the fingerprints
//' @importFrom Rcpp cpp_object_initializer
//' @export
class MorganFPS {

public:

  // Constructor accepts a named character vector of hex strings
  // either in full hexadecimal format or in packed RDKIT format
  MorganFPS(const Rcpp::CharacterVector& fps_hex) {
    convert_fps(fps_hex, fp_names, fps);
  }

  // Constructor accepts a file path to load fingerprints from binary file
  MorganFPS(const std::string& filename, const bool from_file) {
    read_file(filename);
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

  // Tanimoto similarity of all NxN combinations of fingerprints
  Rcpp::DataFrame tanimoto_threshold(double threshold) {
    Rcpp::IntegerVector id_1;
    Rcpp::IntegerVector id_2;
    Rcpp::NumericVector sims;
    for (int i = 0; i < fps.size(); i++) {
      if (i % 10000 == 0)
        Rcpp::Rcout << i << " done" << std::endl;
      Rcpp::checkUserInterrupt();
      for (int j = i + 1; j < fps.size(); j++) {
        auto sim = jaccard_fp(fps[i], fps[j]);
        if (sim > threshold) {
          id_1.push_back(fp_names[i]);
          id_2.push_back(fp_names[j]);
          sims.push_back(sim);
        }
      }
    }
    return Rcpp::DataFrame::create(
      Rcpp::Named("id_1") = id_1,
      Rcpp::Named("id_2") = id_2,
      Rcpp::Named("structural_similarity") = sims
    );
  }

  // Tanimoto similarity of drug list vs the same or another drug list
  Rcpp::DataFrame tanimoto_subset(Rcpp::RObject& x, Rcpp::RObject& y) {
    auto x_names = convert_sort_name_vec(x);
    auto x_fps = fp_index(x_names);
    std::vector<FingerprintName> x_name;
    std::vector<FingerprintName> y_name;
    std::vector<double> similarity;
    size_t n_total;
    if (y.isNULL()) {
      n_total = x_names.size() * n();
      x_name.reserve(n_total);
      y_name.reserve(n_total);
      similarity.reserve(n_total);
      for (int i = 0; i < x_names.size(); i++) {
        for (int j = 0; j < n(); j++) {
          x_name.push_back(x_names.at(i));
          y_name.push_back(fp_names.at(j));
          similarity.push_back(jaccard_fp(x_fps.at(i), fps.at(j)));
        }
      }
    } else {
      auto y_names = convert_sort_name_vec(y);
      auto y_fps = fp_index(y_names);
      n_total = x_names.size() * y_names.size();
      x_name.reserve(n_total);
      y_name.reserve(n_total);
      similarity.reserve(n_total);
      for (int i = 0; i < x_names.size(); i++) {
        for (int j = 0; j < y_names.size(); j++) {
          x_name.push_back(x_names.at(i));
          y_name.push_back(y_names.at(j));
          similarity.push_back(jaccard_fp(x_fps.at(i), y_fps.at(j)));
        }
      }
    }
    return Rcpp::DataFrame::create(
      Rcpp::Named("id_1") = x_name,
      Rcpp::Named("id_2") = y_name,
      Rcpp::Named("structural_similarity") = similarity
    );
  }

  // Tanimoto similarity of an external drug to every other drug
  //   in the collection
  Rcpp::DataFrame tanimoto_ext(const Rcpp::CharacterVector& others) {
    std::vector<FingerprintName> other_names;
    std::vector<Fingerprint> other_fps;
    convert_fps(others, other_names, other_fps);
    size_t nn = other_fps.size() * n();
    Rcpp::IntegerVector id_1(nn);
    Rcpp::IntegerVector id_2(nn);
    Rcpp::NumericVector sim(nn);
    size_t idx = 0;
    for (size_t i = 0; i < n(); i++) {
      for (size_t j = 0; j < other_fps.size(); j++) {
        sim[idx] = jaccard_fp(fps[i], other_fps[j]);
        id_1[idx] = other_names[j];
        id_2[idx] = fp_names[i];
        idx++;
      }
    }
    return Rcpp::DataFrame::create(
      Rcpp::Named("id_1") = id_1,
      Rcpp::Named("id_2") = id_2,
      Rcpp::Named("structural_similarity") = sim
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
    Rcpp::Rcout << "Wrinting " << n << " fingerprints\n";

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

    Rcpp::Rcout << "Fingerprints compressed " << fingerprints_compressed << " bytes\n";
    // Save number of bytes of the compressed data. Important for finding
    // second block with names for decompression
    out_stream.write(reinterpret_cast<const char*>(&fingerprints_compressed), sizeof(size_t));
    out_stream.write(out_buffer.data(), fingerprints_compressed);
    Rcpp::Rcout << "Wrote fingerprints\n";

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

    Rcpp::Rcout << "Names compressed " << names_compressed << " bytes\n";
    // Save number of bytes of the compressed data. Important for finding
    // second block with names for decompression
    out_stream.write(reinterpret_cast<const char*>(&names_compressed), sizeof(size_t));
    out_stream.write(out_buffer.data(), names_compressed);
    Rcpp::Rcout << "Wrote Names\n";

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
      Rcpp::stop("Fingerprint %i not found", x_name);
    return fps.at(fp_pt - fp_names.begin());
  }

  std::vector<std::reference_wrapper<Fingerprint>> fp_index(std::vector<FingerprintName>& names) {
    std::vector<std::reference_wrapper<Fingerprint>> hits;
    hits.reserve(names.size());
    auto last_idx = fp_names.begin();
    std::vector<FingerprintName>::iterator cur_idx;
    for (auto x: names) {
      cur_idx = std::lower_bound(last_idx, fp_names.end(), x);
      if (*cur_idx != x)
        Rcpp::stop("Fingerprint %i not found", x);
      else
        hits.push_back(std::ref(fps.at(cur_idx - fp_names.begin())));
        last_idx = cur_idx;
    }
    return hits;
  }

  void read_file(std::string filename) {
    std::ifstream in_stream;
    in_stream.open(filename, std::ios::in | std::ios::binary);

    // std::in_stream in_stream(filename, std::ios::in | std::ios::binary);
    // std::vector<char> file_buffer(std::istreambuf_iterator<char>(input), {});
    // std::vector<char> decompressed_buffer;

    char magic[] = "xORGANFPS";
    in_stream.read(magic, 9);

    if (strcmp(magic, "MORGANFPS") != 0) {
      Rcpp::stop("File is incompatible, doesn't start with 'MORGANFPS': '%s'", magic);
    }

    FingerprintN n;
    in_stream.read(reinterpret_cast<char*>(&n), sizeof(FingerprintN));
    Rcpp::Rcout << "Reading " << n << " fingerprints from file\n";
    fps.resize(n);
    fp_names.resize(n);

    size_t size_next_block;
    in_stream.read(reinterpret_cast<char*>(&size_next_block), sizeof(size_t));
    Rcpp::Rcout << "Fingerprint block has " << size_next_block << " bytes\n";

    size_t expected_decompressed_size = fps.size() * sizeof(Fingerprint);

    zstd_frame_decompress(
      in_stream, size_next_block, reinterpret_cast<char*>(fps.data()),
      expected_decompressed_size
    );

    Rcpp::Rcout << "Fingerprints decompressed\n";

    in_stream.read(reinterpret_cast<char*>(&size_next_block), sizeof(size_t));
    Rcpp::Rcout << "Names block has " << size_next_block << " bytes\n";

    expected_decompressed_size = fp_names.size() * sizeof(FingerprintName);

    zstd_frame_decompress(
      in_stream, size_next_block, reinterpret_cast<char*>(fp_names.data()),
      expected_decompressed_size
    );

    Rcpp::Rcout << "Names decompressed\n";
  }

};

// https://stackoverflow.com/a/42585733/4603385
template <typename T0, typename T1>
bool typed_valid(SEXP* args, int nargs){
  return nargs == 2 && Rcpp::is<T0>(args[0]) && Rcpp::is<T1>(args[1]);
}

// Expose all relevant classes through an Rcpp module
RCPP_EXPOSED_CLASS(MorganFPS)
RCPP_MODULE(morgan_cpp) {

  using namespace Rcpp;

  class_<MorganFPS>( "MorganFPS" )
    .constructor<Rcpp::CharacterVector>("Construct fingerprint collection from vector of fingerprints")
    .constructor<std::string, bool>("Construct fingerprint collection from binary file", &typed_valid<std::string, bool>)
    .method("size", &MorganFPS::size, "Size of the data in bytes")
    .method("n", &MorganFPS::n, "Number of elements")
    .method("tanimoto", &MorganFPS::tanimoto,
	    "Similarity between two fingerprints in the collection")
    .method("tanimoto_all", &MorganFPS::tanimoto_all,
	    "Similarity of a fingerprint against all other fingerprints in the collection")
    .method("tanimoto_threshold", &MorganFPS::tanimoto_threshold,
	    "@field tanimoto_threshold (threshold) similarity of all NxN combinations of fingerprints above the set threshold")
    .method("tanimoto_subset", &MorganFPS::tanimoto_subset,
	    "Similarity of a set of fingerprints against another set or all other fingerprints in the collection.")
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

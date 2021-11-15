#include <Rcpp.h>
#include <array>
#include <vector>
#include <string>

#ifndef MORGANCPP_UTILS_H
#define MORGANCPP_UTILS_H


using Fingerprint = std::array<std::uint64_t, 32>;
using FingerprintName = std::int32_t;
using FingerprintN = std::uint64_t;

using namespace Rcpp;

struct FingerprintHasher
{
  size_t operator()(const Fingerprint& fp) const noexcept
  {
    size_t h = 0;
    for (auto x: fp)
      h ^= x << 1;
    return h;
  }
};

using FingerprintMap = std::unordered_map<Fingerprint, FingerprintName, FingerprintHasher>;

FingerprintName convert_name(std::string x);
FingerprintName convert_name(RObject& x);
std::vector<FingerprintName> convert_name_vec(RObject& names);
std::vector<size_t> sort_indices(std::vector<FingerprintName>& unsorted_names);
std::vector<FingerprintName> convert_sort_name_vec(RObject& names);
int parse_hex_char(const char& c);
Fingerprint raw2fp(const std::string& raw);
Fingerprint hex2fp(const std::string& hex);
Fingerprint rdkit2fp(const std::string& hex);
std::string guess_fp_format(const CharacterVector& fps_hex);
std::function<Fingerprint (const std::string&)> select_fp_reader(const std::string& format);
size_t zstd_frame_decompress(
    std::ifstream &in_stream, size_t &compressed_size, char* out_buffer,
    size_t &out_buffer_size
);

#endif

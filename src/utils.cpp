#include <Rcpp.h>
#include "zstd/zstd.h"
#include <array>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

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
// std::string hex2raw(const std::string& hex) {
//   if (hex.length() % 2 != 0) {
//     ::Rf_error("Hex input length must be multiple of 2");
//   }
//   // Convert hex to raw bytes
//   std::string raw(hex.length() / 2, '\0');
//   auto hi = std::cbegin(hex);
//   for (auto& ri : raw) {
//     ri = parse_hex_char(*hi++) | (parse_hex_char(*hi++) << 4);
//   }
//   return raw;
// }

std::string hex2raw(const std::string& hex_string) {
  if (hex_string.length() % 2 != 0) {
    ::Rf_error("Hex input length must be multiple of 2");
  }
  std::string byte_string;
  byte_string.reserve(hex_string.length() / 2);
  for (size_t i = 0; i < hex_string.length(); i += 2) {
    byte_string.push_back(
      static_cast<char>(
        (parse_hex_char(hex_string[i]) << 4) | parse_hex_char(hex_string[i + 1])
      )
    );
  }
  return byte_string;
}

// Convert ASCII hex string to fingerprint.
Fingerprint hex2fp(const std::string& hex) {
  if (hex.length() != 512) {
    ::Rf_error("Input hex string must be of length 512");
  }
  return raw2fp(hex2raw(hex));
}

// From https://github.com/rdkit/rdkit/blob/78aac3c1bcc8f652053fdab26e5fe835fdaea53b/Code/RDGeneral/StreamOps.h#L143
uint32_t readPackedIntFromStream(std::stringstream &ss) {
  uint32_t val, num;
  int shift, offset;
  char tmp;
  ss.read(&tmp, sizeof(tmp));
  if (ss.fail()) {
    throw std::runtime_error("failed to read from stream");
  }

  val = static_cast<unsigned char>(tmp);
  offset = 0;
  if ((val & 1) == 0) {
    shift = 1;
  } else if ((val & 3) == 1) {
    ss.read((char *)&tmp, sizeof(tmp));
    if (ss.fail()) {
      throw std::runtime_error("failed to read from stream");
    }

    val |= (static_cast<unsigned char>(tmp) << 8);
    shift = 2;
    offset = (1 << 7);
  } else if ((val & 7) == 3) {
    ss.read((char *)&tmp, sizeof(tmp));
    if (ss.fail()) {
      throw std::runtime_error("failed to read from stream");
    }

    val |= (static_cast<unsigned char>(tmp) << 8);
    ss.read((char *)&tmp, sizeof(tmp));
    if (ss.fail()) {
      throw std::runtime_error("failed to read from stream");
    }

    val |= (static_cast<unsigned char>(tmp) << 16);
    shift = 3;
    offset = (1 << 7) + (1 << 14);
  } else {
    ss.read((char *)&tmp, sizeof(tmp));
    if (ss.fail()) {
      throw std::runtime_error("failed to read from stream");
    }

    val |= (static_cast<unsigned char>(tmp) << 8);
    ss.read((char *)&tmp, sizeof(tmp));
    if (ss.fail()) {
      throw std::runtime_error("failed to read from stream");
    }

    val |= (static_cast<unsigned char>(tmp) << 16);
    ss.read((char *)&tmp, sizeof(tmp));
    if (ss.fail()) {
      throw std::runtime_error("failed to read from stream");
    }

    val |= (static_cast<unsigned char>(tmp) << 24);
    shift = 3;
    offset = (1 << 7) + (1 << 14) + (1 << 21);
  }
  num = (val >> shift) + offset;
  // num = EndianSwapBytes<LITTLE_ENDIAN_ORDER,HOST_ENDIAN_ORDER>(num);
  return num;
}

inline void fp_set_bit(Fingerprint& fp, size_t i) {
  fp[i / 32] |= 1 << i % 32;
}

// From https://github.com/rdkit/rdkit/blob/06027dcd05674787b61f27ba46ec0d42a6037540/Code/DataStructs/BitVect.cpp#L23
Fingerprint rdkit2fp(const std::string& hex) {
  auto raw = hex2raw(hex);
  std::stringstream ss(raw);

  Fingerprint fp;
  fp.fill(0);
  std::int32_t format = 0;
  std::uint32_t nOn = 0;
  std::int32_t size;
  std::int32_t version = 0;

  // earlier versions of the code did not have the version number encoded, so
  //  we'll use that to distinguish version 0
  ss.read(reinterpret_cast<char*>(&size), sizeof(size));
  if (size < 0) {
    version = -1 * size;
    // Rcpp::Rcout << "Size: " << size << " Version: " << version << std::endl;
    if (version == 16) {
      format = 1;
    } else if (version == 32) {
      format = 2;
    } else {
      throw std::runtime_error("bad version in BitVect pickle");
    }
    ss.read(reinterpret_cast<char*>(&size), sizeof(size));
  } else {
    throw std::runtime_error("invalid BitVect pickle");
  }

  ss.read(reinterpret_cast<char*>(&nOn), sizeof(nOn));
  // Rcpp::Rcout << "Bits set: " << nOn << std::endl;

  // if the either have older version or or version 16 with ints for on bits
  if ((format == 0) ||
      ((format == 1) && (size >= std::numeric_limits<unsigned short>::max()))) {
    std::uint32_t tmp;
    for (unsigned int i = 0; i < nOn; i++) {
      ss.read(reinterpret_cast<char*>(&tmp), sizeof(tmp));
      fp_set_bit(fp, tmp);
    }
  } else if (format == 1) {  // version 16 and on bits stored as short ints
    std::uint16_t tmp;
    for (unsigned int i = 0; i < nOn; i++) {
      ss.read(reinterpret_cast<char*>(&tmp), sizeof(tmp));
      fp_set_bit(fp, tmp);
    }
  } else if (format == 2) {  // run length encoded format
    std::uint32_t curr = 0;
    for (unsigned int i = 0; i < nOn; i++) {
      curr += readPackedIntFromStream(ss);
      // Rcpp::Rcout << "Offset: " << curr << std::endl;
      fp_set_bit(fp, curr);
      curr++;
    }
  }
  return fp;
}

std::string guess_fp_format(const Rcpp::CharacterVector& fps_hex) {
  std::string format;
  if (!fps_hex.hasAttribute("format")) {
    format = "full";
  } else {
    format = Rcpp::as<std::string>(fps_hex.attr("format"));
  }
  if (format != "full" && format != "rle")
    Rcpp::stop("Unknown format");
  return format;
}

std::function<Fingerprint (const std::string&)> select_fp_reader(const std::string& format) {
  if (format == "full") {
    return hex2fp;
  } else if (format == "rle") {
    return rdkit2fp;
  } else {
    Rcpp::stop("Unknown format");
  }
}


size_t zstd_frame_decompress(
    std::ifstream &in_stream, size_t &compressed_size, char* out_buffer,
    size_t &out_buffer_size
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
  if (decompressed_size != out_buffer_size) {
    Rcpp::stop("Decompressed size differs from output buffer size: %i and %i", decompressed_size, out_buffer_size);
  }
  size_t const decompressed_bytes = ZSTD_decompress(
    out_buffer, decompressed_size,
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

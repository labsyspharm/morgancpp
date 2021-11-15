#include <Rcpp.h>
#include <vector>
#include <array>
#include <unordered_map>

#include "utils.hpp"

using namespace Rcpp;

//' @name MorganMap
//' @title Morgan fingerprint collection for identity checking
//' @description Efficient structure for checking identity of Morgan fingerprints
//' @importFrom Rcpp cpp_object_initializer
//' @export
class MorganMap {

public:
  MorganMap(const CharacterVector& fps_hex) {
    auto n = fps_hex.length();
    auto format = guess_fp_format(fps_hex);
    auto string_to_fp = select_fp_reader(format);
    RObject passed_names = fps_hex.names();
    fps.reserve(n);
    if(passed_names.isNULL()) {
      for (int i = 0; i < n; i++) {
        fps.emplace(string_to_fp(as<std::string>(fps_hex[i])), i + 1);
      }
    } else {
      auto unsorted_names = convert_name_vec(passed_names);
      for (int i = 0; i < n; i++) {
        fps.emplace(string_to_fp(as<std::string>(fps_hex[i])), unsorted_names[i]);
      }
    }
  }

  DataFrame find_matches(const CharacterVector& fps_hex) {
    auto n = fps_hex.length();
    auto format = guess_fp_format(fps_hex);
    auto string_to_fp = select_fp_reader(format);
    RObject passed_names = fps_hex.names();
    std::vector<FingerprintName> id_1;
    std::vector<FingerprintName> id_2;
    if(passed_names.isNULL()) {
      for (int i = 0; i < n; i++) {
        auto search = fps.find(string_to_fp(as<std::string>(fps_hex[i])));
        if (search != fps.end()) {
          id_1.push_back(i);
          id_2.push_back(search->second);
        }
      }
    } else {
      auto unsorted_names = convert_name_vec(passed_names);
      for (int i = 0; i < n; i++) {
        auto search = fps.find(string_to_fp(as<std::string>(fps_hex[i])));
        if (search != fps.end()) {
          id_1.push_back(unsorted_names[i]);
          id_2.push_back(search->second);
        }
      }
    }
    return DataFrame::create(
      Named("id_1") = id_1,
      Named("id_2") = id_2
    );
  }

  FingerprintMap fps;

};

// Expose all relevant classes through an Rcpp module
RCPP_EXPOSED_CLASS(MorganMap)
RCPP_MODULE(morgan_identity_cpp) {

  using namespace Rcpp;

  class_<MorganMap>( "MorganMap" )
    .constructor<CharacterVector>("Construct fingerprint collection from vector of fingerprints")
    .method("find_matches", &MorganMap::find_matches,
         "Find identical fingerprints");
}

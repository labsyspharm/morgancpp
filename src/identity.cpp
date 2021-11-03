#include <Rcpp.h>
#include <vector>
#include <array>
#include <unordered_map>

#include "utils.hpp"

//' @name MorganMap
//' @title Morgan fingerprint collection for identity checking
//' @description Efficient structure for checking identity of Morgan fingerprints
//' @importFrom Rcpp cpp_object_initializer
//' @export
class MorganMap {

public:
  MorganMap(const Rcpp::CharacterVector& fps_hex) {
    auto n = fps_hex.length();
    Rcpp::RObject passed_names = fps_hex.names();
    fps.reserve(n);
    if(passed_names.isNULL()) {
      for (int i = 0; i < n; i++) {
        fps.emplace(hex2fp(Rcpp::as<std::string>(fps_hex[i])), i + 1);
      }
    } else {
      auto unsorted_names = convert_name_vec(passed_names);
      for (int i = 0; i < n; i++) {
        fps.emplace(hex2fp(Rcpp::as<std::string>(fps_hex[i])), unsorted_names[i]);
      }
    }
  }

  Rcpp::DataFrame find_matches(const Rcpp::CharacterVector& fps_hex) {
    auto n = fps_hex.length();
    Rcpp::RObject passed_names = fps_hex.names();
    std::vector<FingerprintName> id_1;
    std::vector<FingerprintName> id_2;
    if(passed_names.isNULL()) {
      for (int i = 0; i < n; i++) {
        auto search = fps.find(hex2fp(Rcpp::as<std::string>(fps_hex[i])));
        if (search != fps.end()) {
          id_1.push_back(i);
          id_2.push_back(search->second);
        }
      }
    } else {
      auto unsorted_names = convert_name_vec(passed_names);
      for (int i = 0; i < n; i++) {
        auto search = fps.find(hex2fp(Rcpp::as<std::string>(fps_hex[i])));
        if (search != fps.end()) {
          id_1.push_back(unsorted_names[i]);
          id_2.push_back(search->second);
        }
      }
    }
    return Rcpp::DataFrame::create(
      Rcpp::Named("id_1") = id_1,
      Rcpp::Named("id_2") = id_2
    );
  }

  FingerprintMap fps;

};

// Expose all relevant classes through an Rcpp module
RCPP_EXPOSED_CLASS(MorganMap)
RCPP_MODULE(morgan_identity_cpp) {

  using namespace Rcpp;

  class_<MorganMap>( "MorganMap" )
    .constructor<Rcpp::CharacterVector>("Construct fingerprint collection from vector of fingerprints")
    .method("find_matches", &MorganMap::find_matches,
         "Find identical fingerprints");
}

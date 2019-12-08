#include <Rcpp.h>
#include <bitset>

// Converts character hex to character binary
std::bitset<2048> hex2bin( const std::string& hex )
{
  // Verify that the string contains 2048 bits of data
  if( hex.length() != 512 )
    ::Rf_error("Input hex string must be of length 512");

  // Convert hex to binary
  // Perform the conversion by hand to keep things optimized
  std::string bin;
  for( unsigned int i = 0; i != hex.length(); ++i )
    {
      switch(toupper(hex[i]))
	{
        case '0': bin += "0000"; break;
        case '1': bin += "0001"; break;
        case '2': bin += "0010"; break;
        case '3': bin += "0011"; break;
        case '4': bin += "0100"; break;
        case '5': bin += "0101"; break;
        case '6': bin += "0110"; break;
        case '7': bin += "0111"; break;
        case '8': bin += "1000"; break;
        case '9': bin += "1001"; break;
        case 'A': bin += "1010"; break;
        case 'B': bin += "1011"; break;
        case 'C': bin += "1100"; break;
        case 'D': bin += "1101"; break;
        case 'E': bin += "1110"; break;
        case 'F': bin += "1111"; break;
	default: bin += "0000"; break;
	}
    }

  // Instantiate the bitset from the binary string
  return std::bitset<2048>(bin);
}

// Jaccard similarity of two binary strings of length 2048
double bin_jaccard( const std::bitset<2048>& b1, const std::bitset<2048>& b2 )
{
  int val1 = (b1 & b2).count();    // intersection
  int val2 = (b1 | b2).count();    // union
  return static_cast<double>(val1) / val2;
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
double tanimoto( const std::string& s1, const std::string& s2 )
{
  std::bitset<2048> b1 = hex2bin(s1);
  std::bitset<2048> b2 = hex2bin(s2);
  return bin_jaccard( b1, b2 );
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
  MorganFPS( const std::vector<std::string>& vhx )
  {
    for( unsigned int i = 0; i < vhx.size(); ++i )
      {
	std::bitset<2048> bs = hex2bin(vhx[i]);
	vbs.push_back(bs);
      }
  }

  // Tanimoto similarity between drugs i and j
  // Adjust for 0-based indexing
  double tanimoto( unsigned int i, unsigned int j )
  {
    if( i >= vbs.size() || j >= vbs.size() )
      ::Rf_error("Index out of range");
    return bin_jaccard( vbs[i-1], vbs[j-1] );
  }

  // Tanimoto similarity of drug i to every other drug
  std::vector<double> tanimoto_all( unsigned int i )
  {
    // Argument verification
    if( i >= vbs.size() )
      ::Rf_error("Index out of range");

    // Iterate over the collection
    std::vector<double> res( vbs.size() );
    for( unsigned int ii = 0; ii < vbs.size(); ++ii )
      res[ii] = bin_jaccard( vbs[i-1], vbs[ii] );
    return res;
  }

  // Tanimoto similarity of an external drug to every other drug
  //   in the collection
  std::vector<double> tanimoto_ext( const std::string& other )
  {
    std::bitset<2048> b = hex2bin(other);
    std::vector<double> res( vbs.size() );
    for( unsigned int ii = 0; ii < vbs.size(); ++ii )
      res[ii] = bin_jaccard( b, vbs[ii] );
    return res;
  }

  // Size of the dataset
  unsigned int size() {return vbs.size() * sizeof(std::bitset<2048>);}
  
private:
  std::vector< std::bitset<2048> > vbs;
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

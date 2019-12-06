#include <Rcpp.h>
#include <bitset>

// Converts character hex to character binary
std::string hex2bin( const std::string& hex )
{
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
  return bin;
}

// Jaccard similarity of two binary strings of length 2048
double bin_jaccard( const std::bitset<2048>& b1, const std::bitset<2048>& b2 )
{
  int val1 = (b1 & b2).count();    // intersection
  int val2 = (b1 | b2).count();    // union
  return static_cast<double>(val1) / val2;
}

// Jaccard similarity of two hexadecimal strings of length 2048/4
// [[Rcpp::export]]
double hex_jaccard( const std::string& s1, const std::string& s2 )
{
  std::bitset<2048> b1( hex2bin(s1) );
  std::bitset<2048> b2( hex2bin(s2) );
  return bin_jaccard( b1, b2 );
}

// Morgan fingerprint collection
class MorganFPS {
public:
  // Constructor
  MorganFPS( const std::vector<std::string>& vhx )
  {
    for( unsigned int i = 0; i < vhx.size(); ++i )
      {
	std::string s = vhx[i];
	if( s.length() != 512 )
	  ::Rf_error("Input hex strings must be of length 512");

	// Create a new bit set and add it to collection
	std::bitset<2048> bs( hex2bin(s) );
	vbs.push_back(bs);
      }
  }

  // Tanimoto similarity between drugs i and j
  double tanimoto( unsigned int i, unsigned int j ) {return bin_jaccard( vbs[i-1], vbs[j-1] );}

  // Tanimoto similarity of drug i to every other drug
  std::vector<double> tanimoto_all( unsigned int i )
  {
    std::vector<double> res( vbs.size() );
    for( unsigned int ii = 0; ii < vbs.size(); ++ii )
      res[ii] = bin_jaccard( vbs[i-1], vbs[ii] );
    return res;
  }
  
private:
  std::vector< std::bitset<2048> > vbs;
};

using namespace Rcpp;

// Expose all relevant classes through an Rcpp module
RCPP_EXPOSED_CLASS(MorganFPS)
RCPP_MODULE(morgan_cpp) {
  
  class_<MorganFPS>( "MorganFPS" )
    .constructor< std::vector<std::string> >()
    .method("tanimoto", &MorganFPS::tanimoto )
    .method("tanimoto_all", &MorganFPS::tanimoto_all )
    ;
}

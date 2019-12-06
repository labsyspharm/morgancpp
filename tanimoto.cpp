#include <Rcpp.h>

#include <bitset>

// Converts character hex to character binary
// [[Rcpp::export]]
std::string hex2bin( std::string hex )
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
// [[Rcpp::export]]
double bin_jaccard( std::string s1, std::string s2 )
{
  if( s1.size() != 2048 || s2.size() != 2048 )
    {
      Rcpp::Rcerr << "One or more inputs do not contain 2048 bits" << std::endl;
      return 0.0;
    }
  std::bitset<2048> b1(s1);
  std::bitset<2048> b2(s2);
  int val1 = (b1 & b2).count();    // intersection
  int val2 = (b1 | b2).count();    // union
  return static_cast<double>(val1) / val2;
}

// Jaccard similarity of two hexadecimal strings of length 2048/4
// [[Rcpp::export]]
double hex_jaccard( std::string s1, std::string s2 )
{
  return bin_jaccard( hex2bin(s1), hex2bin(s2) );
}

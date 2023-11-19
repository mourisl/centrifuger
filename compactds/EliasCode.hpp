#ifndef _MOURISL_COMPACTDS_ELIASCODE
#define _MOURISL_COMPACTDS_ELIASCODE

#include "Utils.hpp"

namespace compactds {
class EliasCode
{
public:
  EliasCode() {}
  ~EliasCode() {}
  
  // These function will output
  // These methods can only encode positive numbers.
  // The bits are also reversed so accessing them is easier.
  // Even though the input value is 32-bit, the encoded bits can be greater than 32-bit.
  static WORD Unary(int in, int &l) 
  {
    l = in ;
    return 1ull << (in - 1);
  }
  
  // Elias gamma
  static WORD Gamma(int in, int &l)
  {
    int i ;
    const int n = Utils::CountBits(in) ;
    WORD ret = Unary(n, l) ;
    // the rightmost bit of Unary(n) and the leftmost bit of in are both 1, so we only need to shift by once.
    for (i = n - 2 ; i >= 0 ; --i, ++l)
    {
      ret |= (((in>>i)&1ull) << l) ;
    }
    //printf("%s: %d => %d %d %d\n", __func__, in, ret, n, l);
    return ret ; 
  }

  // Elias delta
  static WORD Delta(int in, int &l)
  {
    int i ;
    int n = Utils::CountBits(in) ;
    WORD ret = Gamma(n, l);
    
    for (i = n - 2 ; i >= 0 ; --i, ++l)
      ret |= (((in>>i)&1) << l) ;
    return ret ; // the leftmost bit of in is implicitly 1.
  }
  
  // Read in one Gamma encoded word starting from W's ith bits 
  // return: the value; l - # of processed bits
  static int ReadOneGamma(WORD *W, size_t i, int &l)
  {
    size_t j, k ;
    // Determine the length
    for (j = i ; Utils::BitRead(W, j) == 0 ; ++j)
      ;
    l = j - i + 1 ;
    int ret = 1 ;
    for (k = j + 1 ; k < j + l ; ++k)
      ret = (ret << 1) | Utils::BitRead(W, k) ;
    l = k - i ;
    return ret ;
  }
  
  // TODO: implement this
  static int ReadOneDelta(WORD *W, size_t i, int &out)
  {
    return 0 ;
  }
} ;
}

#endif

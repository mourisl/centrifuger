#ifndef _MOURISL_COMPACTDS_VARIABLESIZEELEM_ARRAY_DIRECTACCESS
#define _MOURISL_COMPACTDS_VARIABLESIZEELEM_ARRAY_DIRECTACCESS

#include <vector>

#include "Utils.hpp"
#include "Encode.hpp"
#include "FixedSizeElemArray.hpp"
#include "VariableSizeElemArray.hpp"


/*
 * The class for the array where each element has variable size
 */
namespace compactds {
class VariableSizeElemArray_DirectAccess : public VariableSizeElemArray
{
private:
  WORD **M ; // mark whether this is the last piece  
  WORD **P ; // the piece of block size b
 
  int b ; // block size
  int levelCnt ; // the number of dimensions for M and P
public:
  VariableSizeElemArray() 
  {
  }

  ~VariableSizeElemArray()
  {
    Free() ;
  }

  void Free()
  {
  }
  
  // Create the variable size element array  
  // b - block size
  // in - input array
  // n - the length of input array
  // 
  void InitFromArray(int blockSize, const unsigned int *in, const size_t &n)
  {
    int totalL = 0 ; // total bit length
    int maxL = 0 ; 
    int i ;
    for (i = 0 ; i < n ; ++i)
    {
      int bcnt = Utils::CountBits(in[i]) ;
      totalL += bcnt ;
      if (bcnt > maxL)
        maxL = bcnt ;
    }

    if (b <= 0)
      b = log(n) / log(2) ; //TODO: check 
    levelCount = DIV_CEIL(maxL, b) ;
   
    M = (WORD *)malloc(sizeof(WORD *) * levelCount) ;
    P = (WORD *)malloc(sizeof(WORD *) * levelCount) ;

    for (i = 0 ; i < n ; ++i)
    {
          
    }
  }

  unsigned int Read(int i)
  {
    return 0 ;
  }
} ;
}

#endif

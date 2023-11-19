#ifndef _MOURISL_COMPACTDS_VARIABLESIZEELEM_ARRAY_SAMPLEDPOINTERS
#define _MOURISL_COMPACTDS_VARIABLESIZEELEM_ARRAY_SAMPLEDPOINTERS

#include <vector>

#include "Utils.hpp"
#include "EliasCode.hpp"
#include "FixedSizeElemArray.hpp"

#include "VariableSizeElemArray.hpp"

/*
 * The class for the array where each element has variable size
 * Implement with sampled pointers (Section 3.2.1)
 */
namespace compactds {
class VariableSizeElemArray_SampledPointers: public VariableSizeElemArray
{
private:
  WORD *M ; // the compressed data
  size_t *P ; // sampled pointer 
  int b ;
  
  int space ;
public:
  VariableSizeElemArray_SampledPointers() 
  {
    M = NULL ;
    P = NULL ;
    space = 0 ;
  }

  ~VariableSizeElemArray_SampledPointers() 
  {
    Free() ;
  }

  void Free()
  {
    if (M != NULL)
      free(M) ;
    if (P != NULL)
      free(P) ;
    M = NULL ;
    P = NULL ;
  }
  
  // Create the variable size element array  
  // b - block size. has to be > 1
  // in - input array
  // n - the length of input array
  void InitFromArray(int blockSize, const unsigned int *in, const size_t &n) 
  {
    size_t i ;
    size_t totalEncodeBits = 0 ;
    for (i = 0 ; i < n ; ++i)
    {
      int bcnt = Utils::CountBits(in[i] + 1) ; // need to shift 1 to allow 0. 
      totalEncodeBits += 2 * bcnt - 1 ; // we use gamma encoding because the number of bits for each number is less than 32 in general, which makes it more efficient than delta encoding 
    }
    
    b = blockSize ;
    
    if (b <= 1)
    {
      b = sizeof(WORD) * 8 ; // extra overhead 1 bit per element 
    }

    size_t blockCnt = DIV_CEIL(n, b) ;
    P = (size_t *)malloc(blockCnt * sizeof(size_t)) ;
    space = blockCnt * sizeof(size_t) ;
    
    M = Utils::MallocByBits(totalEncodeBits) ; // We don't store the highest bit so -n
    space += Utils::BitsToWordBytes(totalEncodeBits) ;

    // Encode the data
    size_t sumL = 0 ;
    for (i = 0 ; i < n ; ++i)
    {
      if (i % b == 0)
      {
        P[i/b] = sumL ;
      }

      int l ; 
      WORD x = EliasCode::Gamma(in[i] + 1, l) ;
      Utils::BitsWrite(M, sumL, sumL + l - 1, x) ;

      //int tmpl ;
      //printf("i=%d: in[i]=%d encode=%lld l=%d. sumL=%d. decode=%d\n", i, in[i], x, l, sumL,
      //    EliasCode::ReadOneGamma(M, sumL, tmpl));
      sumL += l ;
    }
  }

  unsigned int Read(size_t i) 
  {
    size_t pi = i / b ; // index in P
    size_t j = pi * b ; 
    size_t offset = P[pi] ;
    int ret = 1 ;
    int l ;
    for (; j <= i ; ++j)
    {
      ret = EliasCode::ReadOneGamma(M, offset, l) ;
      offset += l ;
    }

    return ret - 1;
  }

  int GetSpace()
  {
    return space + sizeof(*this) ;  
  }
} ;
}
#endif

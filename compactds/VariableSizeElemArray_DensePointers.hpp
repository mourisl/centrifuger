#ifndef _MOURISL_COMPACTDS_VARIABLESIZEELEM_ARRAY_DENSEPOINTERS
#define _MOURISL_COMPACTDS_VARIABLESIZEELEM_ARRAY_DENSEPOINTERS

#include <vector>

#include "Utils.hpp"
#include "FixedSizeElemArray.hpp"

#include "VariableSizeElemArray.hpp"

/*
 * The class for the array where each element has variable size
 * Implement with dense pointers for constant time access (Section 3.2.2)
 */
namespace compactds {
class VariableSizeElemArray_DensePointers: public VariableSizeElemArray
{
private:
  WORD *M ; // the compressed data
  size_t *P ; // sampled pointer 
  FixedSizeElemArray offsets ; // the offset within each block
  int b ;
  size_t n ;
  int lastPosInM ; // the last position used in M 
  int space ;
public:
  VariableSizeElemArray_DensePointers() 
  {
    M = NULL ;
    P = NULL ;
    space = 0 ;
  }

  ~VariableSizeElemArray_DensePointers() 
  {
    Free() ;
  }

  void Free()
  {
    if (M != NULL)
      free(M) ;
    if (P != NULL)
      free(P) ;
    offsets.Free() ;
    M = NULL ; 
    P = NULL ;
  }
  
  // Create the variable size element array  
  // b - block size. has to be > 1
  // in - input array (non-negative)
  // n - the length of input array
  void InitFromArray(int blockSize, const unsigned int *in, const size_t &n) 
  {
    size_t i, j ;
    int maxL = 0 ;  
    size_t totalL = 0 ;
    this->n = n ;
    for (i = 0 ; i < n ; ++i)
    {
      int bcnt = Utils::CountBits(in[i] + 1) ; // need to shift 1 to allow 0. 
      if (bcnt > maxL) 
        maxL = bcnt ;
      totalL += bcnt ;
    }
    
    b = blockSize ;
    
    if (b <= 1)
    {
      b = CEIL(sizeof(WORD) * 8 * log(2)) ; // TODO: automatic block size determination
    }

    size_t blockCnt = DIV_CEIL(n, b) ;
    P = (size_t *)malloc(blockCnt * sizeof(size_t)) ;
    offsets.Malloc( Utils::Log2Ceil( (b - 1) * (double)(maxL-1) ), n - blockCnt) ; 
    space = blockCnt * sizeof(size_t) + offsets.GetSpace() - sizeof(offsets);

    M = Utils::MallocByBits(totalL - n) ; // We don't store the highest bit so -n
    space += Utils::BitsToWordBytes(totalL - n) ;

    // Encode the data
    size_t sumL = 0 ;
    size_t withinOffset = 0 ;
    // i for indexing the input array, j for indexing the offsets array
    for (i = 0, j = 0 ; i < n ; ++i)
    {
      int bcnt = Utils::CountBits(in[i] + 1) ;
      if (i % b == 0)
      {
        P[i/b] = sumL ;
        withinOffset = 0 ;      
      }
      else // Only store the offet of the first element
      {
        offsets.Write(j, withinOffset) ;
        ++j ;
      }
      if (in[i] == 0)
        continue ;
     
      Utils::BitsWrite(M, sumL, sumL + bcnt - 1 - 1, (in[i] + 1) & MASK(bcnt - 1)) ;
      //if (j > 0) printf("i=%d j=%d: in[i]=%d bcnt-1=%d. sumL=%d withinOffset=%d. offsets[j]=%d |elem|=%d\n", i, j, in[i], bcnt - 1, sumL, withinOffset, offsets.Read(j - 1), offsets.GetElemLength());
      sumL += bcnt - 1 ;
      withinOffset += bcnt - 1 ;
    }
    lastPosInM = sumL - 1 ; 
  }

  unsigned int Read(size_t i) 
  {
    int pi = i / b ; // index in P
    int presidual = i % b ;
    int nextpi = (i + 1) / b ;
    size_t ms, me ; // start and end in M.
    
    ms = P[pi] ;
    if (presidual > 0)
      ms = P[pi] + offsets.Read(i - pi - 1) ; // -pi because each block skip one elemtn in offsets
    if (i + 1 < n)
    {
      if (pi == nextpi) 
        me = P[pi] + offsets.Read(i + 1 - pi - 1) - 1;
      else
        me = P[nextpi] - 1 ;
    }
    else
      me = lastPosInM ; 

    //printf("\ni=%d: ms=%d me=%d. pi=%d offset=%d\n", i, ms, me, pi, offsets.Read(i - pi)) ;
    if (ms > me || (i == 0 && me == (size_t)-1)) 
      return 0 ;
   return (Utils::BitsRead(M, ms, me) | (1<<(me - ms + 1))) - 1;
  }

  int GetSpace()
  {
    return space + sizeof(*this) ;  
  }
} ;
}

#endif

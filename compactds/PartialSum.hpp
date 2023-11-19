#ifndef _MOURISL_COMPACTDS_PARTIALSUM
#define _MOURISL_COMPACTDS_PARTIALSUM

#include "Utils.hpp"
#include "Bitvector_Sparse.hpp"

namespace compactds {
class PartialSum
{
private:
  Bitvector_Sparse _B ; // underlying sparse bit vector
  size_t _n ;
  uint64_t _totalSum ;
public:
  PartialSum()
  {
    _n = _totalSum = 0 ;
  }

  ~PartialSum()
  {
    Free() ;
  }
  
  int GetSpace() 
  {
    return _B.GetSpace() + sizeof(*this) ;
  }

  void Free()
  {
    _B.Free() ;
  }
  
  void SetSupportSearch(bool supportSearch) 
  {
    _B.SetSupportRank(supportSearch) ;
  }

  void SetSpeed(int speed)
  {
    _B.SetSpeed(speed) ;
  }
  
  void Init(const int *array, const size_t n)
  {
    size_t i ;
    uint64_t *psum ;
    psum = (uint64_t *)malloc(sizeof(*psum) * (n+1)) ; // We store an extra element for all the length in sum

    psum[0] = 0 ;
    for (i = 1 ; i < n + 1 ; ++i) 
      psum[i] = psum[i - 1] + array[i - 1] ;
    InitFromPartialSum(psum, n) ;
    free(psum) ;
  }
  
  void Init(const size_t *array, const size_t n)
  {
    size_t i ;
    uint64_t *psum ;
    psum = (uint64_t *)malloc(sizeof(*psum) * (n+1)) ; // We store an extra element for all the length in sum

    psum[0] = 0 ;
    for (i = 1 ; i < n + 1 ; ++i) 
      psum[i] = psum[i - 1] + array[i - 1] ;
    InitFromPartialSum(psum, n) ;
    free(psum) ;
  }

  // n is the number of elements
  // psum records the partial sum before the i-th element
  // the last element should be the total sum, so need to store psum[n]  
  void InitFromPartialSum(const uint64_t *psum, const size_t n)
  {
    this->_n = n ;
    this->_totalSum = psum[n] ;
    _B.InitFromOnes(psum, n + 1, _totalSum) ;
  }

  // Initalize where the numbers are marked on bit vector
  //    i.e., the partial sum is the index on the bit array
  //    It assumes the lowest bit of W[0] is 1, and the last 
  //    index corresponds to the total sum
  void InitFromBitvector(WORD *W, const size_t wsize)
  {
    _B.Init(W, wsize) ;
    _n = _B.GetOneCnt() - 1 ;
    _totalSum = _B.GetLastOneIdx() ;
  }
  
  // Get the partial sum for index i
  // sum_0^[i-1] A[j]
  // Another interpretation is the summation for the first i elements.
  uint64_t Sum(size_t i) const
  {
    if (i == 0)
      return 0 ;
    else if (i >= _n)
      return _totalSum ;
    else
      // The input to Select is 1-based
      return _B.Select(i + 1) ;    
  }

  // Return the max i that Sum(i) <= the value of v
  size_t Search(const uint64_t v) const
  {
    if (v >= _totalSum)
      return _n ;
    return _B.Rank(1, (size_t)v) - 1 ;
  }
  
  // Read the value of an element
  int AccessValue(size_t i) const
  {
    if (i >= _n)
      return -1 ;
    return (int)(Sum(i + 1) - Sum(i)) ; 
  }

  void Save(FILE *fp)
  {
    SAVE_VAR(fp, _n);
    SAVE_VAR(fp, _totalSum);
    _B.Save(fp) ;
  }

  void Load(FILE *fp)
  {
    Free() ;

    LOAD_VAR(fp, _n);
    LOAD_VAR(fp, _totalSum);
    _B.Load(fp) ;
  }
} ;
}

#endif

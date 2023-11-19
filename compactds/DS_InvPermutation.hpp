#ifndef _MOURISL_COMPACTDS_DS_INVPERMUTATION
#define _MOURISL_COMPACTDS_DS_INVPERMUTATION

#include "Utils.hpp"
#include "Bitvector_Sparse.hpp"
#include "Bitvector_Plain.hpp"

// The standalone data structure for inverse query on a plain permutation using the idea of short churt.
// Time complexity: O(t)
// Space complexity: O(n/t * logn)
// Based Chapter 5.1. Difference is that the book samples is one off than this implementation
// This could be also useful for encoding the inverse function of an 1-to-1 mapping

namespace compactds {
class DS_InvPermutation
{
private:
  size_t _t ; // step size
  size_t _space ;
  Bitvector_Plain _B ; // mark whether a position is sampled
  FixedSizeElemArray _S ; // sampled pointer with value Pi^{-t}[x]
  size_t _sampledCnt ; // |_S| 
public:
  DS_InvPermutation()
  {
    _space = 0 ;
    _t = 0 ;
  }

  ~DS_InvPermutation()
  {
    Free() ;
  }

  void Free()
  {
    _S.Free() ;
  }

  size_t GetSpace()
  {
    return _space + _B.GetSpace() - sizeof(_B) + sizeof(*this) ;
  }

  void SetSampleRate(size_t t)
  {
    _t = t ;
  }

  void Init(size_t *Pi, size_t n)
  {
    if (_t == 0)
      _t = Utils::Log2Ceil(n) ;
    WORD *B = Utils::MallocByBits(n) ; // the label sampled positions
    WORD *V = Utils::MallocByBits(n) ; // the bits mark the cycles
    
    size_t i, j, k ;
    
    _sampledCnt = 0 ;
    for (i = 0 ; i < n ; ++i)
    {
      if (Utils::BitRead(V, i))
        continue ;
      Utils::BitSet(V, i) ;
      j = Pi[i] ;
      k = 1 ;
      while (j != i)
      {
        Utils::BitSet(V, j) ;
        if (k % _t == 0)
        {
          Utils::BitSet(B, j) ;
          ++_sampledCnt ;
        }
        j = Pi[j] ;
        ++k ;
      }

      if (k > _t) // may exist dangling part, without this, the time could be 2*_t-1
      {
        Utils::BitSet(B, i) ;
        ++_sampledCnt ;
      }
    }

    _B.Init(B, n) ;
    _S.Malloc(Utils::Log2Ceil(n), _sampledCnt) ;
    for (i = 0 ; i < n ; ++i)
    {
      if (!Utils::BitRead(V, i))
        continue ;

      Utils::BitFlip(V, i) ;
      j = Pi[i] ;
      while (Utils::BitRead(V, j))
      {
        if (Utils::BitRead(B, j))
        {
          // Since B[j]==1, use inclusive==0 automatically subtract the rank value by 1
          _S.Write( _B.Rank(1, j, 0), i) ;
          i = j ;
        }
        Utils::BitFlip(V, j) ;
        j = Pi[j] ;
      }
      if (Utils::BitRead(B, j))
      {
        _S.Write( _B.Rank(1, j, 0), i) ;
      }
      i = j ;
    }

    free(B) ;
    free(V) ;
  }

  //@return: Pi^{-1}[i]
  size_t Query(size_t *Pi, size_t i)
  {
    size_t j = i ;
    bool jumped = false ;
    while (Pi[j] != i)
    {
      if (!jumped && _B.Access(j))
      {
        j = _S.Read(_B.Rank(1, j, 0)) ;
        jumped = true ;
      }
      else
        j = Pi[j] ;
    }
    return j ;
  }
} ;
}

#endif

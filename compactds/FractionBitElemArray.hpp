#ifndef _MOURISL_COMPACTDS_FRACTIONBITELEM_ARRAY
#define _MOURISL_COMPACTDS_FRACTIONBITELEM_ARRAY

#include <stdint.h>
#include <stdlib.h>

#include <vector>

#include "Utils.hpp"

/*
 * The class for the array where each element is in the range of [0..d-1] and _d is far from the power of 2
 * The idea is that each WORD is a d-ary number (Section 3.1)
 */

namespace compactds {
class FractionBitElemArray
{
private:
  WORD *_W ;
  const int _w ;
  size_t _size ;
  size_t _d ; // element is in the range of [0..d-1]
  size_t _n ;
  int _k ; // number of elements per word
public:
  FractionBitElemArray():_w(8 * sizeof(WORD)) 
  {
    _W = NULL ;
  }

  ~FractionBitElemArray() 
  {
    Free() ;
  }

  // Allocate the memory for _n elements, where each element is in the range of [0..d-1]
  void Malloc(size_t d, size_t n)
  {
    this->_n = n ;
    this->_d = d ;
    _k = (int)(_w / ((double)log((double)_d) / (double)log(2.0))) ;
    _size = DIV_CEIL(n, _k) ;
    _W = Utils::MallocByBits(_size * WORDBITS) ;
  }
  
  // in - input array
  // n - the length of input array
  void InitFromArray(size_t d, const unsigned int *in, const size_t &n)
  {
    size_t i ;
    if (d == 0)
    {
      // We determine the best fixed size
      d = 0 ;
      for (i = 0 ; i < n ; ++i)
      {
        if (in[i] > d)
          d = in[i] ;
      }
      ++d ;
    }

    Malloc(d, n) ;
    for (i = 0 ; i < _n ; ++i)
      Write(i, in[i]) ;
  }

  void Free()
  {
    if (_W != NULL)
      free(_W) ;
    _W = NULL ;
  }
  
  // Get the i-th element
  unsigned Read(size_t i) const 
  {
    return (_W[i/_k] / Utils::PowerInt(_d, i%_k)) % _d ;
  }

  void Write(size_t i, int x)
  {
    size_t j = i / _k ;
    size_t p = Utils::PowerInt(_d, i%_k) ;
    _W[j] = _W[j] - ((_W[j] / p) %_d) * p + x * p ;
  }

  size_t GetSpace() const
  {
    return sizeof(_W[0]) * _size + sizeof(*this) ; 
  }

  int GetElemRange() const
  {
    return _d ;
  }
  
  size_t GetSize() const
  {
    return _n ;
  }

  const WORD* GetData() const
  {
    return _W ;
  }

  void Resize(size_t newn)
  {
    _n = newn ;
    _size = Utils::BitsToWords(DIV_CEIL(_n, _k)) ;
    _W = (WORD *)realloc(_W, _size * sizeof(WORD)) ;
  }
} ;
}

#endif

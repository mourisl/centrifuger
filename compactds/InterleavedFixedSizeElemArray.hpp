#ifndef _MOURISL_COMPACTDS_INTERLEAVEDFIXEDSIZEELEM_ARRAY
#define _MOURISL_COMPACTDS_INTERLEAVEDFIXEDSIZEELEM_ARRAY

// The class handles two levels of arrays. 
// Also a class where the first level is 64bit.

#include "Utils.hpp"

namespace compactds {
class InterleavedFixedSizeElemArray
{
private:
  size_t _l0, _l1 ; // length of element 0 and 1
  size_t _n0 ;
  size_t _f1 ; // frequency of element 1 after each element 0
  size_t _size ; //memory size, in words 
  WORD *_W ;
public:
  InterleavedFixedSizeElemArray()
  {
    _W = NULL ;
    _size = 0 ;
    _n0 = _f1 = 0 ;
    _l0 = _l1 = 0 ;
  }

  ~InterleavedFixedSizeElemArray()
  {
    Free() ;
  }

  void Free()
  {
    if (_n0 > 0)
    {
      free(_W) ;
      _W = NULL ;
      _size = 0 ;
      _n0 = _f1 = 0 ;
      _l0 = _l1 = 0 ;
    }
  }
  
  size_t GetSpace() const
  {
    return sizeof(_W[0]) * _size + sizeof(*this) ; 
  }

  void Malloc(size_t l0, size_t n0, int l1, size_t f1) 
  {
    Free() ;

    _l0 = l0 ;
    _n0 = n0 ;
    _l1 = l1 ;
    _f1 = f1 ;
    _size = Utils::BitsToWords(l0 * n0 + l1 * n0 * f1) ;
    _W = (WORD *)malloc(_size * sizeof(WORD)) ;
  }

  void Resize(size_t newn1)
  {
    _n0 = newn1 ;
    _size = Utils::BitsToWords(_l0 * _n0 + _l1 * _n0 * _f1) ;
    _W = (WORD *)realloc(_W, _size * sizeof(WORD)) ;
  }

  int GetElem0Length() const
  {
    return _l0 ;
  }

  int GetElem1Length() const
  {
    return _l1 ;
  }

  size_t GetSize0() const
  {
    return _n0 ;
  }

  size_t GetSize1() const
  {
    return _n0 * _f1 ;
  }

  void SetSize(size_t n0)
  {
    _n0 = n0 ;
  }

  WORD Read(int type, size_t i) const
  {
    if (type == 0)
    {
      const size_t offset = i * (_l0 + _f1 * _l1) ;
      return Utils::BitsRead(_W, offset, offset + _l0 - 1) ;
    }
    else
    {
      const size_t offset = (i / _f1) * (_l0 + _f1 * _l1) + _l0 + _l1 * (i%_f1);
      return Utils::BitsRead(_W, offset, offset + _l1 - 1) ;
    }
  }

  void Write(size_t type, size_t i, int x)
  {
    if (type == 0)
    {
      const size_t offset = i * (_l0 + _f1 * _l1) ;
      Utils::BitsWrite(_W, offset, offset + _l0 - 1, x) ;
    }
    else
    {
      const size_t offset = (i / _f1) * (_l0 + _f1 * _l1) + _l0 + _l1 * (i%_f1);
      Utils::BitsWrite(_W, offset, offset + _l1 - 1, x) ;
    }
  }
} ;

// Optimized for level 0 is 64bit integer. 
//  The second level will be paded 
class Interleaved64FixedSizeElemArray
{
private:
  size_t _l1 ; // length of element 0 and 1
  size_t _n0 ;
  size_t _f1 ; // frequency of element 1 after each element 0
  size_t _size ; //memory size 
  WORD *_W ;
  size_t _b ; // block size for each element 0 and attached element 1, in words
public:
  Interleaved64FixedSizeElemArray()
  {
    _W = NULL ;
    _size = 0 ;
    _n0 = _f1 = 0 ;
    _l1 = 0 ;
  }

  ~Interleaved64FixedSizeElemArray()
  {
    Free() ;
  }

  void Free()
  {
    if (_n0 > 0)
    {
      free(_W) ;
      _W = NULL ;
      _size = 0 ;
      _n0 = _f1 = 0 ;
      _l1 = 0 ;
    }
  }
  
  size_t GetSpace() const
  {
    return sizeof(_W[0]) * _size + sizeof(*this) ; 
  }

  void Malloc(size_t n0, int l1, size_t f1) 
  {
    Free() ;

    _n0 = n0 ;
    _l1 = l1 ;
    _f1 = f1 ;
    _b = Utils::BitsToWords(WORDBITS + DIV_CEIL(l1 * f1, WORDBITS) * WORDBITS) ;
    _size = Utils::BitsToWords(_n0 * _b * WORDBITS) ;
    _W = (WORD *)malloc(_size * sizeof(WORD)) ;
  }

  void Resize(size_t newn1)
  {
    _n0 = newn1 ;
    _size = Utils::BitsToWords(_n0 * _b * WORDBITS) ;
    _W = (WORD *)realloc(_W, _size * sizeof(WORD)) ;
  }

  int GetElemr0Length() const
  {
    return 64 ;
  }

  int GetElem2Length() const
  {
    return _l1 ;
  }

  size_t GetSize1() const
  {
    return _n0 ;
  }

  size_t GetSize2() const
  {
    return _n0 * _f1 ;
  }

  void SetSize(size_t n0)
  {
    _n0 = n0 ;
  }

  WORD Read0(size_t i) const
  {
    return _W[i * _b] ; 
  }

  WORD Read1(size_t i) const
  {
    const size_t tmp = i / _f1 ;
    const size_t offset = (tmp * _b + 1)* WORDBITS + (i - tmp * _f1) * _l1 ;
    return Utils::BitsRead(_W, offset, offset + _l1 - 1 ) ;
  }

  void Write0(size_t i, WORD x)
  {
    _W[i * _b] = x ;
  }

  void Write1(size_t i, int x)
  {
    const size_t tmp = i / _f1 ;
    const size_t offset = (tmp * _b + 1)* WORDBITS + (i - tmp * _f1) * _l1 ;
    Utils::BitsWrite(_W, offset, offset + _l1 - 1, x ) ;
  }
} ;


typedef InterleavedFixedSizeElemArray ILArray ;
typedef Interleaved64FixedSizeElemArray IL64Array ;
}

#endif 

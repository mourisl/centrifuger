#ifndef _MOURISL_COMPACTDS_FIXEDSIZEELEM_ARRAY
#define _MOURISL_COMPACTDS_FIXEDSIZEELEM_ARRAY

#include <stdint.h>
#include <stdlib.h>

#include <vector>

#include "Utils.hpp"

/*
 * The class for the array where each element is of fixed size
 * We use a word size w = 64bit to maximize the chance of within word access
 * Externally the index is continuous, but interally they are segmented by the word as the right of Fig 3.3 
 */

namespace compactds {
class FixedSizeElemArray
{
private:
  WORD *_W ;
  size_t _size ; // memory size in word 
  int _l ;
  size_t _n ;
public:
  FixedSizeElemArray() 
  {
    _W = NULL ;
    _size = 0 ;
    _n = 0 ;
    _l = 0 ;
  }

  ~FixedSizeElemArray() 
  {
    Free() ;
  }

  // Allocate the memory for _n elements, where each elements takes _l bits 
  void Malloc(int l, size_t n)
  {
    Free() ;
    this->_n = n ;
    this->_l = l ;
    _size = Utils::BitsToWords(l * n) ;
    _W = Utils::MallocByBits(l * n) ;
  }

  // _l - number of bits for each element. <=0: automatically decide
  // in - input array
  // n - the length of input array
  void InitFromArray(int l, const unsigned int *in, const size_t &n)
  {
    size_t i ;
    if (l <= 0)
    {
      // We determine the best fixed size
      l = 1 ;
      for (i = 0 ; i < n ; ++i)
      {
        int bitCounts = Utils::CountBits(in[i]) ;
        if (bitCounts > l)
          l = bitCounts ;
      }
    }

    Malloc(l, n) ;
    for (i = 0 ; i < n ; ++i)
      Write(i, in[i]) ;
  }
  
  void InitFromArray(int l, const size_t *in, const size_t &n)
  {
    size_t i ;
    if (l <= 0)
    {
      // We determine the best fixed size
      l = 1 ;
      for (i = 0 ; i < n ; ++i)
      {
        int bitCounts = Utils::CountBits(in[i]) ;
        if (bitCounts > l)
          l = bitCounts ;
      }
    }

    Malloc(l, n) ;
    for (i = 0 ; i < n ; ++i)
      Write(i, in[i]) ;
  }

  void Free()
  {
    if (_W != NULL)
      free(_W) ;
    _W = NULL ;
    _n = _size = 0 ;
    _l = 0 ;
  }
  
  // Get the i-th element
  uint64_t Read(size_t i) const 
  {
    return Utils::BitsRead(_W, i * _l, (i + 1)* _l - 1) ;
  }

  uint64_t operator[](size_t i) const
  {
    return Read(i) ;
  }

  void Write(size_t i, int x)
  {
    Utils::BitsWrite(_W, i * _l, (i + 1) * _l - 1, x) ;
  }

  /*uint64_t Read64(size_t i) const
  {
    return Utils::BitsRead(_W, i * _l, (i + 1)* _l - 1) ;
  }*/

  void Write64(size_t i, uint64_t x)
  {
    Utils::BitsWrite(_W, i * _l, (i + 1) * _l - 1, x) ;
  }

  size_t GetSpace() const
  {
    return sizeof(_W[0]) * _size + sizeof(*this) ; 
  }

  int GetElemLength() const
  {
    return _l ;
  }

  void SetElemLength(int l)
  {
    _l = l ;
  }
  
  size_t GetSize() const
  {
    return _n ;
  }

  // Assume we don't need to change the memory size
  void SetSize(size_t n) 
  {
    _n = n ;
  }

  // Set the elements from [start, start+len -1] to 0
  void SetZero(size_t start, size_t len )
  {
    size_t s = start * _l ;
    size_t e = (start + len) * _l - 1 ; 
    //if (start + len >= _n && e % WORDBITS > 0)
    //  e += (WORDBITS - e % WORDBITS - 1) ;

    size_t sk = s / WORDBITS ;
    size_t ek = e / WORDBITS ;
    if (sk < ek)
    {
      _W[sk] &= (MASK(s % WORDBITS)) ;
      _W[ek] &= (~MASK_WCHECK(e % WORDBITS + 1)) ;
      if (sk + 1 < ek)
        memset(_W + sk + 1, 0, (ek - sk - 1) * WORDBYTES) ;
    }
    else
      _W[sk] &= (MASK(s % WORDBITS) | (~MASK_WCHECK(e % WORDBITS + 1))) ;
  }

  const WORD* GetData() const
  {
    return _W ;
  }
  
  // The i-th element starts from the ret-th (0-based) bit in a word
  int GetElemOffsetInWord(size_t i) const
  {
    return (i * _l) % WORDBITS ;
  }

  // The i-th element is in the ret-th word (0-based)
  size_t GetElemWordIndex(size_t i) const
  {
    return (i * _l) / WORDBITS ;
  }
  
  // Return num elements starting from i.
  // @return: bit packed _W[i].._W[i + num - 1]
  WORD PackRead(size_t i, size_t num) const
  {
    return Utils::BitsRead(_W, i * _l, (i + num) * _l - 1) ;
  }

  WORD PackReadRev(size_t i, size_t num) const
  {
    size_t j ;
    WORD ret = 0 ;
    for (j = 0 ; j < num ; ++j)
      ret = (ret << _l) + Read(i + j) ;
    return ret ;
  }

  // Write WORD w contains num elements to the [i, i+num-1]
  void PackWrite(size_t i, WORD x, size_t num)
  {
    Utils::BitsWrite(_W, i * _l, (i + num) * _l - 1, x) ;
  }

  // Find the length of the matching prefix between A[s..e] and B[s..e]
  //  assumes _l is the same
  // If all match, return min(e-s+1, eb-sb+1)
  size_t PrefixMatchLen(size_t s, size_t e, const FixedSizeElemArray &B, size_t sb, size_t eb) const
  {
    if (e >= _n)
      e = _n - 1 ;
    if (eb >= B._n)
      eb = B._n - 1 ;
    size_t ai ;
    size_t bi ;

    int block = WORDBITS / _l ;
    ai = s ;
    bi = sb ;
    int len = MIN(e-s+1, eb-sb+1) ;
    if (len < block)
      block = len ;
    if (block > 1)
    {
      for ( ; ai + block - 1 <= e && bi + block - 1 <= eb ; 
          ai += block, bi += block)
      {
        WORD wa = PackRead(ai, block) ;
        WORD wb = B.PackRead(bi, block) ;
        if (wa == wb)
          continue ;

        return ai + Utils::CountTrailingZeros(wa^wb) / _l - s ;
      }
      /*else // When no element crosses the WORD boundary for A
      {
        size_t remainder = ai % block ;
        // To the WORD boundary
        if (remainder > 0)  
        {
          WORD wa = PackRead(ai, block - remainder) ;
          WORD wb = PackRead(bi, block - remainder) ;

          if (wa != wb)
            return ai + Utils::CountTrailingZeros(wa^wb) / _l - s ;
          else
          {
            ai += (block - remainder) ;
            bi += (block - remainder) ;
          }
        }

        for ( ; ai + block - 1 <= e && bi + block - 1 <= eb ; 
            ai += block, bi += block)
        {
          WORD wa = _W[ai/block] ;
          WORD wb = B.PackRead(bi, block) ;
          if (wa == wb)
            continue ;
          return ai + Utils::CountTrailingZeros(wa^wb) / _l - s ;
        }
      }*/
    }

    for ( ; ai <= e && bi <= eb ; ++ai, ++bi)
    {
      WORD smalla = Read(ai) ;
      WORD smallb = B.Read(bi) ; 
      if (smalla != smallb)
        return ai - s ;
    }

    return MIN(e - s + 1, eb - sb + 1) ;
  }
  
  // Compare A[s..e] and B[sb..eb]
  // @return: sign(A-B)
  int SubrangeCompare(size_t s, size_t e, const FixedSizeElemArray &B, size_t sb, size_t eb) const
  {
    if (_l != B._l)
      return _l - B._l ;
    if (e >= _n)
      e = _n - 1 ;
    if (eb >= B._n)
      eb = B._n - 1 ;
    size_t matchCnt = PrefixMatchLen(s, e, B, sb, eb) ;
    
    if (matchCnt == MIN(e - s + 1, eb - sb + 1))
    {
      if (e - s + 1 == eb - sb + 1)
        return 0 ;
      else if (e - s + 1 < eb - sb + 1)
        return -1 ;
      else
        return 1 ;
    }
    else
    {
      WORD smalla = Read(s + matchCnt) ;
      WORD smallb = B.Read(sb + matchCnt) ; 
      
      if (smalla < smallb) 
        return -1 ;
      else // they have to be different at this point 
        return 1 ;
    }
  }
  
  // Malloc by copying the first p element of B
  void InitFromOtherPrefix(const FixedSizeElemArray &B, size_t p)
  {
    Malloc(B._l, p) ;
    size_t wordBytes = Utils::BitsToWordBytes(_n * _l) ;
    memcpy(_W, B._W, wordBytes) ;
  }

  void Resize(size_t newn)
  {
    _n = newn ;
    _size = Utils::BitsToWords(_l * newn) ;
    _W = (WORD *)realloc(_W, _size * sizeof(WORD)) ;
  }
  
  // Reserve the space for m elements without changing current element
  void Reserve(size_t m)
  {
    if (m <= _n || Utils::BitsToWords(_l * m) <= _size)
      return;

    _size = Utils::BitsToWords(_l * m) ;
    if (_W != NULL)
      _W = (WORD *)realloc(_W, _size * sizeof(WORD)) ;
    else
      _W = Utils::MallocByBits(_l * m) ;
  }

  // push back another element to the end of the array.
  // This function also handles expand the array 
  void PushBack(int x)
  {
    if (Utils::BitsToWords(_l * _n) >= _size)
      Reserve(2 * _n) ;
    Write(_n, x);
    ++_n ;
  }

  // Push the first m element from b 
  void PushBack(const FixedSizeElemArray &b, size_t m)
  {
    size_t i ;
    if (m > b._n)
      m = b._n ;

    Reserve((_n + m + 2) / 2 * 3 + 1) ;
    int block = WORDBITS / _l ;
    for (i = 0 ; i < m ; i += block)
    {
      int actBlock = block ;
      if (i + block - 1 >= b._n)
        actBlock = b._n - i ;
      WORD wb = b.PackRead(i, actBlock) ;
      PackWrite(i + _n, wb, actBlock) ;
    }
    _n += m ;
  }

  void Print(FILE *fp, bool withIndex = false, char sep = ' ') const
  {
    size_t i ; 
    for (i = 0 ; i < _n ; ++i)
    {
      if (withIndex)
        fprintf(fp, "%d:%d%c", (int)i, (int)Read(i), sep) ;
      else
        fprintf(fp, "%d%c", (int)Read(i), sep) ;
    }
    fprintf(fp, "\n") ;
  }

  void Save(FILE *fp) 
  {
    SAVE_VAR(fp, _size) ;
    SAVE_VAR(fp, _l) ;
    SAVE_VAR(fp, _n) ;
    fwrite(_W, sizeof(_W[0]), Utils::BitsToWords(_n * _l), fp) ;
  }

  void Load(FILE *fp)
  {
    Free() ;
    LOAD_VAR(fp, _size) ;
    LOAD_VAR(fp, _l) ;
    LOAD_VAR(fp, _n) ;
    _W = Utils::MallocByBits(WORDBITS * _size) ;
    fread(_W, sizeof(_W[0]), Utils::BitsToWords(_n * _l), fp) ;
  }
} ;
}

#endif

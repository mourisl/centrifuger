#ifndef _MOURISL_COMPACTDS_SEQUENCE_PLAIN
#define _MOURISL_COMPACTDS_SEQUENCE_PLAIN

#include "Utils.hpp"
#include "Alphabet.hpp"
#include "Sequence.hpp"

#include "Bitvector_Plain.hpp"
#include "Bitvector_RunLength.hpp"

// The sequence representation where each alphabet is a bitvector 
namespace compactds {
template <class BvClass>
class Sequence_Plain: public Sequence
{
private:
  BvClass *_Bvs ; // bitvectors
  int _selectSpeed ;
public:
  Sequence_Plain() 
  {
    _selectSpeed = BITVECTOR_DEFAULT_SELECT_SPEED ; 
    _space = 0;
  }

  ~Sequence_Plain() 
  {
    Free() ;
  }

  void Free()
  {
    delete[] _Bvs ;
    Sequence::Free() ;
  }

  size_t GetSpace() 
  {
    return _space + _alphabets.GetSpace() - sizeof(_alphabets) + sizeof(*this) ;
  }

  void Init(const FixedSizeElemArray &S, size_t sequenceLength, const ALPHABET *alphabetMap)
  {
    size_t i, j, k ;
    _space = 0 ;
    
    if (_alphabets.GetSize() == 0)
      _alphabets.InitFromList(alphabetMap, strlen(alphabetMap)) ;
    
    this->_n = sequenceLength ;
    size_t alphabetSize = _alphabets.GetSize() ;
    
    _Bvs = new BvClass[alphabetSize] ;
    WORD *B = Utils::MallocByBits(_n) ;
    for (i = 0 ; i < alphabetSize ; ++i)
    {
      for (j = 0 ; j < _n ; j += WORDBITS)
      {
        WORD w = 0 ;
        for (k = 0 ; k < WORDBITS && j + k < _n ; ++k)
        {
          if (S.Read(j + k) == i)
            w |= (1ull<<k) ;
        }
        B[j / WORDBITS] = w ;
      }
      _Bvs[i].SetSelectSpeed( _selectSpeed ) ;
      _Bvs[i].Init(B, _n) ;
      _space += _Bvs[i].GetSpace() ;
    }
    free(B) ;
  }

  ALPHABET Access(size_t i) const 
  {
    size_t alphabetSize = _alphabets.GetSize() ;
    size_t j ;
    for (j = 0 ; j < alphabetSize ; ++j)
      if (_Bvs[j].Access(i) == 1)
        return _alphabets.Decode(j, 0) ;
    return 0 ;
  }
  
  size_t Rank(ALPHABET c, size_t i, int inclusive = 1) const 
  {
    int l ;
    return _Bvs[_alphabets.Encode(c, l)].Rank(1, i, inclusive) ; 
  }

  size_t Select(ALPHABET c, size_t i) const 
  {
    int l ;
    return _Bvs[_alphabets.Encode(c, l)].Select(i) ; 
  }

  void PrintStats()
  {
  }
} ;
}

#endif

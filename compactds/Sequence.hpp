#ifndef _MOURISL_COMPACTDS_SEQUENCE
#define _MOURISL_COMPACTDS_SEQUENCE

#include "Utils.hpp"
#include "Alphabet.hpp"
#include "FixedSizeElemArray.hpp"

namespace compactds {
class Sequence
{
protected:
  size_t _space ;
  Alphabet _alphabets ;
  size_t _n ; // sequence length
public:
  Sequence() {_space = 0 ; _n = 0 ;}
  ~Sequence() {}

  void SetAlphabet(const Alphabet &a)
  {
    _alphabets = a ;  
  }

  virtual void Save(FILE *fp)
  {
    SAVE_VAR(fp, _space) ;
    SAVE_VAR(fp, _n) ;
    _alphabets.Save(fp) ;
  }

  virtual void Load(FILE *fp)
  {
    LOAD_VAR(fp, _space) ;
    LOAD_VAR(fp, _n) ;
    _alphabets.Load(fp) ;
  }
  
  // Some flexible set function allow each sequence class to set its own special parametes
  //    such as the block size for runblock.
  virtual void SetExtraParameter(void *p) 
  {
  }
 
  virtual void Init(const FixedSizeElemArray &S, size_t sequenceLength, const ALPHABET *alphabetMap) = 0 ;
  virtual void Free() 
  {
    _n = 0 ;
    _space = 0 ;
    _alphabets.Free() ;
  }
  virtual size_t GetSpace() = 0 ; 
  virtual ALPHABET Access(size_t i) const = 0 ;
  virtual size_t Rank(ALPHABET c, size_t i, int inclusive = 1) const = 0 ;
  virtual size_t Select(ALPHABET c, size_t i) const = 0 ;
  virtual void PrintStats() = 0 ;
} ;
}

#endif

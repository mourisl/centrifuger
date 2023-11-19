#ifndef _MOURISL_COMPACTDS_SEQUENCE_PERMUTATION
#define _MOURISL_COMPACTDS_SEQUENCE_PERMUTATION

#include "Utils.hpp"
#include "Alphabet.hpp"
#include "FixedSizeElemArray.hpp"
#include "Sequence.hpp"

namespace compactds {
class Sequence_Permutation: public Sequence
{
private:
public:
  void Init(const FixedSizeElemArray &S, size_t sequenceLength, const ALPHABET *alphabetMap) 
  {
  }

  void Free()
  {
  }
  
  size_t GetSpace() 
  {
  }
  
  ALPHABET Access(size_t i) const 
  {
    return AccessLong(i) ;
  }

  size_t Rank(ALPHABET c, size_t i, int inclusive = 1) const
  {
    return RankLong(c, i, inclusive) ;
  }

  size_t Select(ALPHABET c, size_t i) const
  {
    return SelectLong(c, i) ;
  }

  size_t AccessLong(size_t i) const 
  {
  }
  
  size_t RankLong(size_t c, size_t i, int inclusive = 1) const
  {
  }

  size_t SelectLong(size_t c, size_t i) const
  {
  }


  void Save(FILE *fp)
  {
    Sequence::Save(fp) ;
  }

  void Load(FILE *fp)
  {
    Sequence::Load(fp) ;
  }
  
  void PrintStats() 
  {
  }
} ;
}

#endif

#ifndef _MOURISL_COMPACTDS_BITVECTOR
#define _MOURISL_COMPACTDS_BITVECTOR

#include "Utils.hpp"

#define BITVECTOR_DEFAULT_SELECT_SPEED 3

// The overall functionality of bitvector
namespace compactds {
class Bitvector
{
protected:
  size_t _space ;
public:
  Bitvector() {_space = 0 ;} 
  ~Bitvector() {}
  
  // W is the plain bit vector
  virtual void Init(const WORD *W, const size_t n) = 0 ; 
  virtual void Free() = 0 ;
  virtual size_t GetSpace() = 0;
  
  // Return the ith bits (0-based)
  virtual int Access(size_t i) const = 0 ;
  // Return the number of 1s before i
  virtual size_t Rank1(size_t i, int inclusive = 1) const = 0 ;
  // Return the index of th i-th (i is 1-based, so rank and select are inversible) 1
  // it is for 1 only for now
  virtual size_t Select(size_t i) const = 0 ;
 
  // Return the rightmost 1 in [0..i]
  // TODO: Handle the boundary cases
  size_t Pred(size_t i) const
  {
    return Select( Rank1(i) ) ;
  }

  // Return the leftmost 1 in [i..n-1]
  size_t Succ(size_t i) const
  {
    return Select( Rank1(i, /*inclusive=*/0) + 1 ) ;
  }

  // Return the number of 0s before i
  size_t Rank0(size_t i, int inclusive = 1) const
  {
    // There are i+1 elements in [0..i], and Rank(i) of them are 1's
    return i + inclusive - Rank1(i, inclusive) ;
  }

  size_t Rank(int type, size_t i, int inclusive = 1) const
  {
    if (type == 1)
      return Rank1(i, inclusive) ;
    else
      return Rank0(i, inclusive) ;
  }

  virtual void Save(FILE *fp) 
  {
    SAVE_VAR(fp, _space) ;
  }
  
  virtual void Load(FILE *fp) 
  {
    LOAD_VAR(fp, _space) ;
  }
} ;
}
#endif

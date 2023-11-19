#ifndef _MOURISL_COMPACTDS_BITVECTOR_RUNLENGTH
#define _MOURISL_COMPACTDS_BITVECTOR_RUNLENGTH

#include "Utils.hpp"

#include "Bitvector.hpp"
#include "Bitvector_Sparse.hpp"
#include "PartialSum.hpp"
#include "SimpleVector.hpp"

// The run-length bitvector built upon the sparse bit vector 
// Based on section: 
namespace compactds {
class Bitvector_RunLength: public Bitvector
{
protected:
  bool _zerofirst ; // whether the run-length array starts with 0 or _not.
  int _partialSumSpeed ;
  size_t _n ; // total _number of bits
  size_t _rcnt ; // the _number of runs

  PartialSum _R ; // the partial sum of runs
  PartialSum _O ; // the partial sum of 1s
public:
  Bitvector_RunLength() 
  {
    _partialSumSpeed = BITVECTOR_DEFAULT_SELECT_SPEED ;
  } 
  ~Bitvector_RunLength() {}
  
  void Free() 
  {
    _R.Free() ;
    _O.Free() ;
  }

  size_t GetSpace()
  {
    return _R.GetSpace() - sizeof(_R) 
      + _O.GetSpace() - sizeof(_O) + sizeof(*this) ;
  }

  void SetSelectSpeed(int speed)
  {
  }

  void SetPartialSumSpeed(int _partialSumSpeed)
  {
    _R.SetSpeed(_partialSumSpeed) ;
    _O.SetSpeed(_partialSumSpeed) ;
  }

  void SetSupportSelect(int supportSelect)
  {
    _O.SetSupportSearch(supportSelect) ;
  }
  
  // W is the plain bit vector
  void Init(const WORD *W, const size_t n)
  {
    size_t i ;
    if (n == 0)
      return ;
    _n = n ;
    _zerofirst =false ;
    if (Utils::BitRead(W, 0) == 0)
      _zerofirst = true ;
    
    WORD *B = Utils::MallocByBits(n + 1) ; // bits for the sums of runs
    WORD *oneB = Utils::MallocByBits(n + 1) ; // bits for the sums of runs of 1s
    Utils::BitSet(B, 0) ;
    Utils::BitSet(oneB, 0) ;
    size_t oneLen = 0 ; // run length for ones
    int prevc = 0 ;
    if (!_zerofirst)
    {
      prevc = 1 ;
      oneLen = 1 ;
    }

    for (i = 1 ; i < n ; ++i) 
    {
      int c = Utils::BitRead(W, i) ;
      if (c)
        ++oneLen ;
      if (c != prevc)
      {
        Utils::BitSet(B, i) ;
        if (c == 0) // previous c == 1
          Utils::BitSet(oneB, oneLen) ;
      }
      prevc = c ;
    }

    Utils::BitSet(B, n) ;
    Utils::BitSet(oneB, oneLen) ;
    _R.InitFromBitvector(B, n + 1) ;
    _O.InitFromBitvector(oneB, oneLen + 1) ;
    free(B) ;
    free(oneB) ;
    
    /*size_t len = 1 ;
    SimpleVector<int> rlens ;
    rlens.Reserve(_n / WORDBITS + 1) ;
    for (i = 1 ; i < _n ; ++i)
    {
      if (Utils::BitRead(W, i) != Utils::BitRead(W, i - 1))
      {
        rlens.PushBack(len) ; 
        len = 1 ;
      }
      else
        ++len ;
    }
    rlens.PushBack(len) ;
    InitFromRunLength(rlens.BeginAddress(), rlens.Size(), _n, _zerofirst) ;*/
  }

  void InitFromRunLength(const int *rlens, const size_t rcnt, const size_t n, const bool zerofirst)
  {
    this->_rcnt = rcnt ;
    this->_n = n ;
    this->_zerofirst = zerofirst ;
  
    _R.Init(rlens, _rcnt) ;

    size_t i = 0 ;
    uint64_t *oneSums = (uint64_t *)malloc(sizeof(*oneSums) * (_rcnt / 2 + 2));
  
    if (_zerofirst)
      i = 1 ;
    
    uint64_t psum = 0 ; 
    for ( ; i < _rcnt ; i += 2)
    {
      oneSums[i/2] = psum ;
      psum += rlens[i] ;
    }
    oneSums[i/2] = psum ;
    _O.InitFromPartialSum(oneSums, i/2) ;
    free(oneSums) ;
  }
    
  // Return the ith bits (0-based)
  int Access(size_t i) const
  {
    size_t ri = _R.Search(i) ;
    int inOne = (ri&1) ;  //whether ri is block for 1 or 0
    if (!_zerofirst)
      inOne = 1 - inOne ;
    return inOne ;
  }

  // Return the _number of 1s before i
  size_t Rank1(size_t i, int inclusive = 1) const 
  {
    size_t ri = _R.Search(i) ; // run-length block index 
    size_t oi = ri / 2 ;
    int inOne = (ri&1) ;  //whether ri is block for 1 or 0
    if (!_zerofirst)
      inOne = 1 - inOne ;
    if (!inOne)
    {
      // Each small block is a (00..11..) xx,
      //   or (11..00..) runs,
      //   so we _need to adjust whether we want the sum of 1's include the current small bock or _not.
      if (_zerofirst)
        return _O.Sum(oi) ;
      else
        return _O.Sum(oi + 1) ;
    }
    else
    {
      //The sum of 1s before current run and the _number of 1s in the current block
      return _O.Sum(oi) + (i - _R.Sum(ri) + inclusive) ;   
    }
  }
  
  // Return the index of th i-th (i is 1-based, so rank and select are inversible) 1
  // Did not have a select mode for 0 here.
  size_t Select(size_t i) const  
  {
    if (i == 0)
      return POSITIVE_INF ;

    --i ;
    size_t oi = _O.Search(i) ;
    // Map oi back to the ri 
    size_t ri = 2 * oi ;
    if (_zerofirst)
      ++ri ;
    return _R.Sum(ri) + (i - _O.Sum(oi)) ;
  }

  size_t Select(int type, size_t i) const
  {
    if (type == 1)
      return Select(i) ;
    else
      return POSITIVE_INF ; 
  }

  void Save(FILE *fp)
  {
    Bitvector::Save(fp) ;
    SAVE_VAR(fp, _zerofirst) ; 
    SAVE_VAR(fp, _partialSumSpeed) ; 
    SAVE_VAR(fp, _n) ; 
    SAVE_VAR(fp, _rcnt) ; 
    _R.Save(fp) ;
    _O.Save(fp) ;
  }

  void Load(FILE *fp)
  {
    Free() ;
    Bitvector::Load(fp) ;
    LOAD_VAR(fp, _zerofirst) ; 
    LOAD_VAR(fp, _partialSumSpeed) ; 
    LOAD_VAR(fp, _n) ; 
    LOAD_VAR(fp, _rcnt) ; 
    _R.Load(fp) ; 
    _O.Load(fp) ;
  }
} ;
}

#endif

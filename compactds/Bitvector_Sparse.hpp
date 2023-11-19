#ifndef _MOURISL_COMPACTDS_BITVECTOR_SPARSE
#define _MOURISL_COMPACTDS_BITVECTOR_SPARSE

#include "Utils.hpp"
#include "Bitvector.hpp"
#include "FixedSizeElemArray.hpp"
#include "Bitvector_Plain.hpp"

// The very sparse bitvector based on chaptor 4.4
// This data structure seems can be used for set predecessor query when 
// the elements are increasing?
//
// This is a super clever data structure. Come and think about it from time to time.
namespace compactds {
class Bitvector_Sparse: public Bitvector
{
private:
  size_t _n ; // total length of the bit vector
  size_t _onecnt ; // number of 1s in the bit vector
  size_t _lastOneIdx ; // the index of the last one
  int _lowerBits ; // the split for lower and upper bits

  FixedSizeElemArray _L ; // stores lower bits (of size r)
  Bitvector_Plain _H ; // stores the higher bits

  int _hSelectSpeed ; // this is the speed for H
public:
  Bitvector_Sparse() 
  {
    _lowerBits = 0 ;
    _hSelectSpeed = BITVECTOR_DEFAULT_SELECT_SPEED ;
  } 
  ~Bitvector_Sparse() 
  {
    Free() ;
  }
  
  void Free()
  {
    _n = 0 ;
    _L.Free() ;
    _H.Free() ;
  }

  size_t GetSpace()
  {
    return _space + _L.GetSpace() - sizeof(_L) + 
      _H.GetSpace() - sizeof(_H) + sizeof(*this) ;
  }

  void SetLowerBits(int lowerBits)
  {
    this->_lowerBits = lowerBits ;
  }

  // the speed is for H.select
  void SetSpeed(int speed)
  {
    _H.SetSelectSpeed(speed) ;
  }

  void SetSupportRank(bool supportRank)
  {
    _H.SetSelectTypeSupport(supportRank ? 3 : 2) ;
  }

  size_t GetLastOneIdx()
  {
    return _lastOneIdx ;
  }

  size_t GetOneCnt()
  {
    return _onecnt ;
  }

  // Init directly from the bit vector
  // W is the plain bit vector
  void Init(const WORD *W, const size_t n) 
  {
    size_t i, k ;
    size_t wordCnt = Utils::BitsToWords(n) ;

    _n = n ; 
    _onecnt = 0 ;
    for (i = 0 ; i < wordCnt ; ++i)
      _onecnt += Utils::Popcount(W[i]) ;
    
    if (_onecnt == 0 || wordCnt == 0)
    {
      _lastOneIdx = 0 ;
      return ;
    }
    
    // Get the last 1
    i = wordCnt - 1 ;
    while (1)
    {
      if (W[i] == 0)
      {
        --i ;
        continue ;
      }
      else
      {
        int j ;
        for (j = WORDBITS - 1 ; j >= 0 ; --j)
        {
          if ((W[i] >> j) & 1)
          {
            _lastOneIdx = i * WORDBITS + j ;
            break ;
          }
        }
        break ;
      }
      if (i == 0)
        break ;
      else
        --i ;
    }

    if (_lowerBits == 0)
      _lowerBits = int(log((double)n / _onecnt) / log(2.0)) ;

    if (_lowerBits < 1)
      _lowerBits = 1 ;
    _L.Malloc(_lowerBits, _onecnt) ;

    size_t hsize = (_lastOneIdx >> _lowerBits) + _onecnt + 1; // need +1 here to accommdate the max value
    ++hsize ; // Plus one here is because we want to append a 0 to the last block
    _H.Malloc(hsize) ;

    k = 0 ;
    for (i = 0 ; i < n ; i += WORDBITS)
    {
      WORD w = W[i/WORDBITS] ;
      if (w == 0) 
        continue ;
      size_t j ;
      for (j = 0 ; j < WORDBITS && i + j < n ; ++j)
        if ((w >>j) & 1) 
        {
          _L.Write(k, (i + j) & MASK(_lowerBits)) ;
          _H.BitSet(((i + j) >> _lowerBits) + k) ;
          ++k ;
        }
    }
    _H.Init() ;
  }
  
  // Init from the know positions of 1s
  void InitFromOnes(const uint64_t *S, const size_t onecnt, const size_t n)
  {
    size_t i ;

    _n = n ;
    _onecnt = onecnt ;
    if (_onecnt > 0)
      _lastOneIdx = S[onecnt - 1] ;
    else
    {
      _lastOneIdx = 0 ;
      return ;
    }

    if (_lowerBits == 0)
      _lowerBits = int(log((double)n / onecnt) / log(2.0)) ;
    
    if (_lowerBits < 1)
      _lowerBits = 1 ;
    _L.Malloc(_lowerBits, onecnt) ;

    size_t hsize = (S[onecnt - 1] >> _lowerBits) + onecnt + 1; // need +1 here to accommdate the max value
    ++hsize ; // Plus one here is because we want to append a 0 to the last block
    _H.Malloc(hsize) ; 
    for (i = 0 ; i < onecnt ; ++i)
    {
      _L.Write(i, S[i] & MASK(_lowerBits)) ;
      _H.BitSet((S[i] >> _lowerBits) + i) ;
    }
    _H.Init() ;
  }
  
  // Return the ith bits (0-based)
  int Access(size_t i) const
  {
    if (Pred(i) == i)
      return 1 ;
    else
      return 0 ;
  }

  // Return the number of 1s before i
  size_t Rank1(size_t i, int inclusive = 1) const 
  {
    if (inclusive == 0)
    {
      if (i == 0)
        return 0 ;
      else
        --i ; 
    }

    if (i >= _lastOneIdx) // this should contains the case that i>=n
      return _onecnt ;
    
    size_t iH = i >> _lowerBits ;
    size_t iL = i & MASK(_lowerBits) ;
    size_t l, m, r ;
    
    // We don't want to +1 for iH in select because the
    // 0 marks the beginning the block with starts with iH<<r
    //
    // The difference between Select(0, iH) and iH is the number 
    // of 1s in the range in the range of [0, iH<<r)
    // (This is because each observed 1 will shift Select right by 1)
    // Therefore, this gives the range in L, because L did not miss anything
    //   and has the number 
    size_t selectIH = 0 ;
    if (iH == 0)
      l = 0 ;
    else
    {
      selectIH = _H.Select(0, iH) ;
      l = selectIH - (iH - 1) ; // Note that the input to Select is 1-based, and the shift is 0-based, so we add (iH - 1) here.
    }
    
    // When i >= last one index, H.Select(iH)==H.Select(iH+1),
    //   then l>r.
    //   Fortunately, we handle this case at the beginning.
    if (iH == 0 || _H.Access( selectIH + 1 ) != 0)
      r = _H.Select(0, iH + 1) - iH ;
    else
      r = selectIH + 1 - iH ;

    if (l == r || _L.Read(l) > iL)
    {
      // The current r block is empty
      // or the first element in the block is greater than what we search for.
      // So the number of 1s before current block (l) is the answer
      return l ; 
    }

    // r points to the start of the next r block, so we need -1
    //   to make it match with the end of the current block
    // The l==r test above makes sure r-1 is non-negative here.
    --r ;
    while (l <= r)
    {
      m = (l + r) / 2 ;
      if (_L.Read(m) <= iL)
        l = m + 1 ;
      else
      {
        if (r == 0) 
          break ; // the test before the binary search make sure at least one element
                  // in the block is less than the desired target,
                  // so we can directly termiante the binary search.
        else
          r = m - 1 ;
      }
    }
    return l ; // l-1 is the last element index <= the desired one, so l is the number element   
  }

  // Return the index of th i-th (i is 1-based, so rank and select are inversible) 1
  size_t Select(size_t i) const
  {
    if (i > _onecnt)
      return _lastOneIdx ;
    if (i == 0)
      return POSITIVE_INF ;
    // Use (i-1) instead of i to convert the 1-based to 0-based, which is
    //  the base when creating H.
    return ((_H.Select(1, i) - (i - 1)) << _lowerBits) + _L.Read(i - 1) ;
  }

  size_t Select(int type, size_t i) const
  {
    if (type == 1)
      return Select(i) ;
    else
    {
      // Sadly, we can only do plain binary search 
      // Haven't tested it yet.
      // Don't recommend this operation.
      size_t l = 0 ;
      size_t r = _n - 1 ;
      size_t m ;

      while (l <= r)
      {
        m = (l + r) / 2 ;
        if (m - Select(m) < i)
          l = m + 1 ;
        else
        {
          if (m == 0)
            return 0 ;
          else
            r = m - 1 ;
        }
      }
      return r + 1 ;
    }
  }

  void Save(FILE *fp)
  {
    Bitvector::Save(fp) ;
    SAVE_VAR(fp, _n) ;
    SAVE_VAR(fp, _onecnt) ;
    SAVE_VAR(fp, _lastOneIdx) ;
    SAVE_VAR(fp, _lowerBits) ;
    SAVE_VAR(fp, _hSelectSpeed) ;
    _L.Save(fp) ;
    _H.Save(fp) ;
  }
  
  void Load(FILE *fp)
  {
    Free() ;
    Bitvector::Load(fp) ;
    LOAD_VAR(fp, _n) ;
    LOAD_VAR(fp, _onecnt) ;
    LOAD_VAR(fp, _lastOneIdx) ;
    LOAD_VAR(fp, _lowerBits) ;
    LOAD_VAR(fp, _hSelectSpeed) ;
    _L.Load(fp) ;
    _H.Load(fp) ;
  }
} ;
}

#endif

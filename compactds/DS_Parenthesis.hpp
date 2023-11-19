#ifndef _MOURISL_COMPACTDS_PARENTHESIS
#define _MOURISL_COMPACTDS_PARENTHESIS

#include "Utils.hpp"
#include "DS_RangeMinMaxTree.hpp"
#include "DS_PatternRankSelect.hpp"

namespace compactds {
class DS_Parenthesis
{
private:
  DS_RangeMinMaxTree _rmmTree ;
  DS_PatternRankSelect _patRS ;
  
  void GenerateRandomBalanceParenthesisSegment(WORD *B, size_t n, size_t i, size_t j)
  {
    if (j == i + 1)
    {
      Utils::BitsWrite(B, i, j, 2) ; // write binary 10 
      return ;
    }
    else if (j <= i)
    {
      return ;
    }

    size_t split = i + rand() % (j - i + 1) ;
    while ((split - i + 1) % 2 == 1 )
      split = i + rand() % (j - i + 1) ;
    
    Utils::BitSet(B, i) ;
    Utils::BitClear(B, j) ;
    if (split == i || split == j)
      GenerateRandomBalanceParenthesisSegment(B, n, i + 1, j - 1) ;
    else
    {
      Utils::BitClear(B, split) ;
      Utils::BitSet(B, split + 1) ;
      GenerateRandomBalanceParenthesisSegment(B, n, i + 1, split - 1) ;
      GenerateRandomBalanceParenthesisSegment(B, n, split + 2, j - 1) ;
    }
  }

public:
  DS_Parenthesis() {} 
  ~DS_Parenthesis() {}

  void Free() 
  {
    _rmmTree.Free() ;
    _patRS.Free() ;
  }

  void SetRmmTreeBlockSize(size_t b)
  {
    _rmmTree.SetBlockSize(b) ;
  }

  size_t GetSpace(bool inclusive = true)
  {
    return _rmmTree.GetSpace(false) + (inclusive ? sizeof(*this) : 0) ; 

  }

  void Init(const WORD *B, size_t n, WORD pat, int patLen)
  {
    _rmmTree.Init(B, n) ;
    if (patLen > 0)
      _patRS.Init(B, n, pat, patLen) ;
  }

  // Expose the internal rmmTree.
  const DS_RangeMinMaxTree& GetRmmTree() const
  {
    return _rmmTree ;
  }

  size_t Close(size_t i, const WORD *B, size_t n) const
  {
    // Notice that our FwdSearch include the effect of i, so the d is slightly different than the textbook.
    return _rmmTree.FwdSearch(i, 0, B, n) ;  
  }

  size_t Open(size_t i, const WORD *B, size_t n) const
  {
    return _rmmTree.BwdSearch(i, 0, B, n) ; 
  }

  size_t Enclose(size_t i, const WORD *B, size_t n) const
  {
    return _rmmTree.BwdSearch(i, -1 - Utils::BitRead(B, i), B, n) ;
  }

  bool IsBalance(const WORD *B, size_t n) const
  {
    size_t i ;
    int64_t excess = 0 ;
    for (i = 0 ; i < n ; ++i)
    {
      excess += (2 * Utils::BitRead(B, i) - 1 ) ;
      if (excess < 0)
        return false ;
    }
    return true ;
  }

  size_t PatternRank(size_t i, const WORD *B, size_t n, int inclusive = 1) const
  {
    return _patRS.Rank(i, B, n, inclusive) ;
  }

  size_t PatternSelect(size_t i, const WORD *B, size_t n) const 
  {
    return _patRS.Select(i, B, n) ;
  }

  void GenerateRandomBalanceParenthesis(WORD *B, size_t n, int seed = 17)
  {
    srand(seed) ;
    GenerateRandomBalanceParenthesisSegment(B, n, 0, n - 1) ; 
  }

  void Print(FILE *fp, const WORD *B, size_t n)
  {
    size_t i ;
    for (i = 0 ; i < n ; ++i)
    {
      if (Utils::BitRead(B, i) == 1)
        fprintf(fp, "(") ;
      else
        fprintf(fp, ")") ;
    }
    fprintf(fp, "\n") ;
  }

  void Save(FILE *fp)
  {
    _rmmTree.Save(fp) ;
    _patRS.Save(fp) ;
  }

  void Load(FILE *fp)
  {
    _rmmTree.Load(fp) ;
    _patRS.Load(fp) ;
  }
} ;
}

#endif

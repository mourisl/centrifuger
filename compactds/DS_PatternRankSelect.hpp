#ifndef _MOURISL_COMPACTDS_DS_PATTERN_RANK_SELECT
#define _MOURISL_COMPACTDS_DS_PATTERN_RANK_SELECT

// Binary search based method to calculate  for pattern (not bit, but several bits). 
// The tree structure looks like range min max tree, where we pre-record the information into blocks.
#include "Utils.hpp"

namespace compactds {
class DS_PatternRankSelect
{
private:
  size_t _r ; // number of regions
  size_t _b ; // block size
  size_t _height ; 
  size_t *_counts ; // pattern count in each tree node
  WORD _pat ;
  int _patLen ; // pat is at most 64bits
  size_t _space ;

  // Maps leaf id k to the tree idx
  size_t LeafNum(size_t k) const
  {
    const size_t width = 1ull<<_height ;
    if (k < 2 * _r - width)
      return width - 1 + k ; 
    else
      return width - 1 - _r + k ;
  }

  // Maps tree index to leaf index 
  size_t NumLeaf(size_t v) const
  {
    const size_t width = 1ull<<_height ;
    if (v >= width - 1)
      return v - width + 1 ;
    else
      return v - width + 1 + _r ;
  }

  //@return: the 0-based level of tree index v located at
  int GetLevel(size_t v) const
  {
    return Utils::CountBits(v + 1) - 1 ;
  }

  // The max tree index that is at the same level as v
  size_t LevelMaxNum(size_t v) const
  {
    return (1ull<<(GetLevel(v) + 1)) - 2 ;
  }

  // The min tree index that is at the same level as v
  size_t LevelMinNum(size_t v) const
  {
    return LevelMaxNum(v) / 2 ;
  }

  // Find the tree index containing v, l levels higher
  size_t PromoteLevel(size_t v, int l) const
  {
    return ((v+1) >> l) - 1;
  }

  // Get the leftmost and rightmost leaf id (in leaf idx)
  // Haven't validated yet.
  size_t GetLeftmostLeaf(size_t v) const
  {
    size_t l = GetLevel(v) ;
    size_t diff = v - LevelMinNum(v) ;
    // Each node on this level covers chunk amount of leaves
    size_t chunk = (1ull << (_height - l)) ; 
    size_t ret = (1 << _height) - 1 + diff * chunk ;
    // Pretend this is a complete tree, and then adjust the extra leaves
    if (ret > 2 * _r - 1)
      ret /= 2 ;
    return NumLeaf(ret) ;
  }

  size_t GetRightmostLeaf(size_t v) const
  {
    size_t l = GetLevel(v) ;
    size_t diff = v - LevelMinNum(v) ;
    // Each node on this level covers chunk amount of leaves
    size_t chunk = (1ull << (_height - l)) ; 
    size_t ret = (1 << _height) - 1 + (diff + 1) * chunk - 1 ;
    // Pretend this is a complete tree, and then adjust the extra leaves
    if (ret > 2 * _r - 1)
      ret = (ret - 1) / 2 ;
    return NumLeaf(ret) ;
  }

public:
  DS_PatternRankSelect()
  {
    _space = _b = _r = 0 ;
    _counts = NULL ;
  }

  ~DS_PatternRankSelect()
  {
    Free() ;
  }

  size_t GetSpace(bool inclusive = true)
  {
    return _space + (inclusive ? sizeof(*this) : 0) ;
  }
  
  void Free()
  {
    if (_counts)
    {
      free(_counts) ;
      _counts = NULL ;
      _r = _b = _space = 0 ;
    }
  }

  void SetBlockSize(size_t b)
  {
    _b = b ;
  }

  void Init(const WORD *B, size_t n, WORD pat, int patLen)
  {
    size_t i, j ;
    _pat = pat ;
    _patLen = patLen ;

    if ((int)_b <= patLen)
      // Take the block size as word bits to maintain low space usage
      _b = 16 * WORDBITS  ; 

    _r = DIV_CEIL(n, _b) ; 
    _height = Utils::Log2Ceil(_r) ;
    
    _counts = (size_t *)malloc(sizeof(size_t) * (2 * _r - 1)) ;
    _space = sizeof(*_counts) * (2*_r - 1) ;

    // Fill the leaves 
    for (i = 0 ; i < n ; i += _b)
    {
      size_t l = i / _b ; // leaf id 
      size_t ltid = LeafNum(l) ;
      size_t count = 0 ;
      // the count include the positions that may stretch to the next block
      for (j = i ; j + _patLen - 1 < n && j < i + _b ; ++j)
      {
        WORD w = Utils::BitsRead(B, j, j + _patLen - 1) ;  
        if (w == pat)
          ++count ;
      }
      _counts[ltid] = count ;
    }

    // Fill the internal nodes
    for (i = _r - 2 ; i < _r ; --i)
    {
      size_t count = _counts[2 * i + 1] ;
      if (2 * i + 2 < (2 * _r - 1))
        count += _counts[2 * i + 2] ;
      _counts[i] = count ;
    }
  }

  size_t Rank(size_t i, const WORD *B, size_t n, int inclusive = 1) const
  {
    if (!inclusive)
    {
      if (i == 0)
        return 0 ;
      else
        --i ;
    }
    size_t j ;
    size_t lid = i / _b ;
    size_t rank = 0 ;
    size_t v = 0 ;
  
    while (2 * v + 2 < 2 * _r - 1)
    {
      if (lid <= GetRightmostLeaf(2 * v + 1))
      {
        v = 2 * v + 1 ;
      }
      else
      {
        rank += _counts[2 * v + 1] ;
        v = 2 * v + 2 ;
      }
    }
    
    for (j = i/_b * _b ; j <= i && j + _patLen - 1 < n; ++j)
    {
      // TODO: this part could be optimized by read in one WORD and 
      // use left shift+mask to search the pattern.
      WORD w = Utils::BitsRead(B, j, j + _patLen - 1) ;
      if (w == _pat)
        ++rank ;
    }
    return rank ;
  }

  size_t Select(size_t i, const WORD *B, size_t n) const
  {
    size_t j ;
    size_t v = 0 ;
    size_t count = 0 ;
    while (2 * v + 2 < 2 * _r - 1)
    {
      if (_counts[2 * v + 1] + count >= i)
        v = 2 * v + 1 ;
      else
      {
        count += _counts[2 * v + 1] ;
        v = 2 * v + 2 ;
      }
    }

    size_t k = NumLeaf(v) ;
    for (j = k * _b ; j < (k + 1) * _b ; ++j)
    {
      WORD w = Utils::BitsRead(B, j, j + _patLen - 1) ;
      if (w == _pat)
        ++count ;
      if (count == i)
        return j ; 
    }
    return 0 ; 
  }

  bool IsPattern(size_t i, const WORD *B, size_t n) const
  {
    if (i + _patLen - 1 >= n)
      return false ;
    return (Utils::BitsRead(B, i, i + _patLen - 1) == _pat) ;

  }

  void Save(FILE *fp)
  {
    SAVE_VAR(fp, _r) ;
    SAVE_VAR(fp, _b) ;
    SAVE_VAR(fp, _height) ;
    SAVE_VAR(fp, _pat) ;
    SAVE_VAR(fp, _patLen) ;
    
    if (_r > 0)
      SAVE_ARR(fp, _counts, 2 * _r - 1) ;
  }

  void Load(FILE *fp)
  {
    Free() ;
        
    LOAD_VAR(fp, _r) ;
    LOAD_VAR(fp, _b) ;
    LOAD_VAR(fp, _height) ;
    LOAD_VAR(fp, _pat) ;
    LOAD_VAR(fp, _patLen) ;

    if (_r > 0)
    {
      _counts = (size_t *)malloc(sizeof(size_t) * (2 * _r - 1)) ;
      LOAD_ARR(fp, _counts, 2 * _r - 1) ;
      _space = sizeof(*_counts) * (2*_r - 1) ;
    }
  }
} ;

}
#endif

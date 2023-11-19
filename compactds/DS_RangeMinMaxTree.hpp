#ifndef _MOURISL_COMPACTDS_DS_RANGEMINMMAXTREE
#define _MOURISL_COMPACTDS_DS_RANGEMINMMAXTREE

// Based on section 7.1.1
// Handles the excessive information in bit vector
#include "Utils.hpp"

// Note that the excess can be negative, so we use int64_t in many places.
// The input B is like:
//   1 1 1 0 0 0
// The excess tracking will be
//  0 1 2 3 2 1 0
// (The search for index i is always with respect to between i-th B and (i+1)-th B.)
//  The search for index i is always inclusive.
//   The forward search will return the effect after j-th B.
//   The backward search will return the effect before j-th B.
namespace compactds {
struct _rangeMinMaxTreeNode
{
  // Each number is within the block, so the value range should be small
  // Assume block size is less than 2^15
  int16_t e ; // excess with respect to the beginning of the region
  int16_t min ; // min e
  int16_t max ; // max e
  int16_t n ; // number of times hit min

  void Merge(const struct _rangeMinMaxTreeNode &b)
  {
    if (e + b.min < min)
    {
       min = e + b.min ;
       n = b.n ;
    }
    else if (e + b.min == min)
      n += b.n ;
    
    if (e + b.max > max)
      max = e + b.max ;
    e += b.e ;
  }

  // Wrapper if we need the information from right to left (direction<0)
  int16_t RevE()
  {
    return -e ;
  }

  int16_t RevMin()
  {
    return min <= 0 ? (min - e) : -e ;
  }

  int16_t RevMax()
  {
    return max >= 0 ? (max - e) : -e ;
  }
} ;

class DS_RangeMinMaxTree
{
private:
  size_t _space ;
  size_t _r ; // number of regions
  size_t _b ; // block size
  size_t _n ; 
  size_t _height ;
  struct _rangeMinMaxTreeNode *_tree ; 
  int _cwidth ; // chunk size
  struct _rangeMinMaxTreeNode *_C ; // precomputed chunk

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
    size_t ret = (1ull<<(GetLevel(v) + 1)) - 2 ;
    //if (ret >= 2*_r-1)
    //  ret = 2*_r - 2 ;
    return ret ;
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

  // Update extreme value (min, max)
  int64_t UpdateExtreme(int64_t e, int64_t x, int type) const
  {
    if ((type == 0 && x < e)
        || (type == 1 && x > e))
      return x ;
    return e ;
  }

  void InitPrecomputedChunks()
  {
    size_t i ;
    int j ;
    // Precomputed block
    _C = (struct _rangeMinMaxTreeNode *)malloc(sizeof(*_C) * (1<<_cwidth)) ;
    _space += sizeof(*_C) * (1<<_cwidth) ;

    // Fill precomputed block
    for (i = 0 ; i < (1ull<<_cwidth) ; ++i)
    {
      int16_t excess = 0 ;
      int16_t min = 2, max = -2, minCnt = 0 ;
      for (j = 0 ; j < _cwidth ; ++j)
      {
        excess += (2 * ((i >> j) & 1) - 1) ;
        
        if (excess < min || minCnt == 0)
        {
          min = excess ;
          minCnt = 1 ;
        }
        else if (excess == min)
          ++minCnt ;

        if (excess > max)
          max = excess ; 
      }
      _C[i].e = excess ;
      _C[i].min = min ;
      _C[i].max = max ;
      _C[i].n = minCnt ;
    }
  }

  // Search excess difference after i inside the block containing i (block size is _b)
  // @return: d, or the excess from i to the end of the block when no match. 
  //  retj: the coordinate of the matched position, or just pass the block
  int64_t FwdBlock(size_t i, int64_t d, size_t &retj, const WORD *B) const
  {
    size_t j ;
    
    size_t f = i / _cwidth ; // current chunk. f:from; t: to.
    size_t t = ((i / _b + 1) * _b ) / _cwidth - 1 ; // last chunk in the block 
    
    int64_t excess = 0 ;
    // search the current chunk
    for (j = i ; j < f * _cwidth + _cwidth && j < _n ; ++j)
    {
      excess += (2 * Utils::BitRead(B, j) - 1) ;
      
      if (excess == d)
      {
        retj = j ;
        return d ;
      }
    }
    
    // search the remaining chunks in the block
    size_t p ;
    for (p = f + 1 ; p <= t ; ++p)
    {
      WORD chunk = 0 ;
      if (p * _cwidth >= _n)
        break ;

      if ((p + 1) * _cwidth - 1 <= _n - 1)
        chunk = Utils::BitsRead(B, p * _cwidth, (p+1) * _cwidth -1) ;
      else // the chunk will be padding with 0's (\)), which may cause too small m
        //  but the border case will be handled after searching the hit chunk.
        chunk = Utils::BitsRead(B, p * _cwidth, _n - 1) ;
      
      if ((d <= 0 && excess > d && excess + _C[chunk].min <= d)
          || (d >= 0 && excess < d && excess + _C[chunk].max >= d)) 
          break ;

      excess += _C[chunk].e ; 
    }

    // Could not find it in current block
    if (p > t)
    {
      retj = _cwidth * (t + 1) ; 
      return excess ; 
    }
    if (p * _cwidth >= _n)
    {
      retj = _n ;
      return excess ;
    }
    
    // Search the hit chunk
    for (j = p * _cwidth ; j < p * _cwidth + _cwidth && j < _n ; ++j)
    {
      excess += (2 * Utils::BitRead(B, j) - 1) ;
      
      if (excess == d)
      {
        retj = j ;
        return d ;
      }
    }

    // Can only reach here in border case
    retj = _n ; 
    return excess ;
  }
  
  // Similarly to FwdBlack, but search backwards (to left)
  int64_t BwdBlock(size_t i, int64_t d, size_t &retj, const WORD *B) const
  {
    size_t j ;
    
    size_t f = i / _cwidth ; // current chunk
    size_t t = ((i / _b) * _b ) / _cwidth ; // first chunk in the block 
    
    int64_t excess = 0 ;
    // search the current chunk
    for (j = i ; j >= f * _cwidth && j < _n ; --j)
    {
      excess -= (2 * Utils::BitRead(B, j) - 1) ;
      
      if (excess == d)
      {
        retj = j ;
        return d ;
      }
    }
    
    //search the remaining chunks in the block
    size_t p ;
    for (p = f - 1 ; p >= t && p < _n ; --p)
    {
      WORD chunk = 0 ;

      chunk = Utils::BitsRead(B, p * _cwidth, (p+1) * _cwidth -1) ;
      
      if ((d <= 0 && excess > d && excess + _C[chunk].RevMin() <= d)
          || (d >= 0 && excess < d && excess + _C[chunk].RevMax() >= d)) 
          break ;

      excess += _C[chunk].RevE() ; 
    }

    // Could not find it in current block
    if (p < t)
    {
      retj = _cwidth * t - 1 ; 
      return excess ; 
    }
    if (p >= _n) // not exist 
    {
      retj = _n ;
      return excess ;
    }
    
    // Search the hit chunk
    for (j = p * _cwidth + _cwidth - 1 ; j >= p * _cwidth && j < _n ; --j)
    {
      excess -= (2 * Utils::BitRead(B, j) - 1) ;
      if (excess == d)
      {
        retj = j ;
        return d ;
      }
    }

    // Can only reach here in border case
    retj = _n ; 
    return excess ;
  }

  // Scanning a block for min/max[i,j]
  // Assumes i and j are in the same block
  // As other block searches, it returns the excess in this block,
  //  the extreme value is passed through the reference
  int64_t ExtremeBlock(size_t i, size_t j, int type, int64_t &extreme, const WORD *B) const
  {
    size_t p ; // the index for loop
    
    size_t f = i / _cwidth ; // first chunk 
    size_t t = j / _cwidth ; // last chunk
    
    int64_t excess = 0 ;
    if (type == 0)
      extreme = 2 ;
    else
      extreme = -2 ;
    for (p = i ; p < f * _cwidth + _cwidth && p <= j && p < _n ; ++p)
    {
      excess += (2 * Utils::BitRead(B, p) - 1) ;
      extreme = UpdateExtreme(extreme, excess, type) ;
    }
    if (f == t)
      return excess ;
    
    // Since we search to t-1, there is no need to worry about reading after _n
    for (p = f + 1 ; p <= t - 1 ; ++p) 
    {
      WORD chunk = Utils::BitsRead(B, p * _cwidth, (p+1) * _cwidth - 1) ; 
      extreme = UpdateExtreme(extreme, excess + 
          (type == 0 ? _C[chunk].min : _C[chunk].max), type) ;
      excess += _C[chunk].e ;
    }

    // last chunk
    for (p = t * _cwidth ; p <= j && p <= _n ; ++p)
    {
      excess += (2 * Utils::BitRead(B, p) - 1) ;
      extreme = UpdateExtreme(extreme, excess, type) ;
    }

    return excess ;
  }

  // Given the global min value between [i,j], try to find its count
  // @return: the excess after j or the block if j is out of the bo
  int64_t MinCountBlock(size_t i, size_t j, int64_t min, size_t &minCnt, const WORD *B) const
  {
    size_t p ; // the index for loop

    size_t f = i / _cwidth ; // first chunk 
    size_t t = j / _cwidth ; // last chunk

    int64_t excess = 0 ;
    minCnt = 0 ;
    
    for (p = i ; p < f * _cwidth + _cwidth && p <= j && p < _n ; ++p)
    {
      excess += (2 * Utils::BitRead(B, p) - 1) ;
      if (excess == min)
        ++minCnt ;
    }
    if (f == t)
      return excess ;

    // Since we search to t-1, there is no need to worry about reading after _n
    for (p = f + 1 ; p <= t - 1 ; ++p) 
    {
      WORD chunk = Utils::BitsRead(B, p * _cwidth, (p+1) * _cwidth - 1) ;
      if (excess + _C[chunk].min == min)
        minCnt += _C[chunk].n ;
      excess += _C[chunk].e ;
    }

    // last chunk
    for (p = t * _cwidth ; p <= j && p <= _n ; ++p)
    {
      excess += (2 * Utils::BitRead(B, p) - 1) ;
      if (excess == min)
        ++minCnt ;
    }

    return excess ;
  }

  // Return the excess. The coordinate of the k-th min is returned through selectk. set _n if not found
  int64_t MinSelectBlock(size_t i, size_t j, int64_t min, size_t kthMin, size_t &selectk, size_t &minCnt, const WORD *B) const
  {
    size_t p ; // the index for loop

    size_t f = i / _cwidth ; // first chunk 
    size_t t = j / _cwidth ; // last chunk

    int64_t excess = 0 ;
    minCnt = 0 ;
    
    for (p = i ; p < f * _cwidth + _cwidth && p <= j && p < _n ; ++p)
    {
      excess += (2 * Utils::BitRead(B, p) - 1) ;
      if (excess == min)
      {
        ++minCnt ;
        if (minCnt == kthMin)
        {
          selectk = p ;
          return excess ;
        }
      } 
    }
    if (f == t)
    {
      selectk = _n ;
      return excess ;
    }

    // Since we search to t-1, there is no need to worry about reading after _n
    for (p = f + 1 ; p <= t - 1 ; ++p) 
    {
      WORD chunk = Utils::BitsRead(B, p * _cwidth, (p+1) * _cwidth - 1) ;
      if (excess + _C[chunk].min == min)
      {
        if (kthMin <= minCnt + _C[chunk].n)
        {
          size_t chunki ; 
          for (chunki = p * _cwidth ; chunki < (p + 1) * _cwidth ; ++chunki)
          {
            excess += (2 * Utils::BitRead(B, chunki) - 1) ;
            if (excess == min)
            {
              ++minCnt ;
              if (minCnt == kthMin)
              {
                selectk = chunki ;
                return excess ;
              }
            }
          }
        }
        minCnt += _C[chunk].n ;
      }
      excess += _C[chunk].e ;
    }

    // last chunk
    for (p = t * _cwidth ; p <= j && p <= _n ; ++p)
    {
      excess += (2 * Utils::BitRead(B, p) - 1) ;
      if (excess == min)
      {
        ++minCnt ;
        if (minCnt == kthMin)
        {
          selectk = p ;
          return excess ;
        }
      }
    }

    return excess ;
  }

public:
  DS_RangeMinMaxTree()
  {
    _space = _b = 0 ;
    _cwidth = 8 ;
    _tree = _C = NULL ;
  }

  ~DS_RangeMinMaxTree() 
  {
    Free() ;
  }

  size_t GetSpace(bool inclusive = true)
  {
    return _space + (inclusive ? sizeof(*this) : 0) ; 
  }

  void Free()
  {
    if (_tree != 0)
    {
      free(_tree) ;
      free(_C) ;
      _tree = NULL ;
      _C = NULL ;
    }
  }

  void SetBlockSize(size_t b)
  {
    _b = b ;
  }

  // s: some special to track the count.
  // slen: length of the special character
  void Init(const WORD *B, size_t n)
  {
    size_t i, j ;
    if (_b <= 8) // block size has to be larger than a byte, and should be power of 2
      _b = 1024 ;
        
    _n = n ;
    _r = DIV_CEIL(n, _b) ; 
    _height = Utils::Log2Ceil(_r) ;
    _tree = (struct _rangeMinMaxTreeNode *)malloc(sizeof(*_tree) * (2 * _r - 1)) ;
    _space += sizeof(*_tree) * (2*_r-1) ;

    InitPrecomputedChunks() ;
    
    // Initialize the leafs
    for (i = 0 ; i < n ; i += _b)
    {
      size_t treeIdx = LeafNum(i / _b) ;
      _tree[treeIdx].e = 0 ;
      _tree[treeIdx].min = 2 ;
      _tree[treeIdx].max = -2 ;
      _tree[treeIdx].n = 0 ;
      for (j = i ; j < n && j < i + _b ; j += _cwidth)
      {
        uint64_t chunk = Utils::BitsRead(B, j, j + _cwidth - 1) ;
        _tree[treeIdx].Merge( _C[chunk] ) ;
      }
    }

    // Initialize internal nodes
    for (i = _r - 2 ; i < _r ; --i)
    {
      _tree[i] = _tree[2 * i + 1] ;
      if (2 * i + 2 < (2 * _r - 1))
        _tree[i].Merge(_tree[2 * i + 2]) ;
    }
  }

  //It's a bit different from the textbook that the FwdSearch 
  //  in our implementation include the effect from i.
  //  This makes 0-based indexing have better definition.
  //@return: the index j >= i that has excess different d 
  //        return _n if not found
  size_t FwdSearch(size_t i, int64_t d, const WORD *B, size_t n) const
  {
    size_t j ;
    int64_t excess ;
    excess = FwdBlock(i, d, j, B) ;
    
    if (excess == d)
      return j ;
    if (j == _n)
      return _n ;
    // Not in current block, so we need to search the tree to find the block
    //  v is the tree node index.
    size_t v = LeafNum(i / _b) ;
    
    // Go up the tree first
    // After the iterations, v+1 should be the node containing the target j.
    // The test for v+1>=2*_r-1 is for the case where the last level is not full (leaf level)
    //   if it is rightmost leaf on the last level, we can directly go to parent
    while (v + 1 <= LevelMaxNum(v) && 
        (v + 1 >= 2*_r-1 || (d <= 0 && excess > d && excess + _tree[v + 1].min > d)
         || (d >= 0 && excess < d && excess + _tree[v + 1].max < d) )) 
      // next node block is not enough
      // The next node not necessarily the brother node, but the covered region is adjacent
      //    based on the numbering system.
    {
      if ((v & 1) == 1) // v is left child. Note that our index is 0-based, so it is not the 2*xxx relation in 1-based index. 
        excess += _tree[v + 1].e ; 
      v = (v-1) / 2 ; // parent node
    }
    
    if (v == LevelMaxNum(v)) // Not found. v is the rightmost block on the level
      return _n ;
    
    // Go down the tree to locate the block
    ++v ;
    while (2 * v + 2 < 2 * _r - 1)
    {
      if ((d <= 0 && excess > d && excess + _tree[2 * v + 1].min <= d) 
          || (d >= 0 && excess < d && excess + _tree[2 * v + 1].max >= d))
        v = 2 * v + 1 ;
      else 
      {
        excess += _tree[2 * v + 1].e ;
        v = 2 * v + 2 ;
      }
    }
    
    // The else branch above may go beyond the number of blocks
    if (v >= 2 * _r - 1)  
      return _n ;

    // Search the target leaf block
    // The FwdBlock searches things after the index, so we need put -1 in
    //  NumLeaf(v) * _b
    excess = FwdBlock(NumLeaf(v) * _b, d - excess, j, B) ;
    return j ;
  }

  //@return: the index j <= i that has excess different d comparing with i from right to left 
  //        return _n if not found
  size_t BwdSearch(size_t i, int64_t d, const WORD *B, size_t n) const
  {
    size_t j ;
    int64_t excess ;
    /*if (i == 0)
    {
      // d == 1, b0 == 0
      // or d == -1, b0==1
      if (d == -2 * Utils::BitRead(B, 0) + 1 )
        return 0 ;
      return _n ;
    }*/

    excess = BwdBlock(i, d, j, B) ;
    if (excess == d)
      return j ;
    if (j == _n)
      return _n ;
    size_t v = LeafNum(i / _b) ;
    
    // Go up the tree first
    // After the iterations, v-1 should be the node containing the target j. 
    while (v != 0 && v - 1 >= LevelMinNum(v) && 
        ((d <= 0 && excess > d && excess + _tree[v - 1].RevMin() > d)
         || (d >= 0 && excess < d && excess + _tree[v - 1].RevMax() < d) )) 
      // next node block is not enough
      // The next node not necessarily the brother node, but the covered region is adjacent
      //    based on the numbering system.
    {
      if ((v & 1) == 0) // v is right child. Note that our index is 0-based, so it is not the 2*xxx relation in 1-based index. 
        excess += _tree[v - 1].RevE() ; 
      v = (v-1) / 2 ; // parent node
    }
    
    if (v == LevelMinNum(v)) // Not found. v is the leftmost block on the level
      return _n ;
    
    // Go down the tree to locate the block
    --v ;
    while (2 * v + 2 < 2 * _r - 1)
    {
      if ((d <= 0 && excess > d && excess + _tree[2 * v + 2].RevMin() <= d) 
          || (d >= 0 && excess < d && excess + _tree[2 * v + 2].RevMax() >= d))
        v = 2 * v + 2 ;
      else 
      {
        excess += _tree[2 * v + 2].RevE() ;
        v = 2 * v + 1 ;
      }
    }
    
    // The else branch above may go beyond the number of blocks
    if (v >= 2 * _r - 1)  
      return _n ;

    
    // Search the target leaf block
    // BwdBlock is inclusive
    excess = BwdBlock(NumLeaf(v) * _b + _b - 1, d - excess, j, B) ;
    return j ;
  }

  // type: 0-min, 1-max
  //@return: the min/max value in B[i,j]: included the effects from B[i] and B[j]
  int64_t ExtremeExcess(size_t i, size_t j, int type, const WORD *B, size_t n) const
  {
    int64_t extreme = 0 ;
    int64_t excess = 0 ;

    excess = ExtremeBlock(i, MIN(j, (i / _b) * _b + _b - 1), type, extreme, B);
    if (j/_b <= i / _b) // in the same block
      return extreme ;

    //printf("%d %d: %d %d\n", i, j, extreme, excess) ;
    // Search the tree
    size_t v = LeafNum(i / _b) ;
    size_t l = LeafNum(j / _b) ; 
    int levelv = GetLevel(v) ; 
    int levell = GetLevel(l) ;

    // Upward search
    // v+1 > l: l is in the upper level 
    //  or l is still to the right of v+1
    while (v + 1 > l || v+1 != PromoteLevel(l, levell - levelv)) 
    {
      if ( (v & 1) == 1 && v+1 < 2*_r - 1) //left children
      {
        extreme = UpdateExtreme(extreme, excess + 
            (type == 0 ? _tree[v + 1].min : _tree[v + 1].max), type) ;
        excess += _tree[v + 1].e ;
      }
      v = (v - 1) / 2 ;
      --levelv ;
    }
    //printf("%d %d: %d %d %d\n", i, j, v, extreme, excess) ;

    // Downward search. Now l should be in v+1
    ++v ;
    while (v < _r - 1) // internal nodes 
    {
      if ((type == 0 && extreme <= excess + _tree[v].min)
          || (type == 1 && extreme >= excess + _tree[v].max))
        return extreme ;
      
      if (2 * v + 1 != PromoteLevel(l, levell - (levelv + 1)))
      {
        extreme = UpdateExtreme(extreme, excess + 
            (type == 0 ? _tree[2*v + 1].min : _tree[2*v + 1].max), type) ;
        excess += _tree[2*v + 1].e ;
        v = 2 * v + 2 ;
      }
      else
        v = 2 * v + 1 ;
      ++levelv ;
    }
    //printf("%d %d: %d %d. %d %d\n", i, j, v, extreme, excess, _tree[v].min) ;

    if ((type == 0 && extreme <= excess + _tree[v].min)
        || (type == 1 && extreme >= excess + _tree[v].max))
      return extreme ;

    // last block
    int64_t lastExtreme = 0 ;
    ExtremeBlock((j / _b) * _b, j, type, lastExtreme, B) ;
    //printf("%d %d: %d vs %d %d\n", i, j, extreme, excess, lastExtreme) ;

    return UpdateExtreme(extreme, excess + lastExtreme, type) ;
  }
  
  // Leftmost position of a minimum in excess(B, i, j)
  size_t Rmq(size_t i, size_t j, const WORD *B, size_t n) const
  {
    int64_t min = ExtremeExcess(i, j, 0, B, n) ; 
    return FwdSearch(i, min, B, n) ;
  }

  // Leftmost position of a maximum in excess(B, i, j)
  size_t RMq(size_t i, size_t j, const WORD *B, size_t n) const
  {
    int64_t max = ExtremeExcess(i, j, 1, B, n) ; 
    return FwdSearch(i, max, B, n) ;
  }

  // Need .maxn in node structure to support maxcount but not implemented now,
  //   depends on future application to decide whether implement this feature.
  size_t MinCount(size_t i, size_t j, const WORD *B, size_t n) const
  {
    // The min in this whole range
    int64_t min = ExtremeExcess(i, j, 0, B, n) ; 
    
    // Get first block
    size_t minCnt = 0 ;
    int64_t excess = 0 ;

    excess = MinCountBlock(i, MIN(j, (i / _b) * _b + _b - 1), min, minCnt, B);
    if (j/_b <= i / _b) // in the same block
      return minCnt ;
    //printf("%d: %d %d\n", i, min, minCnt) ;
    // Search the tree
    size_t v = LeafNum(i / _b) ;
    size_t l = LeafNum(j / _b) ; 
    int levelv = GetLevel(v) ; 
    int levell = GetLevel(l) ;

    // Upward search
    // v+1 > l: l is in the upper level 
    //  or l is still to the right of v+1
    while (v + 1 > l || v+1 != PromoteLevel(l, levell - levelv)) 
    {
      if ( (v & 1) == 1 && v+1 < 2*_r - 1) //left children
      {
        if (excess + _tree[v + 1].min == min)
          minCnt += _tree[v + 1].n ;
        excess += _tree[v + 1].e ;
      }
      v = (v - 1) / 2 ;
      --levelv ;
    }
    //printf("%d: %d %d\n", i, v, excess) ;
    // Downward search. Now l should be in v+1
    ++v ;
    while (v < _r - 1) // internal nodes 
    {
      if (min < excess + _tree[v].min)
        return minCnt ;
      
      if (2 * v + 1 != PromoteLevel(l, levell - (levelv + 1)))
      {
        if (excess + _tree[2*v + 1].min == min)
          minCnt += _tree[2*v + 1].n ;
        excess += _tree[2*v + 1].e ;
        v = 2 * v + 2 ;
      }
      else
        v = 2 * v + 1 ;
      ++levelv ;
    }

    //printf("%d: %d %d. %d %d %d\n", i, min, minCnt, excess, v, _tree[v].min) ;
    if (min < excess + _tree[v].min)
      return minCnt ;
    
    // last block
    size_t lastMinCnt = 0 ;
    // Notice the we need to use min-excess here to adjust the excess so far.
    MinCountBlock((j / _b) * _b, j, min - excess, lastMinCnt, B) ;

    //printf("%d: %d %d %d\n", i, min, minCnt, lastMinCnt) ;
    return minCnt + lastMinCnt ;
  }

  // Select the t-th (1-based) minimum element in B[i..j]
  size_t MinSelect(size_t i, size_t j, size_t t, const WORD *B, size_t n) const
  {
    // The min in this whole range
    int64_t min = ExtremeExcess(i, j, 0, B, n) ; 
    
    // Get first block
    size_t minCnt = 0 ;
    int64_t excess = 0 ;
    size_t ret = _n ;

    excess = MinSelectBlock(i, MIN(j, (i / _b) * _b + _b - 1), min, t, ret, minCnt, B);
    //printf("%d: %d %d\n", i, ret, minCnt) ;
    if (j/_b <= i / _b // in the same block
        || ret < _n ) // already found
      return ret ;
    
    // Search the tree
    size_t v = LeafNum(i / _b) ;
    size_t l = LeafNum(j / _b) ; 
    int levelv = GetLevel(v) ; 
    int levell = GetLevel(l) ;

    //printf("%d: %d %d. %d %d\n", i, min, minCnt, v, _r) ;
    // Upward search
    // v+1 > l: l is in the upper level 
    //  or l is still to the right of v+1
    while (v + 1 > l || v+1 != PromoteLevel(l, levell - levelv)) 
    {
      if (v + 1 < 2 * _r - 1 && excess + _tree[v + 1].min == min
          && minCnt + _tree[v + 1].n >= t)
          break ;

      if ( (v & 1) == 1 && v+1 < 2*_r - 1) //left children
      {
        if (excess + _tree[v + 1].min == min)
          minCnt += _tree[v + 1].n ;
        excess += _tree[v + 1].e ;
      }
      v = (v - 1) / 2 ;
      --levelv ;
    }
    
    if (v == LevelMaxNum(v)) // Not found. v is the rightmost block on the level
      return _n ;

    //printf("%d: %d %d %d\n", i, v, excess, minCnt) ;
    // Downward search. 
    ++v ;
    while (v < _r - 1) // internal nodes 
    {
      if (min < excess + _tree[v].min)
        return ret ;
      
      if ( 2 * v + 1 != PromoteLevel(l, levell - (levelv + 1)) //j is in the right chilad 
          && (excess + _tree[2 * v + 1].min != min 
            || minCnt + _tree[2 * v+1].n < t)) // left child could not reach t.
      {
        if (excess + _tree[2*v + 1].min == min)
          minCnt += _tree[2*v + 1].n ;
        excess += _tree[2*v + 1].e ;
        v = 2 * v + 2 ;
      }
      else
        v = 2 * v + 1 ;

      ++levelv ;
    }

    //printf("%d: %d %d. %d %d %d\n", i, min, minCnt, excess, v, _tree[v].min) ;
    if (min < excess + _tree[v].min)
      return _n ;
    
    // last block
    size_t lastMinCnt = 0 ;
    // Notice the we need to use min-excess here to adjust the excess so far.
    v = NumLeaf(v) ;
    MinSelectBlock(v * _b, MIN(j, (v+1) * _b - 1), min - excess, t - minCnt, ret, lastMinCnt, B) ;

    //printf("%d: %d %d %d %d %d\n", i, v, min, minCnt, lastMinCnt, ret) ;
    return ret ;
  }
  
  void Save(FILE *fp)
  {
    SAVE_VAR(fp, _r) ;
    SAVE_VAR(fp, _b) ;
    SAVE_VAR(fp, _n) ;
    SAVE_VAR(fp, _height) ;
    SAVE_VAR(fp, _cwidth) ;
    SAVE_ARR(fp, _tree, 2 * _r - 1) ;
    SAVE_ARR(fp, _C, 1 << _cwidth) ;
  }

  void Load(FILE *fp)
  {
    Free() ;
    LOAD_VAR(fp, _r) ;
    LOAD_VAR(fp, _b) ;
    LOAD_VAR(fp, _n) ;
    LOAD_VAR(fp, _height) ;
    LOAD_VAR(fp, _cwidth) ;
    _tree = (struct _rangeMinMaxTreeNode *)malloc(sizeof(*_tree) * (2 * _r - 1)) ;
    LOAD_ARR(fp, _tree, 2 * _r - 1) ;
    _C = (struct _rangeMinMaxTreeNode *)malloc(sizeof(*_C) * (1<<_cwidth)) ;
    LOAD_ARR(fp, _C, 1 << _cwidth) ;
    _space = sizeof(struct _rangeMinMaxTreeNode) * (2 * _r - 1 + (1<<_cwidth)) ;
  }
} ;
}

#endif

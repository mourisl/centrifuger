#ifndef _MOURISL_COMPACTDS_DIFFERENCECOVER
#define _MOURISL_COMPACTDS_DIFFERENCECOVER

#include "Utils.hpp"
#include "SimpleVector.hpp"

#include <map>
#include <algorithm>

// The class handling difference covers
// Difference cover is a set of numbers D={a_0, ... a_{m-1}} in range [0, v)
//  such that every i in [0,1) there is some a_j, a_k in D s.t. i=(a_j-a_k)%v.
//  So the name comes from the differences of a set can cover all the element.
// This class also handles when query elements larger than _v (cyclic difference cover)
namespace compactds {
class DifferenceCover
{
private:
  int _v ; // period size  
  int *_dcs ; // DCs
  int _m ; // number of DCs 
  std::map<int, int> _dcMap ; // maybe replace this with a bit vector later 
  int *_precomputedD ; // precomputed information for Delta query

  int GetB(int i, int r)
  {
    if (i < r)
      return 1 ;
    else if (i < r + 1)
      return r + 1 ;
    else if (i < 2 * r + 1)
      return 2 * r + 1;
    else if (i < 4 * r + 2)
      return 4 * r + 3 ;
    else if (i < 5 * r + 3)
      return 2 * r + 2 ;
    else if (i < 6 * r + 3)
      return 1 ;
    else
      return 0 ; // ERROR
  }
public:
  DifferenceCover() 
  {
    _v = 4096 ;
    _dcs = NULL ;
    _m = 0 ;
  }

  ~DifferenceCover() 
  {
    Free() ;
  }

  void Free()
  {
    if (_dcs)
    {
      free(_dcs) ;
      free(_precomputedD) ;
      _dcs = NULL ;
      _precomputedD = NULL ;
      _m = 0 ;
    }
  }

  // The construction is based on Colbourn, Ling 2000
  void Init(int v)
  {
    int i ;
    if (_v <= 13)
      _v = 14 ;

    this->_v = _v ;
    // Use the Colbourn, Ling method to find the cover 
    int r = CEIL((-36 + sqrt(1296 - 96*(13 - v)))/48.0) ;
    SimpleVector<int> rawdcs ; 
    rawdcs.Reserve(6 * r + 4) ;
    rawdcs.PushBack(0) ;
    for (i = 1 ; i <= 6 * r + 3 ; ++i)
      rawdcs.PushBack( rawdcs[i - 1] + GetB(i - 1, r)) ;
    
    // Put the finalized difference cover
    _m = 0 ;
    for (i = 0 ; i < 6 * r + 4 ; ++i)
    {
      int dc = rawdcs[i] % _v ;
      if (_dcMap.find(dc) == _dcMap.end())
      {
        _dcMap[dc] = _m ;
        ++_m ;
      }
    }
    
    _dcs = (int *)malloc(sizeof(_dcs[0]) * _m) ;
    i = 0 ;
    for (std::map<int, int>::iterator it = _dcMap.begin() ; it != _dcMap.end() ; ++it, ++i)
    {
      _dcs[i] = it->first ;
    }
    
    // Reorder them into increasing order
    std::sort(_dcs, _dcs + _m) ;
    for (i = 0 ; i < _m ; ++i)
    {
      _dcMap[_dcs[i]] = i ;
    }
    
    // Precompute the look up table d for Delta query
    // Lemma 4 in Fast Lightweight Suffix Array Construction and Checking 
    // We can enumerate all the differences from D 
    int j ;
    _precomputedD = (int *)malloc(sizeof(_precomputedD[0]) * _v) ;
    memset(_precomputedD, -1, sizeof(_precomputedD[0]) * _v) ;
    _precomputedD[0] = 0 ;
    for (i = 0 ; i < _m ; ++i)
    {
      for (j = 0 ; j < _m ; ++j)
      {
        int d = _dcs[j] - _dcs[i] ;
        if (d < 0)
          d += _v ;
        _precomputedD[d] = _dcs[i] ; 
      }
    }
  }

  static size_t EstimateCoverSize(int v)
  {
    if (v <= 13)
      return POSITIVE_INF ;
    int r = CEIL((-36 + sqrt(1296 - 96*(13 - v)))/48.0) ;
    return 6 * r + 4 ;
  }

  // Check whether an element is in diff-cover
  bool IsInDC(size_t i)
  {
    if (_dcMap.find(i%_v) != _dcMap.end())
      return true ;
    return false ;
  }

  int GetV()
  {
    return _v ;
  }

  // Get the size of the DC that can cover [0, n)
  size_t GetSize(size_t n)
  {
    int i ;
    for (i = 0 ; i < _m ; ++i)
    {
      if (_dcs[i] >= (int)(n % _v))
        break ;
    }
    return n / _v * _m + i ;
  }

  // Return the difference cover in a list to cover [0, n)
  size_t GetDiffCoverList(size_t n, size_t *dcList)
  {
    int i ;
    size_t c ;
    size_t cycleCnt = DIV_CEIL(n, _v) ;
    size_t ret = 0 ;
    for (c = 0 ; c < cycleCnt ; ++c)
    {
      for (i = 0 ; i < _m ; ++i)
      {
        size_t x = c * _v + _dcs[i] ;
        if (x >= n)
          break ;
        dcList[ret] = x ;
        ++ret ;
      }
    }
    return ret ;
  }

  // Return the index when skipping the non-DC elements in the list
  //  Assume i is in the difference cover.
  size_t CompactIndex(size_t i)
  {
    return i / _v * _m + _dcMap[i % _v] ;
    //int k = _dcMap[i%v] ;
    //return (n / v) * k + (k < coverCntInLastCycle ? k : coverCntInLastCycle) + i / _v ; 
  }

  // Return the offset delta that (i+delta)%_v and (j+delta)%_v is in the difference cover
  // There are two such offset, depending on the order of i,j, and we select the smaller one. This also makes the selection symmetric
  int Delta(size_t i, size_t j)
  {
    int ri = i % _v ;
    int rj = j % _v ;

    int d = (rj - ri)%_v ;
    if (d < 0)
      d += _v ;
    d = (_precomputedD[d] - ri)%_v ;
    if (d < 0)
      d += _v ; 

    int d2 = (ri - rj)%_v ;
    if (d2 < 0)
      d2 += _v ;
    d2 = (_precomputedD[d2] - rj)%_v ;
    if (d2 < 0)
      d2 += _v ; 

    if (d2 < d)
      return d2 ;
    else
      return d ;     
  }

  void Save(FILE *fp)
  {
    SAVE_VAR(fp, _v) ; 
    SAVE_VAR(fp, _m) ; // number of DCs 
    SAVE_ARR(fp, _dcs, _m) ; // DCs
    SAVE_ARR(fp, _precomputedD, _v) ;
  }

  void Load(FILE *fp)
  {
    Free() ;
    
    LOAD_VAR(fp, _v) ; 
    LOAD_VAR(fp, _m) ; // number of DCs 
    _dcs = (int *)malloc(sizeof(_dcs[0]) * _m) ;
    LOAD_ARR(fp, _dcs, _m) ; // DCs
    _precomputedD = (int *)malloc(sizeof(_precomputedD[0]) * _v) ;
    LOAD_ARR(fp, _precomputedD, _v) ;
    
    _dcMap.clear() ;
    for (int i = 0 ; i < _m ; ++i)
    {
      _dcMap[_dcs[i]] = i ;
    }
  }
} ;
}

#endif

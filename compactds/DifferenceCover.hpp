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
// This class also handles when query elements larger than v (cyclic difference cover)
namespace compactds {
class DifferenceCover
{
private:
  int v ; // period size  
  int *dcs ; // DCs
  int m ; // number of DCs 
  std::map<int, int> dcMap ; // maybe replace this with a bit vector later 
  int *precomputedD ; // precomputed information for Delta query

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
    v = 4096 ;
    dcs = NULL ;
    m = 0 ;
  }

  ~DifferenceCover() 
  {
    if (dcs)
    {
      free(dcs) ;
      free(precomputedD) ;
    }
  }

  // The construction is based on Colbourn, Ling 2000
  void Init(int v)
  {
    int i ;
    if (v <= 13)
      v = 14 ;

    this->v = v ;
    // Use the Colbourn, Ling method to find the cover 
    int r = CEIL((-36 + sqrt(1296 - 96*(13 - v)))/48.0) ;
    SimpleVector<int> rawdcs ; 
    rawdcs.Reserve(6 * r + 4) ;
    rawdcs.PushBack(0) ;
    for (i = 1 ; i <= 6 * r + 3 ; ++i)
      rawdcs.PushBack( rawdcs[i - 1] + GetB(i - 1, r)) ;
    
    // Put the finalized difference cover
    m = 0 ;
    for (i = 0 ; i < 6 * r + 4 ; ++i)
    {
      int dc = rawdcs[i] % v ;
      if (dcMap.find(dc) == dcMap.end())
      {
        dcMap[dc] = m ;
        ++m ;
      }
    }
    
    dcs = (int *)malloc(sizeof(dcs[0]) * m) ;
    i = 0 ;
    for (std::map<int, int>::iterator it = dcMap.begin() ; it != dcMap.end() ; ++it, ++i)
    {
      dcs[i] = it->first ;
    }
    
    // Reorder them into increasing order
    std::sort(dcs, dcs + m) ;
    for (i = 0 ; i < m ; ++i)
    {
      dcMap[dcs[i]] = i ;
    }
    
    // Precompute the look up table d for Delta query
    // Lemma 4 in Fast Lightweight Suffix Array Construction and Checking 
    // We can enumerate all the differences from D 
    int j ;
    precomputedD = (int *)malloc(sizeof(precomputedD[0]) * v) ;
    memset(precomputedD, -1, sizeof(precomputedD[0]) * v) ;
    precomputedD[0] = 0 ;
    for (i = 0 ; i < m ; ++i)
    {
      for (j = 0 ; j < m ; ++j)
      {
        int d = dcs[j] - dcs[i] ;
        if (d < 0)
          d += v ;
        precomputedD[d] = dcs[i] ; 
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
    if (dcMap.find(i%v) != dcMap.end())
      return true ;
    return false ;
  }

  int GetV()
  {
    return v ;
  }

  // Get the size of the DC that can cover [0, n)
  size_t GetSize(size_t n)
  {
    int i ;
    for (i = 0 ; i < m ; ++i)
    {
      if (dcs[i] >= (int)(n % v))
        break ;
    }
    return n / v * m + i ;
  }

  // Return the difference cover in a list to cover [0, n)
  size_t GetDiffCoverList(size_t n, size_t *dcList)
  {
    int i ;
    size_t c ;
    size_t cycleCnt = DIV_CEIL(n, v) ;
    size_t ret = 0 ;
    for (c = 0 ; c < cycleCnt ; ++c)
    {
      for (i = 0 ; i < m ; ++i)
      {
        size_t x = c * v + dcs[i] ;
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
    return i / v * m + dcMap[i % v] ;
    //int k = dcMap[i%v] ;
    //return (n / v) * k + (k < coverCntInLastCycle ? k : coverCntInLastCycle) + i / v ; 
  }

  // Return the offset delta that (i+delta)%v and (j+delta)%v is in the difference cover
  // There are two such offset, depending on the order of i,j, and we select the smaller one. This also makes the selection symmetric
  int Delta(size_t i, size_t j)
  {
    int ri = i % v ;
    int rj = j % v ;

    int d = (rj - ri)%v ;
    if (d < 0)
      d += v ;
    d = (precomputedD[d] - ri)%v ;
    if (d < 0)
      d += v ; 

    int d2 = (ri - rj)%v ;
    if (d2 < 0)
      d2 += v ;
    d2 = (precomputedD[d2] - rj)%v ;
    if (d2 < 0)
      d2 += v ; 

    if (d2 < d)
      return d2 ;
    else
      return d ;     
  }
} ;
}

#endif

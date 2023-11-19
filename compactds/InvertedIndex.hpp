#ifndef _MOURISL_COMPACTDS_INVERTEDINDEX
#define _MOURISL_COMPACTDS_INVERTEDINDEX

// Use permutation to represent inverted index

#include "Utils.hpp"
#include "FixedSizeElemArray.hpp"
#include "Permutation.hpp"
#include "Bitvector_Plain.hpp"
#include "CompactMapper.hpp"

namespace compactds {
class InvertedIndex
{
private:
  size_t _n ;
  Permutation _pi ; 
  Bitvector_Plain _D ; // marker of the start position for each number/alphabet in the concatendated permutation list.
  CompactMapper _map ;
  size_t _space ;

public:
  InvertedIndex() 
  {
  }

  ~InvertedIndex()
  {
  }

  size_t GetSpace(bool inclusive = true)
  {
    return _space + (inclusive ? sizeof(*this) : 0) ; 
  }

  void Init(const FixedSizeElemArray &list, size_t n, bool sparseMap)
  {
    size_t i ;
    _n = n ;

    _map.Init(list, n, sparseMap) ;
     
    size_t *pi = (size_t *)malloc(sizeof(*pi) * _n) ;
    size_t *psum = (size_t *)calloc(sizeof(*psum), _n) ;
    for (i = 0 ; i < _n ; ++i)
    {
      ++psum[ _map.Map(list.Read(i)) ] ;
    }
  
    size_t m = _map.GetCompactSize() ;

    _D.Malloc(_n) ;
    _D.BitSet(0) ;
    for (i = 1 ; i < m ; ++i)
    {
      psum[i] += psum[i - 1] ;
      _D.BitSet(psum[i - 1]) ;
    }
    for (i = m - 1 ; i > 0 ; --i)
      psum[i] = psum[i - 1] ;
    psum[0] = 0 ;
    _D.Init() ;

    for (i = 0 ; i < _n ; ++i)
    {
      size_t tmp = _map.Map(list.Read(i)) ;
      pi[ psum[tmp] ] = i ;
      ++psum[tmp] ;
    }
    _pi.Init(pi, n) ;

    free(pi) ;
    free(psum) ;
  }

  // Search the ith occurence label l (0-based)
  size_t Search(size_t l, size_t i) const 
  {
    size_t mapl = _map.Map(l) ;
    return _pi.Next( _D.Select(mapl + 1) + i) ;  
  }

  // @return: the number of positions for label l
  size_t Positions(size_t l, std::vector<size_t> &pos) const
  {
    size_t mapl = _map.Map(l) ;
    size_t i, cnt ;
    if (mapl == _map.GetCompactSize() - 1)
      cnt = _n - _D.Select(mapl + 1) ;
    else
      cnt = _D.Select(mapl + 2) - _D.Select(mapl) ;

    size_t start = _D.Select(mapl + 1) ;
    for (i = 0 ; i < cnt ; ++i)
      pos.push_back( _pi.Next(start + i) ) ;

    return cnt ;
  }

  // Count the number of label l in the sequences
  size_t Count(size_t l) const
  {
    size_t mapl = _map.Map(l) ;
    if (mapl == _map.GetCompactSize() - 1)
      return _n - _D.Select(mapl + 1) ;
    else
      return _D.Select(mapl + 2) - _D.Select(mapl) ;
  }

  void Save(FILE *fp)
  {
    SAVE_VAR(fp, _n) ;    
    SAVE_VAR(fp, _space) ;   
    _pi.Save(fp) ;
    _D.Save(fp) ;
    _map.Save(fp) ;
  }

  void Load(FILE *fp)
  {
    LOAD_VAR(fp, _n) ;
    LOAD_VAR(fp, _space) ;
    _pi.Load(fp) ;
    _D.Load(fp) ;
    _map.Load(fp) ;
  }
} ;

}

#endif

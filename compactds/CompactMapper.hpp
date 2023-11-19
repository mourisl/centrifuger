#ifndef _MOURISL_COMPACTDS_COMPACTMAPPER
#define _MOURISL_COMPACTDS_COMPACTMAPPER

// Map a set of m distinct elements to [0,m-1]
#include <algorithm>
#include <map>
#include <vector>

#include "FixedSizeElemArray.hpp"
#include "Bitvector_Plain.hpp"
#include "Bitvector_Sparse.hpp"

namespace compactds {
class CompactMapper
{
private:
  bool _sparse ; // whether use sparse representation
  Bitvector_Plain _P ;
  Bitvector_Sparse _S ;
  size_t _m ;
public:
  CompactMapper()
  {
  }

  ~CompactMapper()
  {
    Free() ;
  }

  size_t GetSpace(int inclusive = true)
  {
    return _P.GetSpace() - sizeof(_P) + _S.GetSpace() - sizeof(_S) + (inclusive ? sizeof(*this) : 0) ;
  }
  
  void Free()
  {
    _P.Free() ;
    _S.Free() ;
  }

  void Init(const FixedSizeElemArray &a, size_t n, bool sparse)
  {
    size_t i ;
    _sparse = sparse ;
    if (sparse)
    {
      std::map<size_t, size_t> reduceMap ;
      std::vector<size_t> elems ;
      size_t max = 0 ;
      for (i = 0 ; i < n ; ++i)
      {
        size_t tmp = a.Read(i) ;
        if (reduceMap.find(tmp) == reduceMap.end())
        {
          reduceMap[tmp] = i ;
          elems.push_back(tmp) ;
          if (tmp > max)
            max = tmp ;
        }
      }
      std::sort(elems.begin(), elems.end()) ;
      _m = elems.size() ;
      _S.InitFromOnes(elems.data(), max + 1, _m) ;
    }
    else
    {
      size_t max = 0 ;
      for (i = 0 ; i < n ; ++i)
      {
        size_t tmp = a.Read(i) ;
        if (tmp > max)
          max = tmp ;
      }
      
      _P.Malloc(max + 1) ;
      for (i = 0 ; i < n ; ++i)
      {
        size_t tmp = a.Read(i) ;
        _P.BitSet(tmp) ;
      }
      _P.Init() ;
      _m = _P.Rank1(max) ;
    }
  }

  size_t GetCompactSize() const
  {
    return _m ;
  }

  size_t Map(size_t v) const
  {
    if (_sparse)
      return _S.Rank1(v, 0) ;
    else
      return _P.Rank1(v, 0) ;
  }

  size_t MapBack(size_t i) const
  {
    if (_sparse)
      return _S.Select(i + 1) ;
    else
      return _P.Select(i + 1) ;
  }

  void Save(FILE *fp)
  {
    SAVE_VAR(fp, _sparse) ;
    SAVE_VAR(fp, _m) ;
    if (_sparse)
      _S.Save(fp) ;
    else
      _P.Save(fp) ;
  }

  void Load(FILE *fp)
  {
    LOAD_VAR(fp, _sparse) ;
    LOAD_VAR(fp, _m) ;
    if (_sparse)
      _S.Load(fp) ;
    else
      _P.Load(fp) ;
  }
} ;
} 

#endif

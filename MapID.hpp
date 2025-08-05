#ifndef _MOURISL_MAPID
#define _MOURISL_MAPID

#include <map>
#include <vector>

#include <stdint.h>

// This class that maps of arbitrary object to numeric ID in range of [0, n)
template <class T>
class MapID
{
private:
  std::map<T, size_t> _toNumId ; 
  std::vector<T> _toOrigElem ;
public:
  MapID() {}
  ~MapID() {}

  // @return: mapped id. Return the 
  size_t Add(const T &elem)
  {
    if (_toNumId.find(elem) == _toNumId.end()) 
    {
      uint64_t id = _toNumId.size() ;
      _toNumId[elem] = id ;
      _toOrigElem.push_back(elem) ;
      return id ;
    }
    else
      return _toNumId[elem] ;
  }

  void Clear()
  {
    _toNumId.clear() ;
    _toOrigElem.clear() ;
  }

  size_t Map(const T &elem)
  {
    return _toNumId[elem] ;
  }

  bool IsIn(const T &elem)
  {
    return _toNumId.find(elem) != _toNumId.end() ; 
  }

  // Map to original value
  T Inverse(uint64_t nid)
  {
    return _toOrigElem[nid] ;
  }

  void GetElemList(std::vector<T> &l)
  {
    l = _toOrigElem ;
  }

  size_t GetSize()
  {
    return _toOrigElem.size() ;
  }

  void Save(FILE *fp)
  {
    size_t n = GetSize() ;
    fwrite(&n, sizeof(n), 1, fp) ;
    fwrite(_toOrigElem.data(), sizeof(T), n, fp) ;
  }

  void Load(FILE *fp)
  {
    size_t n ;
    uint64_t *list ;
    fread(&n, sizeof(n), 1, fp) ;
    list = new T[n] ;
    fread(list, sizeof(T), n, fp) ;

    _toOrigElem.clear() ;
    _toOrigElem.insert(_toOrigElem.end(), list, list + n) ;
    
    _toNumId.clear() ;
    for (size_t i = 0 ; i < n ; ++i)
      _toNumId[ _toOrigElem[i] ] = i ;
    delete[] list ;
  }
} ;

#endif


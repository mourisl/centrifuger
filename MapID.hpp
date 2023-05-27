#ifndef _MOURISL_MAPID
#define _MOURISL_MAPID

#include <map>
#include <vector>
#include <stdint.h>
#include <fstream>

// This class that maps of arbitrary object to numeric ID in range of [0, n)
template <class T>
class MapID
{
private:
  std::map<T, uint64_t> _toNumId ; 
  std::vector<T> _toOrigElem ;
public:
  MapID() {}
  ~MapID() {}
  
  void Clear()
  {
    _toNumId.clear() ;
    _toOrigElem.clear() ;
  }

  // @return: mapped id
  uint64_t Add(const T &elem)
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

  uint64_t Map(const T &elem)
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

  size_t GetSize()
  {
    return _toOrigElem.size() ;
  }

  void Save(std::ofstream &ofs)
  {
    size_t n = GetSize() ;
    ofs << n ;
    ofs.write((char *)_toOrigElem.data(), sizeof(T) * n) ;
  }

  void Load(std::ifstream &ifs)
  {
    size_t n ;
    T *list ;
    ifs >> n ;
    list = new T[n] ;
    ifs.read((char *)list, sizeof(T) * n) ;

    _toOrigElem.clear() ;
    _toOrigElem.insert(_toOrigElem.end(), list, list + n) ;
    
    _toNumId.clear() ;
    for (size_t i = 0 ; i < n ; ++i)
      _toNumId[ _toOrigElem[i] ] = i ;
    delete[] list ;
  }
} ;

#endif


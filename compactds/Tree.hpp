#ifndef _MOURISL_COMPACTDS_TREE
#define _MOURISL_COMPACTDS_TREE

#include "Utils.hpp"

namespace compactds {
  
class Tree 
{
protected:
  size_t _space ;
  size_t _n ; // number of nodes in tree

public:
  Tree() 
  {
    _space = 0 ;
    _n = 0 ;
  }
  ~Tree() {}

  virtual size_t GetSpace(bool inclusive) = 0 ;

  virtual size_t Root() const = 0 ;
  virtual size_t ChildSelect(size_t v, size_t t) const = 0 ; // Get the t-th (1-based) child of v
  virtual size_t FirstChild(size_t v) const = 0 ;
  virtual size_t LastChild(size_t v) const = 0 ;
  virtual size_t ChildrenCount(size_t v) const = 0 ;
  virtual size_t ChildRank(size_t v) const = 0 ; // Rank is always 1-based

  virtual size_t NextSibling(size_t v) const = 0 ;
  virtual size_t PrevSibling(size_t v) const = 0 ;

	virtual size_t Parent(size_t v) const = 0 ;
  virtual bool IsLeaf(size_t v) const = 0 ;
  
	virtual size_t NodeMap(size_t v) const = 0 ;
  virtual size_t NodeSelect(size_t i) const = 0 ;
  
  // Whether u is an ancestor of v.
  virtual bool IsAncestor(size_t u, size_t v) const
  {
    size_t p = v ;
    while (p != Root() && p != u)
      p = Parent(p) ;
    if (p == u)
      return true ;
    else
      return false ;
  }

  virtual size_t Depth(size_t v) const
  {
    if (v == Root())
      return 0 ;
    size_t ret = 1 ;
    size_t p = Parent(v) ;
    while (p != Root())
    {
      p = Parent(p) ;
      ++ret ;
    }
    return ret ;
  }

  virtual size_t LeafCountInSubTree(size_t v) const
  {  
    if (IsLeaf(v))
      return 1 ;
    size_t i ;
    size_t childCnt = ChildrenCount(v) ;
    size_t ret = 0 ;
    for (i = 0 ; i < childCnt ; ++i)
      ret += LeafCountInSubTree( ChildSelect(v, i + 1) ) ;
    return ret ;
  }

  virtual size_t SubTreeSize(size_t v) const
  {  
    if (IsLeaf(v))
      return 1 ;
    size_t i ;
    size_t childCnt = ChildrenCount(v) ;
    size_t ret = 0 ;
    for (i = 0 ; i < childCnt ; ++i)
      ret += SubTreeSize( ChildSelect(v, i + 1) ) ;
    return ret + 1 ;
  }

  virtual bool IsFirstChild(size_t v) const
  {
    if (v == Root())
      return true ;
    if (ChildRank(v) == 1)
      return true ;
    return false ;
  }

  virtual bool IsLastChild(size_t v) const
  {
    if (v == Root())
      return true ;

    size_t p = Parent(v) ;
    size_t pChildCnt = ChildrenCount(p) ;
    if (ChildRank(v) == pChildCnt)
      return true ;
    return false ;
  }
  
  virtual size_t LCA(size_t u, size_t v) const
  {
    SimpleVector<size_t, size_t> upath ; 
    SimpleVector<size_t, size_t> vpath ;

    size_t p ;
    
    upath.PushBack(u) ;
    p = Parent(u) ;
    while (p != 0)
    {
      upath.PushBack(p) ;
      p = Parent(p) ;
    }
    upath.PushBack(0) ;
    
    vpath.PushBack(v) ;
    p = Parent(v) ;
    while (p != 0)
    {
      vpath.PushBack(p) ;
      p = Parent(p) ;
    }
    vpath.PushBack(0) ;
    
    upath.Reverse() ;
    vpath.Reverse() ;

    size_t size = MIN(upath.Size(), vpath.Size()) ;
    size_t i ;
    for (i = 0 ; i < size; ++i)
      if (upath[i] != vpath[i])
        break ;
    return upath[i - 1] ;
  }

  size_t GetSize() const
  {
    return _n ;
  }

  virtual void Save(FILE *fp) 
  {
    SAVE_VAR(fp, _space) ;
    SAVE_VAR(fp, _n) ;
  }

  virtual void Load(FILE *fp)
  {
    LOAD_VAR(fp, _space) ;
    LOAD_VAR(fp, _n) ;
  }
} ;

} // end of namespace

#endif

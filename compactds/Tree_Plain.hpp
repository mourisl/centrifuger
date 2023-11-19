#ifndef _MOURISL_COMPACTDS_TREE_PLAIN
#define _MOURISL_COMPACTDS_TREE_PLAIN

#include <vector>

#include "Tree.hpp"
#include "SimpleVector.hpp"

namespace compactds {
struct _plainTreeNode
{
  size_t parent ;
  size_t sibling ;
  size_t child ;
  size_t lastChild ;

  size_t label ; // the label from the parent to itself

  _plainTreeNode(size_t p, size_t s, size_t c, size_t lc)
  {
    parent = p ;
    sibling = s ;
    child = c ;
    lastChild = lc ;
    label = 0 ;
  }

  void Save(FILE *fp) 
  {
    SAVE_VAR(fp, parent) ;
    SAVE_VAR(fp, sibling) ;
    SAVE_VAR(fp, child) ;
    SAVE_VAR(fp, lastChild) ;
  }

  void Load(FILE *fp)
  {
    LOAD_VAR(fp, parent) ;
    LOAD_VAR(fp, sibling) ;
    LOAD_VAR(fp, child) ;
    LOAD_VAR(fp, lastChild) ;
  }
} ;

class Tree_Plain: public Tree 
{
private:
  std::vector<struct _plainTreeNode> _nodes ;
public:
  Tree_Plain() {}
  ~Tree_Plain() {}

  void Init()
  {
    _n = 1 ;
    
    struct _plainTreeNode node(0, 0, 0, 0) ;
    _nodes.push_back(node) ;
  }

  size_t GetSpace(bool inclusive = true)
  {
    return _nodes.capacity() * sizeof(struct _plainTreeNode) + (inclusive ? sizeof(*this) : 0) ;
  }

  // Assumes parent is already in the tree
  //@return: tree index
  size_t AddNode(size_t parent)
  {
    size_t id = _nodes.size() ; 
    struct _plainTreeNode node(parent, 0, 0, 0) ;
    size_t lastSibling = LastChild(parent) ;

    if (lastSibling == 0)
      _nodes[parent].child = id ;
    else
      _nodes[lastSibling].sibling = id ;
    
    _nodes[parent].lastChild = id ;
    _nodes.push_back(node) ;
    ++_n ;

    return id ;
  }

  size_t Root() const 
  {
    return 0 ;
  }

  // t-th(1-based) child ;
  size_t ChildSelect(size_t v, size_t t) const
  {
    --t ;
    
    size_t c = FirstChild(v) ;
    size_t i ;
    for (i = 0 ; i < t ; ++i)
      c = NextSibling(c) ;
    return c ;
  }

  size_t FirstChild(size_t v) const 
  {
    return _nodes[v].child ; 
  }

  size_t LastChild(size_t v) const 
  {
    return _nodes[v].lastChild ;
  }

  size_t ChildrenCount(size_t v) const
  {
    size_t i ;
    size_t c = FirstChild(v) ;
    for (i = 0 ; c != 0 ; ++i)
      c = NextSibling(c) ;
    return i ;
  }
  
  // return: v is the ret-th (1-based) child of the parent.
  size_t ChildRank(size_t v) const
  {
    if (v == Root())
      return 0 ;
    size_t c = FirstChild( Parent(v) ) ; 
    size_t i ;
    for (i = 0 ; c != v ; ++i)
      c = NextSibling(c) ;
    return i + 1 ; // +1: rank is always 1-based
  }

  size_t NextSibling(size_t v) const 
  {
    return _nodes[v].sibling ;
  }

	size_t PrevSibling(size_t v) const
	{
		size_t i ;
		size_t c = FirstChild( Parent(v) ) ;
		for (i = 0 ; v != NextSibling(c) ; ++i)
			c = NextSibling(c) ;
		return c ;
	}

  size_t Parent(size_t v) const
  {
    return _nodes[v].parent ;
  }

  bool IsLeaf(size_t v) const
  {
    if (_nodes[v].child == 0)
      return true ;
    return false ;
  }

  size_t LCA(size_t u, size_t v) const
  {
    SimpleVector<size_t> upath ; 
    SimpleVector<size_t> vpath ;

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

  size_t NodeMap(size_t v) const 
  {
    return v ;
  }

  size_t NodeSelect(size_t i) const 
  {
    return i ;
  }

  const std::vector<struct _plainTreeNode>& GetTreeData() const
  {
    return _nodes ; 
  }

  void SetLabel(size_t v, size_t l)
  {
    _nodes[v].label = l ;
  }

  // Number of children with label l. 1: has such children. 0-don't 
  size_t ChildrenLabeled(size_t v, size_t l) const  
  {
    size_t ret = 0 ;
    int c = _nodes[v].child ;
    while (c != 0)
    {
      if (_nodes[c].label == l)
        ++ret ;
      c = _nodes[c].sibling ;
    }
    return ret ;
  }

  // The child with label l.
  size_t LabeledChildSelect(size_t v, size_t l, size_t t) const 
  {
    size_t cnt = 0 ;
    size_t c = _nodes[v].child ;
    while (c != 0)
    {
      if (_nodes[c].label == l)
        ++cnt ;
      if (cnt >= t)
        break ;
      c = _nodes[c].sibling ;
    }
    return c ;
  }

  // The label of the edge that leads to node v.
  size_t ChildLabel(size_t v) const 
  {
    return _nodes[v].label ;  
  }

  void Save(FILE *fp) 
  {
    Tree::Save(fp) ;
    size_t i ;
    for (i = 0 ; i < _n ; ++i)
      _nodes[i].Save(fp) ;
  }

  void Load(FILE *fp)
  {
    std::vector< struct _plainTreeNode >().swap(_nodes) ;    

    Tree::Load(fp) ;
    size_t i ;
    struct _plainTreeNode node(0, 0, 0, 0) ;
    for (i = 0 ; i < _n ; ++i)
    {
      node.Load(fp) ;
      _nodes.push_back(node) ;
    }
  }
} ;

}

#endif

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
  size_t _root ; // allow flexible root setting. For compact representation, this 
public:
  Tree_Plain() 
  {
    _root = 0 ;
  }

  ~Tree_Plain() {}

  void SetRoot(size_t r)
  {
    _root = r ;
  }

  void Init()
  {
    _n = 1 ;
    
    struct _plainTreeNode node(_root, _root, _root, _root) ;
    _nodes.push_back(node) ;
  }

  void Init(size_t n)
  {
    _nodes.clear() ;

    size_t i ;
    _n = n ;
    for (i = 0 ; i < _n ; ++i)
    {
      struct _plainTreeNode node(_root, _root, _root, _root) ;
      _nodes.push_back(node) ;
    }
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
    struct _plainTreeNode node(parent, _root, _root, _root) ;
    size_t lastSibling = LastChild(parent) ;

    if (lastSibling == _root)
      _nodes[parent].child = id ;
    else
      _nodes[lastSibling].sibling = id ;
    
    _nodes[parent].lastChild = id ;
    _nodes.push_back(node) ;
    ++_n ;

    return id ;
  }

  void AddEdge(size_t c, size_t parent)
  {
    _nodes[c].parent = parent ;
    
    size_t lastSibling = LastChild(parent) ;
    if (lastSibling == _root)
      _nodes[parent].child = c ;
    else
      _nodes[lastSibling].sibling = c ;
    
    _nodes[parent].lastChild = c ;
  }

  size_t Root() const 
  {
    return _root ;
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
    for (i = 0 ; c != _root ; ++i)
      c = NextSibling(c) ;
    return i ;
  }

  std::vector<size_t> GetChildren(size_t v) const
  {
    std::vector<size_t> ret ;
    size_t c = _nodes[v].child ;
    while (c != _root)
    {
      ret.push_back(c) ;
      c = _nodes[c].sibling ;
    }
    return ret ;
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
    if (_nodes[v].child == _root)
      return true ;
    return false ;
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
    size_t c = _nodes[v].child ;
    while (c != _root)
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
    while (c != _root)
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
    SAVE_VAR(fp, _root) ;
    for (i = 0 ; i < _n ; ++i)
      _nodes[i].Save(fp) ;
  }

  void Load(FILE *fp)
  {
    std::vector< struct _plainTreeNode >().swap(_nodes) ;    

    Tree::Load(fp) ;
    LOAD_VAR(fp, _root) ;
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

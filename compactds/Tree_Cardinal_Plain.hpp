// Cardinal tree, plain representation
//
#ifndef _MOURISL_COMPACTDS_TREE_CARDINAL_PLAIN
#define _MOURISL_COMPACTDS_TREE_CARDINAL_PLAIN 

#include "Tree_Cardinal.hpp"

namespace compactds {
struct _plainCardinalTreeNode
{
  size_t k ; // It is the k-th(0-based) children of the parent
  size_t parent ;
  size_t *children ;

  _plainCardinalTreeNode()
  {
    children = NULL ;
  }

  ~_plainCardinalTreeNode() {} ;

  _plainCardinalTreeNode(size_t p, size_t inK, size_t c)
  {
    parent = p ;
    k = inK ;
    children = (size_t *)calloc(c, sizeof(size_t)) ;
  }

  void Free()
  {
    if (children)
    {
      free(children) ;
      children = NULL ;
    }
  }

  void Save(FILE *fp, size_t c)
  {
    SAVE_VAR(fp, parent) ;
    SAVE_ARR(fp, children, c) ;
  }

  void Load(FILE *fp, size_t c)
  {
    Free() ;

    LOAD_VAR(fp, parent) ;
    children = (size_t *)calloc(c, sizeof(size_t)) ;
    LOAD_ARR(fp, children, c) ;
  }
} ;

class Tree_Cardinal_Plain : public Tree_Cardinal
{
private:
  std::vector<struct _plainCardinalTreeNode> _nodes ;
  size_t _c ; //child count
public:
  Tree_Cardinal_Plain() {} ;
  ~Tree_Cardinal_Plain() 
  {
    Free() ;
  } 

  void Init(size_t childCount)
  {
    _c = childCount ;
    struct _plainCardinalTreeNode node(0, 0, _c) ;
    _nodes.push_back(node) ;
    _n = 1 ;
  }

  void Free()
  {
    size_t i ;
    for (i = 0 ; i < _n ; ++i)
      _nodes[i].Free() ;
    _n = 0 ;
  }
  
  size_t GetSpace(bool inclusive = true)
  {
    return _nodes.capacity() * sizeof(_nodes[0]) + (inclusive ? sizeof(*this) : 0) ;
  }

  size_t AddNode(size_t parent, size_t k) 
  {
    size_t id = _nodes.size() ; 
    struct _plainCardinalTreeNode node(parent, k, _c) ;
    
    _nodes[parent].children[k] = id ;
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
    size_t i ;
    size_t cnt = 0 ;
    for (i = 0 ; i < _c ; ++i)
    {
      if (_nodes[v].children[i] != 0)
      {
        if (cnt == t - 1)
          return _nodes[v].children[i] ;
        ++cnt ;
      }
    }
    return 0 ;
  }
  
  size_t FirstChild(size_t v) const 
  {
    size_t i ;
    for (i = 0 ; i < _c ; ++i)
      if (_nodes[v].children[i] != 0)
        return _nodes[v].children[i] ;
    return 0 ;
  }

  size_t LastChild(size_t v) const 
  {
    size_t i ; 
    for (i = _c - 1 ; i < _c ; --i)
      if (_nodes[v].children[i] != 0) 
        return _nodes[v].children[i] ;
    return 0 ;
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
    size_t p = Parent(v) ;
    size_t ret = 0 ;
    size_t i ;
    
    // inclusive, for rank is always 1-based
    for (i = 0 ; i <= _nodes[v].k ; ++i)
      if (_nodes[p].children[i] != 0)
        ++ret ;
    return ret ; 
  }

  size_t NextSibling(size_t v) const 
  {
    size_t i ;
    size_t p = Parent(v) ;
    for (i = _nodes[v].k + 1 ; i < _c ; ++i)
      if (_nodes[p].children[i] != 0)
        return _nodes[p].children[i] ;
    return 0 ;
  }

  size_t PrevSibling(size_t v) const
  {
    size_t i ;
    size_t p = Parent(v) ;
    for (i = _nodes[v].k - 1 ; i < _c ; --i)
      if (_nodes[p].children[i] != 0)
        return _nodes[p].children[i] ;
    return 0 ;
  }

  size_t Parent(size_t v) const 
  {
    return _nodes[v].parent ;
  }

  bool IsLeaf(size_t v) const 
  {
    return (ChildrenCount(v) == 0) ;
  }

  size_t NodeMap(size_t v) const 
  {
    return v ;
  }

  size_t NodeSelect(size_t i) const 
  {
    return i ;
  }
  
  // Number of children with label l. 1: has such children. 0-don't 
  size_t ChildrenLabeled(size_t v, size_t l) const
  {
    if (_nodes[v].children[l] != 0)
      return 1 ;
    else
      return 0 ;
  }

  // The children with label l.
  size_t LabeledChild(size_t v, size_t l) const
  {
    return _nodes[v].children[l] ;
  }

  // The label of the edge that leads to node v.
  size_t ChildLabel(size_t v) const 
  {
    return _nodes[v].k ;
  }

  const std::vector<struct _plainCardinalTreeNode>& GetTreeData() const
  {
    return _nodes ; 
  }

  void Save(FILE *fp)
  {
    Tree_Cardinal::Save(fp) ;
    size_t i ;
    for (i = 0 ; i < _n ; ++i)
      _nodes[i].Save(fp, _c) ;
  }

  void Load(FILE *fp)
  {
    std::vector< struct _plainCardinalTreeNode >().swap(_nodes) ;    

    Tree_Cardinal::Load(fp) ;
    size_t i ;
    struct _plainCardinalTreeNode node ;
    for (i = 0 ; i < _n ; ++i)
    {
      node.Load(fp, _c) ;
      _nodes.push_back(node) ;
    }
  }

} ;
}

#endif

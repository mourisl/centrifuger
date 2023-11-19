#ifndef _MOURISL_COMPACTDS_TREE_DFUDS
#define _MOURISL_COMPACTDS_TREE_DFUDS

// Depth-First Unary Degree Sequence (Chapter 8.3)
// The implementation details might be different as we are using 0-index 
//    tree node id (i) as well.
// In the implementation, v is for the index on the encoded bitvector B.

#include "Tree.hpp"
#include "Tree_Plain.hpp"
#include "Tree_Cardinal_Plain.hpp"
#include "Bitvector_Plain.hpp"

namespace compactds {
class Tree_DFUDS: public Tree
{
private:
  Bitvector_Plain _B ;
  DS_Parenthesis _bp ; // dangling structure, the parentehesis is almost balanced.
  
  size_t _m ; //|_B|, just for coding simplicity
  
  void Build(const struct _plainTreeNode *treeNodes, size_t n, size_t tag, size_t *treeIdMap, size_t &visited, size_t &bi)
  {
    size_t c ;
    c = treeNodes[tag].child ;
    
    treeIdMap[tag] = visited ;
    ++visited ;

    for (c = treeNodes[tag].child ; c != 0 ; c = treeNodes[c].sibling) 
    {
      _B.BitSet(bi) ;
      ++bi ;
    }
    ++bi ; // set as 0

    c = treeNodes[tag].child ;

    for (c = treeNodes[tag].child ; c != 0 ; c = treeNodes[c].sibling) 
      Build(treeNodes, n, c, treeIdMap, visited, bi) ;  
  }

  void BuildFromCardinalTree(const struct _plainCardinalTreeNode *treeNodes, size_t n, size_t childCnt, size_t tag, size_t *treeIdMap, size_t &visited, size_t &bi)
  {
    size_t i ;
    treeIdMap[tag] = visited ;
    ++visited ;

    for (i = 0 ; i < childCnt; ++i) 
    {
      if (treeNodes[tag].children[i] == 0)
        continue ;
      _B.BitSet(bi) ;
      ++bi ;
    }
    ++bi ;

    for (i = 0 ; i < childCnt; ++i) 
    {
      if (treeNodes[tag].children[i] == 0)
        continue ;
      BuildFromCardinalTree(treeNodes, n, childCnt, treeNodes[tag].children[i], treeIdMap, visited, bi) ; 
    }
  }
public:
  Tree_DFUDS() {} 
  ~Tree_DFUDS()
  {
    Free() ;
  }

  void Free()
  {
    if (_n > 0)
    {
      _B.Free() ;
      _bp.Free() ;
      _n = 0 ;
      _m = 0 ;
    }
  }

  size_t GetSpace(bool inclusive = true) 
  {
    return _space + (inclusive ? sizeof(*this) : 0) ;
  }

  void Init(const struct _plainTreeNode *treeNodes, size_t n, size_t *treeIdMap)
  {
    _n = n ;
    _B.Malloc(2 * _n - 1) ;
    size_t bi = 0 ;
    size_t visited = 0 ;
    Build(treeNodes, n, 0, treeIdMap, visited, bi) ;
    _m = bi ;

    _B.Init() ;
    _bp.Init(_B.GetData(), _m, 0, 2) ;
    
    _space = _B.GetSpace() - sizeof(_B) + _bp.GetSpace(false) ;
  }

  void InitFromCardinalTree(const struct _plainCardinalTreeNode *treeNodes, size_t n, 
      size_t childCount, size_t *treeIdMap)
  {
    _n = n ;
    _B.Malloc(2 * _n - 1) ;
    size_t bi = 0 ;
    size_t visited = 0 ;
    BuildFromCardinalTree(treeNodes, n, childCount, 0, treeIdMap, visited, bi) ;
    _m = bi ;

    _B.Init() ;
    _bp.Init(_B.GetData(), _m, 0, 2) ;
    
    _space = _B.GetSpace() - sizeof(_B) + _bp.GetSpace(false) ;
  }

  // The index in B
  size_t Root() const
  {
    return 0 ;
  }

  // @return: the t-th child (1-based) of node v in B vector
  size_t ChildSelect(size_t v, size_t t) const
  {
    size_t childCnt = ChildrenCount(v) ;
    return _bp.Close( v + childCnt - t, _B.GetData(), _m) + 1 ;
  }

  size_t FirstChild(size_t v) const 
  {
    return _B.Succ0(v) + 1 ;
  }

  size_t LastChild(size_t v) const 
  {
    return _bp.Close(v, _B.GetData(), _m) + 1 ;
  }
  
  size_t ChildrenCount(size_t v) const
  {
    return _B.Succ0(v) - v ;
  }
  
  // return: v is the ret-th (1-based) child of the parent.
  size_t ChildRank(size_t v) const
  {
    size_t open = _bp.Open(v - 1, _B.GetData(), _m) ;
    return  _B.Succ0(open) - open ;
  }
  
  // The silbing function assumes v has
  //   those siblings.
  size_t NextSibling(size_t v) const
  {
    return _bp.GetRmmTree().FwdSearch(v, -1, _B.GetData(), _m) + 1 ; 
  }

  size_t PrevSibling(size_t v) const
  {
    // Notice the order here, the closer to the end of the (((), the close
    //   ) corresponds to the earlier children
    return _bp.Close(_bp.Open(v - 1, _B.GetData(), _m) + 1,  
        _B.GetData(), _m) + 1 ;
  }
  
  size_t Parent(size_t v) const
  {
    if (v == 0)
      return 0 ;

    return _B.Pred0(_bp.Open(v - 1, _B.GetData(), _m)) + 1 ;
  }

  // # of nodes in the substree, inclusive. 
  size_t SubTreeSize(size_t v) const
  {
    // For 2*m-1 parenthesis in the range: m ')', m-1 '('
    // The equation below is actually (2*m-2)/2 + 1
    return (_bp.GetRmmTree().FwdSearch(v, -1, _B.GetData(), _m) - v) / 2 + 1 ;
  }
  
  // Whether u is an ancestor of v.
  bool IsAncestor(size_t u, size_t v) const
  {
    size_t uend = _bp.GetRmmTree().FwdSearch(u, -1, _B.GetData(), _m) ;
    if (v >= u && v <= uend)
      return true ;
    return false ;
  }

  bool IsLeaf(size_t v) const 
  {
    return (_B.Access(v) == 0) ;
  }

  size_t LCA(size_t u, size_t v) const
  {
    if (v < u)
    {
      size_t tmp = v ;
      v = u ;
      u = tmp ;
    }

    if (IsAncestor(u, v))
      return u ;
    //printf("%s: %d %d. %d\n", __func__, u, v, 
    //    _bp.GetRmmTree().Rmq(u, v - 1, _B.GetData(), _m)) ;
    
    // Think about this more.
    // Example (())): node with two leaves 
    // v-1 in Rmq then add 1 back handles both the leaf case and internal node case
    return Parent( _bp.GetRmmTree().Rmq(u, v - 1,
          _B.GetData(), _m) + 1)  ;
  }

  size_t LeafCountInSubTree(size_t v) const
  {
    if (IsLeaf(v))
      return 1 ;
    else
    {
      size_t vend = _bp.GetRmmTree().FwdSearch(v, -1, _B.GetData(), _m) ;
      return _bp.PatternRank(vend - 1, _B.GetData(), _m) - _bp.PatternRank(v, _B.GetData(), _m)  ;
    }
  }

  // Rank and select with respect to leaves in _B order
  // Assuming v is the leaf
  size_t LeafRank(size_t v, int inclusive = 1) const
  {
    // The end of 00 corresponds to the leaf
    return _bp.PatternRank(v - 1, _B.GetData(), _m, inclusive) ;
  }

  size_t LeafSelect(size_t i) const
  {
    return _bp.PatternSelect(i, _B.GetData(), _m) + 1 ;
  }
  
  // Maps index in B (v) back up to the actual node id
  size_t NodeMap(size_t v) const
  {
    // Inclusive==0 because leaf node will have ')' on the index
    return _B.Rank(0, v, /*inclusive=*/0) ; 
  }

  // Map actual node id to index in B (v).
  size_t NodeSelect(size_t i) const 
  {
    if (i == 0)
      return 0 ;
    return _B.Select(0, i) + 1 ;
  }

  void Save(FILE *fp)
  {
    Tree::Save(fp) ;

    _B.Save(fp) ;
    _bp.Save(fp) ;
  }

  void Load(FILE *fp)
  {
    Free() ;

    Tree::Load(fp) ;

    _B.Load(fp) ;
    _bp.Load(fp) ;
    
    _m = 2 * _n - 1 ;
  }
} ;

} // end of name space

#endif

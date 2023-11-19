#ifndef _MOURISL_COMPACTDS_TREE_BP
#define _MOURISL_COMPACTDS_TREE_BP

// Represent a tree by balanced parenthesis (Chapter 8.2)
// Each v points to a parenthesis like (...) containing the substree 
// The implementation details might be different as we are using 0-index 
//    tree node id (i) as well.
// In the implementation, v is for the index on the encoded bitvector B.

#include "Tree.hpp"
#include "Tree_Plain.hpp"
#include "Tree_Cardinal_Plain.hpp"
#include "Bitvector_Plain.hpp"
#include "DS_Parenthesis.hpp"

namespace compactds {
class Tree_BP: public Tree
{
private:
  Bitvector_Plain _B ; // bits representation of the parenthesis 
  DS_Parenthesis _bp ; // dangling structure
 
  // DFS to mark the parenthesis as B 
  // tag: tree node id. bi: index on B
  void Build(const struct _plainTreeNode *treeNodes, size_t n, size_t tag, 
      size_t *treeIdMap, size_t &visited, size_t &bi)
  {
    size_t c ;
    _B.BitSet(bi) ;
    
    treeIdMap[tag] = visited ;
    ++visited ;

    ++bi ;
    for (c = treeNodes[tag].child ; c != 0 ; c = treeNodes[c].sibling) 
      Build(treeNodes, n, c, treeIdMap, visited, bi) ;
    //_B.BitClear(bi) ; // close the parentehsis
    ++bi ;
  }
  
  void BuildFromCardinalTree(const struct _plainCardinalTreeNode *treeNodes, size_t n, size_t childCnt, size_t tag, size_t *treeIdMap, size_t &visited, size_t &bi)
  {
    size_t i ;
    _B.BitSet(bi) ;

    treeIdMap[tag] = visited ;
    ++visited ;

    ++bi ;
    for (i = 0 ; i < childCnt ; ++i)
    {
      size_t c = treeNodes[tag].children[i] ;
      if (c == 0)
        continue ;
      BuildFromCardinalTree(treeNodes, n, childCnt, c, treeIdMap, visited, bi) ;
    }
    //_B.BitClear(bi) ; // close the parentehsis
    ++bi ;
  }
public:
  Tree_BP() {} 
  ~Tree_BP()
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
    }
  }

  size_t GetSpace(bool inclusive = true) 
  {
    return _space + (inclusive ? sizeof(*this) : 0) ;
  }

  void Init(const struct _plainTreeNode *treeNodes, size_t n, size_t *treeIdMap)
  {
    _n = n ;
    _B.Malloc(2 * n) ; 
    
    size_t bi = 0 ;
    size_t visited = 0 ;
    Build(treeNodes, n, 0, treeIdMap, visited, bi) ;
    
    _B.Init() ;
    
    // the last 2,2 is for pattern 10 as "()"
    // Note that due to our reading the bits from low to high, "1" will be the first bit
    //  In other word, order is reversed.
    _bp.Init(_B.GetData(), 2 * _n, 1, 2) ; 
    
    _space = _B.GetSpace() - sizeof(_B) + _bp.GetSpace(false) ;
  }

  void InitFromCardinalTree(const struct _plainCardinalTreeNode *treeNodes, size_t n, size_t childCount, size_t *treeIdMap)
  {
    _n = n ;
    _B.Malloc(2 * n) ; 

    size_t bi = 0 ;
    size_t visited = 0 ;
    BuildFromCardinalTree(treeNodes, n, childCount, 0, treeIdMap, visited, bi) ;

    _B.Init() ;

    // the last 2,2 is for pattern 10 as "()"
    // Note that due to our reading the bits from low to high, "1" will be the first bit
    //  In other word, order is reversed.
    _bp.Init(_B.GetData(), 2 * _n, 1, 2) ; 

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
    return _bp.Open(_bp.GetRmmTree().MinSelect(v + 1, _bp.Close(v, _B.GetData(), 2 * _n) - 1,
        t, _B.GetData(), 2 * _n), _B.GetData(), 2*_n) ; 
  } 

  size_t FirstChild(size_t v) const
  {
    return v + 1 ;
  }

  // Parenthesis like
  // (...(...))
  // |   |
  // v   lastchild
  size_t LastChild(size_t v) const
  {
    return _bp.Open(_bp.Close(v, _B.GetData(), 2*_n)-1, _B.GetData(), 2*_n) ;
  }
  
  size_t ChildrenCount(size_t v) const
  {
    if (IsLeaf(v))
      return 0 ;
    // Each child's (...) has excess 0 after the end.
    return _bp.GetRmmTree().MinCount(v + 1, _bp.Close(v, _B.GetData(), 2*_n) - 1, 
        _B.GetData(), 2*_n) ;
  }
  
  // return: v is the ret-th (1-based) child of the parent.
  size_t ChildRank(size_t v) const
  {
    if (v == Root())
      return 0 ;
    size_t p = Parent(v) ;
    if (p + 1 == v)
      return 1 ;
    return _bp.GetRmmTree().MinCount(p + 1, v - 1, _B.GetData(), 2*_n) + 1 ;
  }
  
  // The silbing function assumes v has
  //   those siblings.
  size_t NextSibling(size_t v) const
  {
    return _bp.Close(v, _B.GetData(), 2*_n) + 1 ;
  }

  size_t PrevSibling(size_t v) const
  {
    return _bp.Open(v - 1, _B.GetData(), 2*_n) ;
  }
  
  size_t Parent(size_t v) const
  {
    if (v == Root())
      return 0 ;
    return _bp.Enclose(v, _B.GetData(), 2*_n) ;
  }

  bool IsLeaf(size_t v) const 
  {
    if (_B.Access(v + 1) == 0)
      return true ;
    return false ;
  }

  size_t LCA(size_t u, size_t v) const
  {
    if (u > v)
    {
      size_t tmp = u ;
      u = v ;
      v = tmp ;
    }

    //printf("%d %d: %d %d\n", u, v, _bp.GetRmmTree().Rmq(u, v, _B.GetData(), 2*_n) + 1,
    //    _bp.Enclose(_bp.GetRmmTree().Rmq(u, v, _B.GetData(), 2*_n) + 1, _B.GetData(), 2*_n)) ;
    return _bp.Enclose(_bp.GetRmmTree().Rmq(u, v, _B.GetData(), 2*_n) + 1, _B.GetData(), 2*_n) ;
  }
  
  // Maps index in B (v) back up to the actual node id
  // Pre-order
  size_t NodeMap(size_t v) const
  {
    return _B.Rank(1, v, /*inclusive=*/0) ; 
  }

  //Map actual node id to index in B (v).
  // Pre-order Select
  size_t NodeSelect(size_t i) const 
  {
    return _B.Select(1, i + 1) ; 
  }

  size_t PostOrder(size_t v) const
  {
    return _B.Rank(0, _bp.Close(v, _B.GetData(), 2 * _n), /*inclusive*/0) ;
  }

  size_t PostOrderSelect(size_t i) const
  {
    return _bp.Open(_B.Select(0, i + 1), _B.GetData(), 2 * _n) ;   
  }

  // Root has depth 0
  size_t Depth(size_t v) const
  {
    // Kind of excess
    // inclusive=0 means rank-1 here
    return 2 * _B.Rank(1, v, /*inclusive=*/0) - v ;
  }

  // # of nodes in the substree, inclusive. 
  size_t SubTreeSize(size_t v) const
  {
    return (_bp.Close(v, _B.GetData(), 2 * _n) - v + 1) / 2 ;
  }

  // Whether u is an ancestor of v.
  bool IsAncestor(size_t u, size_t v) const
  {
    size_t uclose = _bp.Close(u, _B.GetData(), 2 * _n) ;
    if (u <= v && v <= uclose)
      return true ;
    return false ;
  }

  // The ancestor at d levels above
  size_t LevelAncestor(size_t v, int64_t d) const
  {
    return _bp.GetRmmTree().BwdSearch(v, -d, _B.GetData(), 2 * _n) ;
  }

  // Would be v it self if v is the leaf.
  size_t DeepestNode(size_t v) const
  {
    return _bp.GetRmmTree().RMq(v, _bp.Close(v, _B.GetData(), 2 * _n), _B.GetData(), 2 * _n) ;
  }

  // The distance from v to the deepest leaf
  // 0 if v is the leaf
  size_t Height(size_t v) const
  {
    size_t depthv = Depth(v) ;
    size_t depthc = Depth( DeepestNode(v) ) ;
    return depthc - depthv ;
  }

  // Number of leaves in the subtree of v
  size_t LeafCountInSubTree(size_t v) const
  {
    if (IsLeaf(v))
      return 1 ;
    else
    {
      // Since close(v) is ")", the LeafRank(close(v)) is automatically exclusive
      return LeafRank( _bp.Close(v, _B.GetData(), 2 *_n)) - LeafRank(v) ;
    }
  }

  // Rank and select with respect to leaves in _B order
  size_t LeafRank(size_t v, int inclusive = 1) const
  {
    return _bp.PatternRank(v, _B.GetData(), 2*_n, inclusive) ;
  }

  size_t LeafSelect(size_t i) const
  {
    return _bp.PatternSelect(i, _B.GetData(), 2*_n) ;
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
  }
} ;

} // end of name space

#endif

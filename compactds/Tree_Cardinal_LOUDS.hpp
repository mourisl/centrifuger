#ifndef _MOURISL_COMPACTDS_TREE_CARDINAL_LOUDS
#define _MOURISL_COMPACTDS_TREE_CARDINAL_LOUDS

// Level-order Unary degree sequence for cardinal tree (Chapter 8.1.1, Algorithm 8.4)
// The implementation details might be different as we are using 0-index 
//    tree node id (i) as well.
// In the implementation, v is for the index on the encoded bitvector B 
//  to be consistent with the genral LOUDS represent. This is DIFFERENT
//  from the textbook.

#include "Tree_Cardinal.hpp"
#include "Tree_Cardinal_Plain.hpp"
#include "Bitvector_Plain.hpp"

namespace compactds {
template <class BvClass = Bitvector_Plain>
class Tree_Cardinal_LOUDS: public Tree_Cardinal
{
private:
  BvClass _B ;
public:
  size_t GetSpace(bool inclusive = true) 
  {
    return _space + (inclusive ? sizeof(*this) : 0) ;
  }
  
  // Use c bits per node.
  void Init(const struct _plainCardinalTreeNode *treeNodes, size_t n, size_t c, size_t *treeIdMap)
  {
    size_t i ;
    WORD *W = Utils::MallocByBits(c * n) ;

    _n = n ;
    _c = c ;
    
    // BFS on tree nodes
    // The algorithm 8.3 will change the original plain tree
    size_t *queue = (size_t *)malloc(sizeof(size_t) * n) ;
    size_t qhead, qtail ;
    queue[0] = 0 ;
    qhead = 0 ;
    qtail = 1 ;
    while (qhead < qtail)
    {
      size_t node = queue[qhead] ;
      if (treeIdMap != NULL)
        treeIdMap[node] = qhead ;

      ++qhead ;
      
      for (i = 0 ; i < _c ; ++i)
      {
        if (treeNodes[node].children[i] != 0)
        {
          Utils::BitSet(W, _c * (qhead - 1) + i) ;
          queue[qtail] = treeNodes[node].children[i] ;
          ++qtail ;
        }
      }
    }
    free(queue) ;

    _B.Init(W, c * n) ;
    _space = _B.GetSpace() - sizeof(_B);
    free(W) ;
  }

  // The index in B
  size_t Root() const
  {
    return 0 ;
  }

  // @return: the t-th child (0-based) of node v in B vector
  size_t ChildSelect(size_t v, size_t t) const
  {
    return (_B.Rank(1, v, /*inclusive=*/0) + t) * _c  ; 
  }

  size_t FirstChild(size_t v) const 
  {
    return (_B.Rank(1, v, /*inclusive=*/0) + 1) * _c ; 
  }

  size_t LastChild(size_t v) const 
  {
    return (_B.Rank(1, v + _c - 1)) * _c ;
  }
  
  size_t ChildrenCount(size_t v) const
  {
    return _B.Rank(1, v + _c - 1) - _B.Rank(1, v, 0) ;
  }
  
  // return: v is the ret-th (1-based) child of the parent.
  size_t ChildRank(size_t v) const
  {
    if ( v == Root())
      return 0 ;
    size_t tid = NodeMap(v) ;
    size_t j = _B.Select(1, tid) ; // edge from parent to v
    return _B.Rank(1, j) - _B.Rank(1, j/_c * _c, /*inclusive=*/0) ;
  }  
  
  // The silbing function assumes v has
  //   those siblings.
  size_t NextSibling(size_t v) const
  {
    return v + 1 ;
  }

  size_t PrevSibling(size_t v) const
  {
    return v - 1 ; 
  }
  
  size_t Parent(size_t v) const
  {
    if (v == 0)
      return Root() ;
    else
    {
     size_t tid = NodeMap(v) ;
     // _B.Select(1, tid) identify the edge from the parent to v
     // Notice that even though tree node id is 0-based, the edge starts
     //   to correspnds to node id 1
     return (_B.Select(1, tid) / _c) * _c ;
    } 
  }

  bool IsLeaf(size_t v) const 
  {
    return ChildrenCount(v) == 0 ;
  }

  size_t LCA(size_t u, size_t v) const
  {
    while (u != v)
    {
      if (u > v)
        u = Parent(u) ;
      else
        v = Parent(v) ;
    }
    return u ;
  }
  
  // Maps index in B (v) back up to the actual node id
  size_t NodeMap(size_t v) const
  {
    return v / _c ;
  }

  //Map actual node id to index in B (v).
  size_t NodeSelect(size_t i) const 
  {
    if (i == 0)
      return Root() ;
    else
      return i * _c ;
  }

  // Number of children with label l. 1: has such children. 0-don't 
  size_t ChildrenLabeled(size_t v, size_t l) const 
  {
    return (_B.Access(v + l) == 1) ? 1 : 0 ;
  }

  // The childr with label l.
  // Assuming label l's child exist
  // Notice the difference from Child()
  size_t LabeledChild(size_t v, size_t l) const 
  {
    return (_B.Rank(1, v + l, /*inclusive=*/1)) * _c  ; 
  }

  // The label of the edge that leads to node v.
  size_t ChildLabel(size_t v) const 
  {
    if (v == Root())
      return 0 ;
    size_t tid = NodeMap(v) ;
    size_t j = _B.Select(1, tid) ; // edge from parent to v
    return j % _c ;
  }

  void Save(FILE *fp)
  {
    Tree_Cardinal::Save(fp) ;
    _B.Save(fp) ;
  }

  void Load(FILE *fp)
  {
    Tree_Cardinal::Load(fp) ;
    _B.Load(fp) ;
  }
 
} ;

} // end of name space

#endif

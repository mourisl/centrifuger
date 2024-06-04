#ifndef _MOURISL_COMPACTDS_TREE_LOUDS
#define _MOURISL_COMPACTDS_TREE_LOUDS

// Level-order Unary degree sequence (Chapter 8.1)
// The implementation details might be different as we are using 0-index 
//    tree node id (i) as well.
// In the implementation, v is for the index on the encoded bitvector B.

#include "Tree.hpp"
#include "Tree_Plain.hpp"
#include "Bitvector_Plain.hpp"

namespace compactds {
class Tree_LOUDS: public Tree
{
private:
  Bitvector_Plain _B ;
public:
  size_t GetSpace(bool inclusive = true) 
  {
    return _space + (inclusive ? sizeof(*this) : 0) ;
  }

  void Init(const struct _plainTreeNode *treeNodes, size_t n, size_t *treeIdMap)
  {
    size_t i ;
    _B.Malloc(2 * n - 1) ; // n nodes, n-1 edges
		_space = _B.GetSpace() - sizeof(_B); 
  
    // BFS on tree nodes
    // The algorithm 8.3 will change the original plain tree
    size_t *queue = (size_t *)malloc(sizeof(size_t) * n) ;
    size_t qhead, qtail ;
    size_t m = 0 ; // position on _B
    queue[0] = 0 ;
    qhead = 0 ;
    qtail = 1 ;
    while (qhead < qtail)
    {
      size_t node = queue[qhead] ;
      if (treeIdMap != NULL)
        treeIdMap[node] = qhead ;

      ++qhead ;
      
      size_t childCnt = 0 ;
      size_t c = treeNodes[node].child ;
      for (childCnt = 0 ; c != 0 ; ++childCnt)
      {
        queue[qtail] = c ;
        ++qtail ;
        c = treeNodes[c].sibling ;
      }
      
      for (i = m ; i < m + childCnt ; ++i)
        _B.BitSet(i) ;
      m += childCnt + 1 ;
    }
    free(queue) ;

    _B.Init() ;
  }

  // The index in B
  size_t Root() const
  {
    return 0 ;
  }

  // @return: the t-th child (1-based) of node v in B vector
  size_t ChildSelect(size_t v, size_t t) const
  {
    return _B.Select(0, _B.Rank(1, v + t - 1)) + 1 ;
  }

  size_t FirstChild(size_t v) const 
  {
    return ChildSelect(v, 1) ;
  }

  size_t LastChild(size_t v) const 
	{
		return ChildSelect(v, ChildrenCount(v)) ;		
	}
  
  size_t ChildrenCount(size_t v) const
  {
		return _B.Succ0(v) - v ;
  }
  
  // return: v is the ret-th (1-based) child of the parent.
  size_t ChildRank(size_t v) const
  {
    if (v == Root())
      return 0 ;
		//_B.Rank(0, v-1): the node id
		//_B.Select(1, _B.Rank(0, v-1)): The edge connect the parent and v
    size_t j = _B.Select(1, _B.Rank(0, v - 1)) ;
		return j - _B.Pred0(j) ; // rank is always 1-based
  }
  
  // The silbing function assumes v has
  //   those siblings.
  size_t NextSibling(size_t v) const
	{
		return _B.Succ0(v) + 1 ;
	}

	size_t PrevSibling(size_t v) const
	{
		return _B.Pred0(v - 2) + 1 ; 
	}
  
  size_t Parent(size_t v) const
  {
    if (v == Root())
      return 0 ;
    size_t j = _B.Select(1, _B.Rank(0, v - 1)) ;
    return _B.Pred0(j) + 1 ;
  }

  bool IsLeaf(size_t v) const
  {
    return _B.Access(v) == 0 ; 
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
    // Exclude current 0 in case v is the leaf.
		return _B.Rank(0, v, /*inclusive=*/0) ; 
  }

	//Map actual node id to index in B (v).
  size_t NodeSelect(size_t i) const
  {
		if (i == 0)
			return Root() ;
		else
			return _B.Select(0, i) + 1 ; 
  }

  void Save(FILE *fp)
  {
		Tree::Save(fp) ;
		_B.Save(fp) ;
  }

  void Load(FILE *fp)
  {
		Tree::Load(fp) ;
		_B.Load(fp) ;
	}
} ;

} // end of name space

#endif

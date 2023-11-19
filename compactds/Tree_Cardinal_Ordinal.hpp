#ifndef _MOURISL_COMPACTDS_TREE_CARDINAL_ORDINAL
#define _MOURISL_COMPACTDS_TREE_CARDINAL_ORDINAL

// Cardinal tree, where we use ordinal compact tree to store the structure
//  and another bit vector to represent the concatenating labels

#include "Tree_Cardinal.hpp"
#include "Tree_Cardinal_Plain.hpp"
#include "Bitvector_Plain.hpp"

#include "Tree_BP.hpp"
#include "Tree_DFUDS.hpp"

namespace compactds {

template <class TreeClass = Tree_DFUDS, class BvClass = Bitvector_Plain>
class Tree_Cardinal_Ordinal: public Tree_Cardinal
{
private:
  TreeClass _t ;
  BvClass _B ; // concatenated labeling showing whether the children for this label exist
public:
  size_t GetSpace(bool inclusive = true) 
  {
    return _space + (inclusive ? sizeof(*this) : 0) ;
  }
  
  // Use c bits per node.
  void Init(const struct _plainCardinalTreeNode *treeNodes, size_t n, size_t c, size_t *treeIdMap)
  {
    size_t i, j ;

    _n = n ;
    _c = c ;
    _t.InitFromCardinalTree(treeNodes, n, c, treeIdMap) ;
     

    WORD *W = Utils::MallocByBits(c * n) ;
    for (i = 0 ; i < n ; ++i)
    {
      size_t k = treeIdMap[i] ;
      for (j = 0 ; j < c ; ++j)
      {
        if (treeNodes[i].children[j] != 0)
          Utils::BitSet(W, _c * k + j) ;
      }
    }
    _B.Init(W, c * n) ;
    free(W) ;

    _space = _t.GetSpace(false) + _B.GetSpace() - sizeof(_B) ;
  }

  // The index in B
  size_t Root() const
  {
    return _t.Root() ;
  }

  // @return: the t-th child (1-based) of node v in B vector
  size_t ChildSelect(size_t v, size_t t) const
  {
    return _t.ChildSelect(v, t) ;
  }

  size_t FirstChild(size_t v) const 
  {
    return _t.FirstChild(v) ;
  }

  size_t LastChild(size_t v) const 
  {
    return _t.LastChild(v) ;
  }
  
  size_t ChildrenCount(size_t v) const
  {
    return _t.ChildrenCount(v) ;
  }
  
  // return: v is the ret-th (1-based) child of the parent.
  size_t ChildRank(size_t v) const
  {
    return _t.ChildRank(v) ;
  }

  // The silbing function assumes v has
  //   those siblings.
  size_t NextSibling(size_t v) const
  {
    return _t.NextSibling(v) ;
  }

  size_t PrevSibling(size_t v) const
  {
    return _t.PrevSibling(v) ;
  }
  
  size_t Parent(size_t v) const
  {
    return _t.Parent(v) ;
  }

  bool IsLeaf(size_t v) const 
  {
    return _t.IsLeaf(v) ;
  }

  size_t LCA(size_t u, size_t v) const
  {
    return _t.LCA(u, v) ;
  }
  
  // Maps index in B (v) back up to the actual node id
  size_t NodeMap(size_t v) const
  {
    return _t.NodeMap(v) ;
  }

  //Map actual node id to index in B (v).
  size_t NodeSelect(size_t i) const 
  {
    return _t.NodeSelect(i) ;
  }

  // Number of children with label l. 1: has such children. 0-don't 
  size_t ChildrenLabeled(size_t v, size_t l) const 
  {
    return (_B.Access(NodeMap(v) * _c + l) == 1) ? 1 : 0 ;
  }

  // The child with label l.
  // Assuming label l's child exist
  // Notice the difference from Child()
  size_t LabeledChild(size_t v, size_t l) const 
  {
    size_t k = NodeMap(v) ;
    size_t r = _B.Rank1(k * _c + l) - _B.Rank1(k * _c) ;
    return ChildSelect(v, r + 1) ;
  }

  // The label of the edge that leads to node v.
  size_t ChildLabel(size_t v) const 
  {
    size_t p = NodeMap(Parent(v)) ;
    size_t r = ChildRank(v) ;
    return _B.Select( _B.Rank1(p * _c, 0) + r) - p * _c ; 
  }

  void Save(FILE *fp)
  {
    Tree_Cardinal::Save(fp) ;
    _t.Save(fp) ;
    _B.Save(fp) ;
  }

  void Load(FILE *fp)
  {
    Tree_Cardinal::Load(fp) ;
    _t.Save(fp) ;
    _B.Load(fp) ;
  }
} ;

} // end of name space

#endif

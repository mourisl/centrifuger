#ifndef _MOURISL_COMPACTDS_TREE_LABELED
#define _MOURISL_COMPACTDS_TREE_LABELED

// Labeled tree. The difference from the text book that we
//   assume general representation of the tree structure, including BP.
//   Therefore, for the concatenated children labels, we have an
//   additional bit vector to indicate the start of each children label series,
//   and we also have a place holder for the start position in the labels
//   (this could be redundant, but having a explicit bit marker is more efficient)

#include "Bitvector_Plain.hpp"
#include "Tree.hpp"

#include "Sequence_WaveletTree.hpp"

#include "Tree_LOUDS.hpp"
#include "Tree_BP.hpp"
#include "Tree_DFUDS.hpp"

#include <map>
#include <vector>

namespace compactds {

template <class TreeClass = Tree_DFUDS, class SequenceClass = Sequence_WaveletTree<>, 
         class SequenceMarkerClass = Bitvector_Plain>
class Tree_Labeled: public Tree
{
private:
  TreeClass _t ;
  SequenceClass _l ;
  SequenceMarkerClass _lmarker ; // markers on l, indicating the start of the labels from a node
  std::map<size_t, ALPHABET> _lmap ; // mapping label from size_to to ALPHABET
  std::vector<ALPHABET> _lmapback ; // map labels from ALPAHBET back to original value
public:
  size_t GetSpace(bool inclusive = true) 
  {
    return _space + _lmap.size() * (sizeof(size_t) + sizeof(ALPHABET)) + _lmapback.capacity() + (inclusive ? sizeof(*this) : 0) ;
  }
  
  // Use c bits per node.
  void Init(const struct _plainTreeNode *treeNodes, size_t n,
      size_t *treeIdMap)
  {
    size_t i ;
    _n = n ;
    _t.Init(treeNodes, n, treeIdMap) ;
    
    for (i = 1 ; i < n ; ++i)
    {
      if (_lmap.find(treeNodes[i].label) == _lmap.end())
      {
        ALPHABET id = _lmap.size() ;
        _lmap[ treeNodes[i].label ] = id ;
        _lmapback.push_back( treeNodes[i].label ) ;
      }
    }

    std::vector<ALPHABET> alphabetList ;
    size_t lmapSize = _lmap.size() ;
    for (i = 0 ; i < lmapSize ; ++i)
      alphabetList.push_back((ALPHABET)i) ;
    
    ALPHABET placeHolder = lmapSize ;
    alphabetList.push_back(placeHolder) ;
    ++lmapSize ;

    FixedSizeElemArray childrenLabels ;
    Alphabet labelAlphabet ;
    int lmapBits = labelAlphabet.InitFromList(alphabetList.data(), lmapSize) ;
    childrenLabels.Malloc(lmapBits, 2 * n - 1) ;
    
    WORD *W = Utils::MallocByBits(2 * n - 1) ;
    size_t lused = 0 ;
    size_t *nodeOrder = (size_t *)malloc(sizeof(*nodeOrder) * n) ;
    for (i = 0 ; i < n ; ++i)
      nodeOrder[ treeIdMap[i] ] = i ;
    for (i = 0 ; i < n ; ++i)
    {
      size_t k = nodeOrder[i] ;
      size_t c = treeNodes[k].child ;
      
      Utils::BitSet(W, lused) ; //mark the start of the child series
      childrenLabels.Write(lused, placeHolder) ;
      ++lused ;
      while (c != 0)
      {
        childrenLabels.Write(lused, _lmap[treeNodes[c].label]) ;
        ++lused ;

        c = treeNodes[c].sibling ;
      }
    }
    free(nodeOrder) ;

    _l.SetAlphabet(labelAlphabet) ;
    _l.Init(childrenLabels, 2 * n - 1, alphabetList.data()) ;
    _lmarker.Init(W, 2 * n - 1) ;

    free(W) ;
    
    _space = _t.GetSpace(false) + _l.GetSpace() - sizeof(_l) + _lmarker.GetSpace() - sizeof(_lmarker) ;
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
    size_t i = NodeMap(v) ;
    ALPHABET lmapped = _lmap.at(l) ;
    
    size_t start = _lmarker.Select(i + 1) ;
    if (i == _n - 1)
    {
      return _l.Rank(lmapped, 2 * _n - 2) - _l.Rank(lmapped, start);
    }
    else
    {
      return _l.Rank(lmapped, _lmarker.Select(i + 2) - 1) - _l.Rank(lmapped, start) ; 
    }
  }

  // The t-th child with label l.
  size_t LabeledChildSelect(size_t v, size_t l, size_t t) const 
  {
    size_t i = NodeMap(v) ;
    ALPHABET lmapped = _lmap.at(l) ;
    size_t start = _lmarker.Select(i + 1) ;

    size_t childRank = _l.Select(lmapped, _l.Rank(lmapped, start) + t) - start ;
    return ChildSelect(v, childRank) ;
  }

  // The label of the edge that leads to node v.
  size_t ChildLabel(size_t v) const 
  {
    if (v == Root())
      return 0 ;

    size_t childRank = ChildRank(v) ;
    size_t p = Parent(v) ;
    return _lmapback.at( _l.Access( _lmarker.Select(NodeMap(p) + 1) + childRank) ) ;
  }

  void Save(FILE *fp)
  {
    Tree::Save(fp) ;
    _t.Save(fp) ;
    _l.Save(fp) ;
    _lmarker.Save(fp) ;

    size_t size = _lmapback.size() ; 
    SAVE_VAR(fp, size) ;
    size_t i ;
    for (i = 0 ; i < size ; ++i) 
    {
      SAVE_VAR(fp, _lmapback[i]) ;
    }
  }

  void Load(FILE *fp)
  {
    Tree::Load(fp) ;
    _t.Save(fp) ;
    _l.Load(fp) ;
    _lmarker.Save(fp) ;

    _lmap.clear() ;
    _lmapback.clear() ;
    size_t lmapSize ;
    LOAD_VAR(fp, lmapSize) ;
    size_t i ;
    for (i = 0 ; i < lmapSize ; ++i)
    {
      size_t l ;
      LOAD_VAR(fp, l) ;
      _lmapback.push_back(l) ;
      _lmap[l] = i ;
    }
  }
} ;

} // end of name space

#endif

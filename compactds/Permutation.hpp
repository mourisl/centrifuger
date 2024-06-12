#ifndef _MOURISL_COMPACTDS_PERMUTATION
#define _MOURISL_COMPACTDS_PERMUTATION

// Compressed permutation representation. Chapter 5.3
// So far it assuems at most 2^31 runs
#include "Utils.hpp"
#include "HuffmanCode.hpp"
#include "Bitvector_Plain.hpp"

namespace compactds {
class Permutation
{
private:
  size_t _space ;
  size_t _n ;
  size_t _rcnt ;
  Bitvector_Plain *_nodeB ; // The left, right child indicator 
  Bitvector_Plain _G ; // Mark the start position for each run in the permuation representation.
  int *_nodePath ; //Buffer to hold the node ids along a path from root to leaf
  HuffmanCode _huffmanTree ;

  // Combine the CreateLeaves and CreateBitvectors of the book into the same function
  // Using S to hold the sequences mimicing merge sort
  // tag: tree id
  void CreateBitvectors(const struct _huffman_node *tree, int tag, size_t *S, size_t offset, size_t *Pi)
  {
    size_t i ;
    if (tree[tag].left == -1)
    {
      // Leaf
      for (i = 0 ; i < tree[tag].freq ; ++i)
        S[offset + i] = Pi[ _G.Select(1, tree[tag].symbol + 1) + i] ;
    }
    else
    {
      CreateBitvectors(tree, tree[tag].left, S, offset, Pi) ;
      CreateBitvectors(tree, tree[tag].right, S, offset + tree[tree[tag].left].freq, Pi) ;

      _nodeB[tag].Malloc(tree[tag].freq) ;
      _space += (_nodeB[tag].GetSpace() - sizeof(_nodeB[tag])) ;

      // Merge
      size_t *buffer = (size_t *)malloc(sizeof(size_t) * tree[tag].freq) ;
      size_t lp = offset, rp = offset + tree[tree[tag].left].freq ; // left/right pointer
      size_t lcnt, rcnt ; // scanned count of left and right children
      lcnt = rcnt = 0 ;
      while (lcnt < tree[tree[tag].left].freq && rcnt < tree[ tree[tag].right ].freq)
      {
        if (S[lp] < S[rp])
        {
          buffer[lcnt + rcnt] = S[lp] ;
          ++lcnt ; ++lp ;
        }
        else if (S[lp] > S[rp])
        {
          buffer[lcnt + rcnt] = S[rp] ;
          _nodeB[tag].BitSet(lcnt + rcnt) ;
          ++rcnt ; ++rp ;
        }
        else
        {
          // ERROR!
        }
      }
      while (lcnt < tree[tree[tag].left].freq)
      {
        buffer[lcnt + rcnt] = S[lp] ;
        ++lcnt ; ++lp ;
      }
      while (rcnt < tree[tree[tag].right].freq)
      {
        buffer[lcnt + rcnt] = S[rp] ;
        _nodeB[tag].BitSet(lcnt + rcnt) ;
        ++rcnt ; ++rp ;
      }
      _nodeB[tag].Init() ;
      for (i = offset ; i < offset + tree[tag].freq ; ++i)
        S[i] = buffer[i - offset] ;
      free(buffer) ;
    }
  }
public:
  Permutation() 
  {
    _space = 0 ;
    _n = 0 ;
    _rcnt = 0 ;
  }

  ~Permutation()
  {
    Free() ;
  }

  void Free()
  {
    if (_n > 0)
    {
      delete[] _nodeB ;
      free(_nodePath) ;
      _n = 0 ;
    }
  }

  size_t GetSpace()
  {
    return _space + sizeof(*this) ; 
  }

  void Init(size_t *Pi, size_t n)
  {
    _n = n ;
    size_t i, j ;
    std::vector<size_t> rstarts ;
    std::vector<uint64_t> rlens ;
    _G.Malloc(n) ;
    for (i = 0 ; i < n ;)
    {
      for (j = i + 1 ; j < n ; ++j)
        if (Pi[j] < Pi[j - 1])
          break ;
      rstarts.push_back(i) ;
      rlens.push_back(j - i) ;
      _G.BitSet(i) ;
      i = j ;
    }
    _rcnt = rstarts.size() ;
    _G.Init() ;
    _space += _G.GetSpace() - sizeof(_G) ;

    _huffmanTree.InitFromFrequency(rlens.data(), _rcnt) ;
    _space += _huffmanTree.GetSpace() - sizeof(_huffmanTree) ;
    
    int depth = _huffmanTree.GetDepth( _huffmanTree.GetRoot() ) ;
    _nodePath = (int *)malloc(sizeof(_nodePath[0]) * (depth + 1)) ;
    _space += sizeof(_nodePath[0]) * (depth + 1) ;
  
    _nodeB = new Bitvector_Plain[2 * _rcnt - 1] ; 
    _space += sizeof(Bitvector_Plain) * (2 * _rcnt - 1) ; 
    size_t *S = (size_t *)malloc(sizeof(size_t) * n) ;
    CreateBitvectors(_huffmanTree.GetTree(), _huffmanTree.GetRoot(), S, 0, Pi) ;
    /*for (i = 0 ; i < _rcnt ; ++i)
    {
      int l ;
      WORD code = _huffmanTree.Encode(i, l) ;
      printf("%d %d %d: %d %d\n", i, rstarts[i], rlens[i], code, l) ;
    }*/
    free(S) ;
  }

  // Pi(i)
  // read() in the book
  size_t Next(size_t i) const
  {
    int j ;
    const struct _huffman_node *tree = _huffmanTree.GetTree() ;
    int len ;
    size_t ri = _G.Rank(1, i) - 1 ; 
    WORD code = _huffmanTree.Encode(ri, len) ;
    
    _nodePath[0] = _huffmanTree.GetRoot() ;
    for (j = 0 ; j < len ; ++j)
    {
      if ((code >> (len - j - 1)) & 1)
        _nodePath[j + 1] = tree[ _nodePath[j] ].right ;
      else
        _nodePath[j + 1] = tree[ _nodePath[j] ].left ;
    }
    
    i = i - _G.Select(1, ri + 1) ;
    for (j = len - 1 ; j >= 0 ; --j)
    {
      if ((code >> (len - 1 - j)) & 1)
      {
        i = _nodeB[_nodePath[j]].Select(1, i + 1) ;
      }
      else
      {
        i = _nodeB[_nodePath[j]].Select(0, i + 1) ;
      }
    }
    return i ;
  }

  // Pi^-1(i)
  // inverse() in the book
  size_t Prev(size_t i) const
  {
    size_t j = i ; // tracking the position of i in a run
    const struct _huffman_node *tree = _huffmanTree.GetTree() ; 
    size_t tag = _huffmanTree.GetRoot() ;
    while (tree[tag].left != -1)
    {
      int b = _nodeB[tag].Access(j) ;
      j = _nodeB[tag].Rank(b, j, /*inclusive=*/0) ;
      if (b == 0)
        tag = tree[tag].left ;
      else
        tag = tree[tag].right ;
    }

    return _G.Select(1, tree[tag].symbol + 1) + j ;
  }

  void Save(FILE *fp)
  {
    SAVE_VAR(fp, _n) ;
    SAVE_VAR(fp, _space) ;
    SAVE_VAR(fp, _rcnt) ;
    _G.Save(fp) ;
    _huffmanTree.Save(fp) ;
    size_t i ;
    for (i = 0 ; i < 2 * _rcnt - 1 ; ++i)
      _nodeB[i].Save(fp) ;
  }

  void Load(FILE *fp)
  {
    Free() ;

    LOAD_VAR(fp, _n) ;
    LOAD_VAR(fp, _space) ;
    LOAD_VAR(fp, _rcnt) ;
    _G.Load(fp) ;
    _huffmanTree.Load(fp) ;
    size_t i ;
    _nodeB = new Bitvector_Plain[2 * _rcnt - 1] ; 
    for (i = 0 ; i < 2 * _rcnt - 1 ; ++i)
      _nodeB[i].Load(fp) ;
    
    int depth = _huffmanTree.GetDepth( _huffmanTree.GetRoot() ) ;
    _nodePath = (int *)malloc(sizeof(_nodePath[0]) * (depth + 1)) ;
  }
} ;
}

#endif

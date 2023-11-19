#ifndef _MOURISL_COMPACTDS_HUFFMANCODE
#define _MOURISL_COMPACTDS_HUFFMANCODE

#include <algorithm> 

#include "Utils.hpp"

namespace compactds {
struct _huffman_node
{
  int symbol ;
  uint64_t freq ;
  int next ; // used in _tree construction as a linked list
  int left, right ; // Left, right children
  bool operator <(const struct _huffman_node &b) const
  {
    return freq < b.freq ; 
  }
} ;

class HuffmanCode
{
private:
  WORD *_codes ; // assume alphabet set is in [0, n-1].
  int *_codeLens ; 
  size_t _n ; // the size of the alphabet 
  struct _huffman_node *_tree ;
  size_t _space ;

  // Algorithm 2.2: building a huffman _tree with linked list instead of heap
  void BuildTree(struct _huffman_node *elems, size_t n)
  {
    std::sort(elems, elems + n) ;
    
    size_t i ;
    size_t nodeCnt ;
    
    _tree = (struct _huffman_node *)malloc(sizeof(*_tree) * (2 * n - 1)) ;
    _space += (sizeof(*_tree) * (2 * n - 1)) ;
    
    for (i = 0 ; i < n ; ++i)
    {
      _tree[i] = elems[i] ;
      if (i + 1 < n)
        _tree[i].next = i + 1 ;
      else
        _tree[i].next = -1 ;
      _tree[i].left = _tree[i].right = -1 ;
    }
    size_t minTag = 0 ; // minTag and minTag+1 is the availble two nodes with minimum
    size_t insertTag = 0 ; // the start position to search the next insert node
                    // this marker is the key for linear time building the _tree after sorting.
    nodeCnt = n ;
    int p ;
    while (1)
    {
      int a = minTag ;
      int b = _tree[minTag].next ;
      if (b == -1)
        break ;
      _tree[nodeCnt].symbol = -1 ;
      _tree[nodeCnt].freq = _tree[a].freq + _tree[b].freq ;
      _tree[nodeCnt].left = a ;
      _tree[nodeCnt].right = b ;

      // Search for the appropriate position to insert the new element
      p = insertTag ;
      while (_tree[p].next != -1 && _tree[ _tree[p].next ].freq <= _tree[nodeCnt].freq)
        p = _tree[p].next ;
      
      _tree[nodeCnt].next = _tree[p].next ;
      _tree[p].next = nodeCnt ;
      
      insertTag = nodeCnt ; 
      ++nodeCnt ;
      minTag = _tree[b].next ;
    }
  }

  // Recurisvely traverse the Huffman _tree to put the code
  void CreateCodes(int tag, WORD c, int l)
  {
    if (_tree[tag].left == -1 && _tree[tag].right == -1)
    {
      _codes[_tree[tag].symbol] = c ;
      _codeLens[_tree[tag].symbol] = l ;
      return ;
    }

    CreateCodes(_tree[tag].left, c<<1, l + 1) ;
    CreateCodes(_tree[tag].right, (c<<1) + 1, l + 1) ;
  }

  void InternalInit(struct _huffman_node *elems, size_t n)
  {
    this->_n = n ;
    _space = 0 ;

    BuildTree(elems, n) ;
    _codes = (WORD *)malloc(sizeof(*_codes) * n) ;
    _codeLens = (int *)malloc(sizeof(*_codeLens) * n) ;
    CreateCodes(2*n - 2, 0, 0) ; 
  }

public:
  HuffmanCode() 
  {
    _n = _space = 0 ;
    _codes = NULL ;
    _codeLens = NULL ;
    _tree = NULL ;
  }
  ~HuffmanCode() {Free();}

  void Free()
  {
    _n = _space = 0 ;
    if (_codes != NULL)
    {
      free(_codes) ; free(_codeLens) ; free(_tree) ;
      _codes = NULL ;
      _codeLens = NULL ;
      _tree = NULL ;
    }
  }

  size_t GetSpace()
  {
    return _space + sizeof(*this) ;
  }
  
  int GetSize()
  {
    return _n ;
  }

  struct _huffman_node *GetTree() const
  {
    return _tree ;
  }

  size_t GetRoot() const
  {
    return 2 * _n - 2 ; 
  }

  HuffmanCode &operator =(const HuffmanCode &in)
  {
    Free() ;
    
    if (in._n == 0)
      return *this;
    _n = in._n ;
    _space = in._space ;

    _codes = (WORD *)malloc(sizeof(*_codes) * _n) ;
    _codeLens = (int *)malloc(sizeof(*_codeLens) * _n) ;
    _tree = (struct _huffman_node *)malloc(sizeof(*_tree) * (2*_n-1)) ;

    memcpy(_codes, in._codes, sizeof(*_codes) * _n) ;
    memcpy(_codeLens, in._codeLens, sizeof(*_codeLens) * _n) ;
    memcpy(_tree, in._tree, sizeof(*_tree)) ;

    return *this ;
  }

  void InitFromFrequency(const uint64_t *freq, const size_t n)
  {
    size_t i ;
    struct _huffman_node *elems = (struct _huffman_node*)malloc(sizeof(*elems) * n);
    
    for (i = 0 ; i < n ; ++i)
    {
      elems[i].symbol = i ;
      elems[i].freq = freq[i] ;
    }
    InternalInit(elems, n) ;

    free(elems) ;
  }

  int GetDepth(int tag)
  {
    if (_tree[tag].left == -1)
      return 0 ;
    int ldepth = GetDepth(_tree[tag].left) ; 
    int rdepth = GetDepth(_tree[tag].right) ;
    return 1 + (ldepth > rdepth ? ldepth : rdepth) ;
  }

  WORD Encode(int x, int &l) const 
  {
    l = _codeLens[x] ;
    return _codes[x] ;
  }

  int Decode(WORD c, int l) const
  {
    int i ;
    int p = 2 * _n - 2 ; // root
    for (i = 0 ; i < l ; ++i)
    {
      if ((c >> (l - i - 1)) & 1)
        p = _tree[p].right ;
      else
        p = _tree[p].left ;
    }
    return _tree[p].symbol ;
  }

  void Save(FILE *fp)
  {
    fwrite(this, sizeof(this), 1, fp) ;
    fwrite(_tree, sizeof(_tree[0]), 2 * _n - 1, fp) ;
  }

  void Load(FILE *fp)
  {
    Free() ;

    fread(this, sizeof(this), 1, fp) ;
    _codes = (WORD *)malloc(sizeof(*_codes) * _n) ;
    _codeLens = (int *)malloc(sizeof(*_codeLens) * _n) ;
    _tree = (struct _huffman_node *)malloc(sizeof(*_tree) * (2*_n-1)) ;
    fwrite(_tree, sizeof(_tree[0]), 2 * _n - 1, fp) ;
    CreateCodes(2*_n - 2, 0, 0) ; 
  }
} ;
}
#endif

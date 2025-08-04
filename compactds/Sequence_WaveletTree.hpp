#ifndef _MOURISL_COMPACTDS_SEQUENCE_WAVELETTREE
#define _MOURISL_COMPACTDS_SEQUENCE_WAVELETTREE

#include "Utils.hpp"
#include "Sequence.hpp"

#include <string.h>

#include "Bitvector_Plain.hpp"
#include "Bitvector_RunLength.hpp"

namespace compactds {
template <class BvClass>
struct _sequence_wavelettree_node
{
  BvClass v ;
  WORD prefix ; // the 
  int prefixLen ; // bits in prefix
  int children[2] ;

  void Save(FILE *fp) 
  {
    SAVE_VAR(fp, prefix) ;
    SAVE_VAR(fp, prefixLen) ;
    SAVE_ARR(fp, children, 2) ;
    v.Save(fp) ;
  }

  void Load(FILE *fp)
  {
    LOAD_VAR(fp, prefix) ;
    LOAD_VAR(fp, prefixLen) ;
    LOAD_ARR(fp, children, 2) ;
    v.Load(fp) ;
  }
} ;

// The implementation of wavelet tree in either
// perfect balanced or huffman shape,
// depending on the choice of alphabet.
template <class BvClass = Bitvector_Plain>
class Sequence_WaveletTree: public Sequence
{
private:
  struct _sequence_wavelettree_node<BvClass> *_T ; 
  int _tNodeCnt ;
  int _selectSpeed ;

  // Based on the pos-th bits (0-index, count from leftside)
  // maxPosToRight: record the maximum distance from pos to right side. 
  // return: the number of 1s
  uint64_t ConvertSequenceToBits(const FixedSizeElemArray &S, const ALPHABET *alphabetMap, int pos, WORD *v, int &maxPosToRight)
  {
    size_t i ;
    size_t n = S.GetSize() ;
    uint64_t ret = 0;
    maxPosToRight = 0 ;

    for (i = 0 ; i < n ; ++i)
    { 
      int codeLen = 0 ;
      int b = _alphabets.Encode( alphabetMap[S.Read(i)], codeLen ) ;
      if (codeLen - pos > maxPosToRight)
        maxPosToRight = codeLen - pos ;
      if (b & (1<<(codeLen - pos - 1)))
      {
        Utils::BitSet(v, i) ; // Assume the array is already initated to be all 0
        ++ret ;
      }
    }
    return ret ;
  }
  
  // Assume left and right's memory has been allocated.
  void SplitSequence(const FixedSizeElemArray &orig, WORD *v,
      FixedSizeElemArray &left, FixedSizeElemArray &right)
  {
    size_t i ;
    size_t len = orig.GetSize() ;
   
    size_t leftLen = 0 ;
    size_t rightLen = 0 ;
    for (i = 0 ; i < len ; ++i)
    {
      if (!Utils::BitRead(v, i))
      {
        left.Write(leftLen, orig.Read(i)) ; 
        ++leftLen ;
      }
      else
      {
        right.Write(rightLen, orig.Read(i)) ;
        ++rightLen ;
      }
    }
  }

  // The recursive function that construct the tree.
  // Assumes that the leaf node always has the brother.
  // depth: how many bits has been processed so far
  // tused: the number of wavelet tree node used so far
  // bufferv: preallcoated memory to holding temporary bit array
  // return: node id (index in T)
  int BuildTree(const FixedSizeElemArray &S, const ALPHABET *alphabetMap, int depth, WORD prefix, WORD *bufferv)
  {
    size_t len = S.GetSize() ;
    int ti = _tNodeCnt ;
    ++_tNodeCnt ;
    int remainingBits ;
    
    memset(bufferv, 0, Utils::BitsToWordBytes(len)) ;
    uint64_t onecnt = ConvertSequenceToBits(S, alphabetMap, depth, bufferv, remainingBits) ;
    _T[ti].v.SetSelectSpeed(_selectSpeed) ;
    _T[ti].v.Init(bufferv, len) ;
    _space += _T[ti].v.GetSpace() - sizeof(_T[ti].v) ;
    _T[ti].prefix = prefix ;
    _T[ti].prefixLen = depth ;
    if (remainingBits == 1 || S.GetSize() == 0)
    {
      // Reach leaf.
      _T[ti].children[0] = _T[ti].children[1] = -1 ;
      return ti ;
    }
    FixedSizeElemArray leftS, rightS ; // the memory should be automatically released
    leftS.Malloc(S.GetElemLength(), len - onecnt) ;
    rightS.Malloc(S.GetElemLength(), onecnt) ;
    SplitSequence(S, bufferv, leftS, rightS) ;

    _T[ti].children[0] = BuildTree(leftS, alphabetMap, depth + 1, prefix << 1, bufferv) ;
    _T[ti].children[1] = BuildTree(rightS, alphabetMap, depth + 1, (prefix << 1) | 1ull, bufferv) ;

    return ti ;
  }
  
  int AccessInNode(int ti, size_t i) const
  {
    return _T[ti].v.Access(i) ;
  }

  // Calculate the rank(type, i) in T[ti]
  size_t RankInNode(int ti, int type, size_t i, int inclusive = 1) const
  {
    return _T[ti].v.Rank(type, i, inclusive) ;
  }

  size_t SelectInNode(int ti, int type, size_t i) const
  {
    return _T[ti].v.Select(type, i);
  }

  // Recursive function for Select
  // c: the code for the alphabet.
  // l: the length of the code
  // i: the select we want to query. The ith chracter c.
  // ti: tree node idx
  // depth: the recursive depth
  size_t RecursiveSelect(WORD c, int l, size_t i, int ti, int depth ) const
  {
    int b = (c >> (l-depth-1)) & 1 ;
    if (depth >= l - 1)
    {
      return SelectInNode(ti, b, i) ;
    }
    
    // Need the +1 to convert the index from Select to the rank as the input to Select
    return SelectInNode(ti, b, RecursiveSelect(c, l, i, _T[ti].children[b], depth + 1) + 1 ) ;
  }
public:
  Sequence_WaveletTree() 
  {
    _tNodeCnt = 0 ;
    _selectSpeed = BITVECTOR_DEFAULT_SELECT_SPEED ;
  }

  ~Sequence_WaveletTree() {Free() ;}
  
  void Free() 
  {
    if (_tNodeCnt)
    {
      delete[] _T ;
      _T = NULL ;
      _tNodeCnt = 0 ;
      Sequence::Free() ;
    }
  }

  void SetSelectSpeed(int speed)
  {
    _selectSpeed = speed ;
  }

  size_t GetSpace() {return _space + _alphabets.GetSpace() - sizeof(_alphabets) + sizeof(this) ;} 
  
  // We compactly represent the input sequence as fixed-size element array in a plain fashion
  // just to save some memory when construct the tree.
  void Init(const FixedSizeElemArray &S, size_t sequenceLength, const ALPHABET *alphabetMap)
  {
    _space = 0 ;
    this->_n = sequenceLength ;
    
    if (_alphabets.GetSize() == 0)
      _alphabets.InitFromList(alphabetMap, strlen(alphabetMap)) ;
     
    _T = new struct _sequence_wavelettree_node<BvClass>[_alphabets.GetAlphabetCapacity() - 1] ; 
    _tNodeCnt = 0 ;
    _space += sizeof(*_T) * (_alphabets.GetAlphabetCapacity() - 1) ;
    
    WORD *bufferv = Utils::MallocByBits(sequenceLength) ; 
    BuildTree(S, alphabetMap, 0, 0, bufferv) ;
    free(bufferv) ;
  }

  // Return: the alphabet at position i.
  ALPHABET Access(size_t i) const 
  {
    int l = 0 ;
    WORD code = 0 ;
    int ti = 0 ;
    for (l = 0 ; ti != -1 ; ++l)
    {
      int b = AccessInNode(ti, i) ;
      code = (code << 1) | b ;
      // Need -1 to convert the rank number to array index.
      // There is no need to check the negativity from -1,
      // because we know the current bit is 0, so rank>=1.
      //i = _T[ti].v.Rank(b, i) - 1 ;
      i = RankInNode(ti, b, i) - 1 ;
      ti = _T[ti].children[b] ;
    }
    return _alphabets.Decode(code, l) ;
  }

  // Return: the number of alphabet c's in [0..i]  
  size_t Rank(ALPHABET c, size_t i, int inclusive = 1) const 
  {
    int l = 0 ; // l: the length of the code
    WORD code = _alphabets.Encode(c, l) ;
    int depth = 0 ;
    int ti = 0 ;
    if (!inclusive) // Since in the wavelet tree, the non-inclusive operation should only
        // happen in the leaf node, so we directly modify i here for simplicity.
    {
      if (i == 0)
        return 0 ;
      else
        --i ;
    }
    for (depth = 0 ; depth < l ; ++depth)
    {
      int b = (code >> (l - depth - 1)) & 1 ;
      
      //i = _T[ti].v.Rank(b, i) ;
      i = RankInNode(ti, b, i) ;
      
      if (i == 0 || depth == l - 1)
        break ;
      // R count the number of 1's (or 0's), so we need to -1 to change it to 0-based index
      //  as in the bitvector.
      --i ; 
      ti = _T[ti].children[b] ;
    }
    return i ;
  }

  // Return: rank of c in [0..i] (inclusive), 
  //  also test whether T[i]==c, return through isC
  size_t RankAndTest(ALPHABET c, size_t i, bool &isC) const
  { 
    int l = 0 ; // l: the length of the code
    WORD code = _alphabets.Encode(c, l) ;
    int depth = 0 ;
    int ti = 0 ;

    isC = true ;
    for (depth = 0 ; depth < l ; ++depth)
    {
      int b = (code >> (l - depth - 1)) & 1 ;
      if (isC && b != AccessInNode(ti, i))
        isC = false ;

      //i = _T[ti].v.Rank(b, i) ;
      i = RankInNode(ti, b, i) ;

      if (i == 0 || depth == l - 1)
        break ;
      // R count the number of 1's (or 0's), so we need to -1 to change it to 0-based index
      //  as in the bitvector.
      --i ; 
      ti = _T[ti].children[b] ;
    }
    return i ;
  }

  // return: the index of the ith (1-based) c 
  size_t Select(ALPHABET c, size_t i) const 
  {
    int l = 0 ;
    WORD code = _alphabets.Encode(c, l) ;
    return RecursiveSelect(code, l, i, 0, 0) ;
  }

  void Save(FILE *fp)
  {
    Sequence::Save(fp) ;
    SAVE_VAR(fp, _tNodeCnt) ;
    SAVE_VAR(fp, _selectSpeed) ;
    int i ;
    for (i = 0 ; i < _tNodeCnt ; ++i)
      _T[i].Save(fp) ;
  }

  void Load(FILE *fp)
  {
    Free() ;
    Sequence::Load(fp) ;
    LOAD_VAR(fp, _tNodeCnt) ;
    LOAD_VAR(fp, _selectSpeed) ;

    if (_alphabets.GetSize() == 0) //empty tree
      return ;

    _T = new struct _sequence_wavelettree_node<BvClass>[_tNodeCnt] ; 
    int i ;
    for (i = 0 ; i < _tNodeCnt ; ++i)
      _T[i].Load(fp) ;
  }

  void PrintStats()
  {
    Utils::PrintLog("Sequence_WaveletTree: total_length: %lu node_cnt: %lu", _n, _tNodeCnt) ;
  }

} ;
}

#endif

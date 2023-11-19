#ifndef _MOURISL_COMPACTDS_BITVECTOR_COMPRESSED
#define _MOURISL_COMPACTDS_BITVECTOR_COMPRESSED

#include "Utils.hpp"
#include "Bitvector.hpp"

#include "FixedSizeElemArray.hpp"

// The compressed bitvector based on chaptor 4
// This seems to be the RRR bitvector
namespace compactds {
class Bitvector_Compressed : public Bitvector
{
private:
  int _b ; // block size for bit vector
  int _pb ; // block size for partial sum array _P
  size_t _n ; // the total raw length of the bits  
 
  // Variables for compress the bit vector
  FixedSizeElemArray _C ; // the array for the count of bits
  WORD *_O ; // encoded offsets within each block
  size_t *_P ; // partial sum on offset array _O

  uint64_t **_choose ; // _C(i, j)
  int *_L ; // the required bits for _O for each _Ci

  // Variables for ranking query
  uint64_t *_R ; // precomputed rank (right-exclusive). R and _P are aligned

  // Variables for selection
  size_t *_S ; 
  int _sb ; // the block size for selection
  size_t _sBlockCnt ; 

  int _selectSpeed ;

  void EncodeBits(const WORD &B, int &c, size_t &o) const
  {
    int i ;
    WORD maskedB = B & MASK(_b) ;
    int onecnt = Utils::Popcount(maskedB) ;
    o = 0 ;
    c = 0 ;
    for (i = _b - 1 ; i >= 0 ; --i)
    {
      if ((maskedB >> i) & 1)
      {
        o += _choose[i][onecnt - c] ;
        ++c ;
      }
    }
  }

  WORD DecodeBits(int c, size_t o) const
  {
    WORD ret = 0 ;
    int usedOnes = 0 ;
    int i ;
    for (i = _b - 1 ; i >= 0 ; --i)
    {
      ret <<= 1 ;
      if (o >= _choose[i][c - usedOnes]) 
      {
        ret |= 1 ;
        o -= _choose[i][c - usedOnes];
        ++usedOnes ;
      }
    }

    return ret ;
  }

  void InitChoose(int b)
  {
    int i, j ;

    // Build the _choose array
    _choose = (uint64_t**)malloc(sizeof(*_choose) * (b+1) ) ;
    _space += sizeof(*_choose) * (b + 1) ;
    for (i = 0 ; i <= b ; ++i)
    {
      _choose[i] = (uint64_t*)malloc(sizeof(**_choose) * (i + 2)) ;
      _space += sizeof(**_choose) * (i + 2) ;
    }
    for (i = 0 ; i <= b ; ++i)
    {
      _choose[i][0] = 1 ;
      for (j = 1 ; j < i ; ++j)
      {
        _choose[i][j] = _choose[i - 1][j - 1] + _choose[i - 1][j] ;
      }
      _choose[i][i] = 1 ;
      _choose[i][i + 1] = 0 ;
    }

    _L = (int *)malloc(sizeof(*_L) * (b + 1)) ;
    _space += sizeof(*_L) * (b + 1) ;
    for (i = 0 ; i <= b ; ++i)
    {
      // There are _choose[b][i] different combinations for each _C_i
      _L[i] = Utils::Log2Ceil(_choose[b][i]) ;
    }
  }

public:
  Bitvector_Compressed() 
  {
    _n = _b = _pb = _sb = 0 ;
    _selectSpeed = BITVECTOR_DEFAULT_SELECT_SPEED ;
  }
  ~Bitvector_Compressed() {Free();}
  
  // blockSize should be 2^x - 1.so _C can be fully utilized
  void SetBlockSize(int blockSize)
  {
    _b = blockSize ;
  }

  void SetPsumBlockSize(int psumBlockSize)
  {
    _pb = psumBlockSize ;
  }

  void SetSelectBlockSize(int selectBlockSize)
  {
    _sb = selectBlockSize ;
  }

  void SetSelectSpeed(int in)
  {
    _selectSpeed = in ;
  }

  // W is the plain bit vector
  void Init(const WORD *W, const size_t n)
  {
    size_t i, j ;
    this->_n = n ;
    _space = 0 ;

    if (_b <= 0)
      _b = 8 * sizeof(WORD) - 1 ;

    if (_pb <= 0)
      _pb = 8 * sizeof(size_t) ;
   
    if (_sb <= _b)
      _sb = 8 * sizeof(size_t) * 8 * sizeof(size_t) ;
    size_t blockCnt = DIV_CEIL(n, _b) ;
    InitChoose(_b) ; // Initialize _choose and _L
   
    _C.Malloc(Utils::Log2Ceil(_b + 1), blockCnt) ; 
    _space += _C.GetSpace() ;

    // _Calculate the size for _O 
    size_t offsetsSize = 0 ;
    int maxOneCntInBlock = 0 ;
    uint64_t totalOneCnt = 0 ;
    for (i = 0 ; i < _n ; i += _b)
    {
      int onecnt = Utils::Popcount( Utils::BitsRead(W, i, (i + _b < _n ? i + _b - 1 : _n - 1)) ) ;

      offsetsSize += _L[onecnt] ;
      totalOneCnt += onecnt ;
      if (onecnt > maxOneCntInBlock)
        maxOneCntInBlock = onecnt ;
    }
    
    _O = Utils::MallocByBits(offsetsSize) ;
    _space += Utils::BitsToWordBytes(offsetsSize) ; 
  
    size_t psumBlockCnt = DIV_CEIL(blockCnt, _pb) ;
    _P = (size_t *)malloc(sizeof(size_t) * psumBlockCnt) ;
    _space += sizeof(size_t) * psumBlockCnt ;
    
    _R = (uint64_t *)malloc(sizeof(uint64_t) * psumBlockCnt) ;
    _space += sizeof(uint64_t) * psumBlockCnt ;

    if (_selectSpeed > 0)
    {
      _sBlockCnt = DIV_CEIL(totalOneCnt, _sb) ;
      _S = (size_t *)malloc(sizeof(size_t) * _sBlockCnt) ;
      _space += sizeof(uint64_t) * _sBlockCnt ;
    }

    // Build the _C, _O, _P that compress the bit vector
    // Also build the data structures for rank and selections
    // j is used to index _O
    size_t blocki ;
    uint64_t onecntSum = 0 ;
    bool locateFirstOne = false ;
    for (i = 0, j = 0, blocki = 0 ; i < _n ; i += _b, ++blocki)
    {
      WORD bits = Utils::BitsRead(W, i, (i + _b < _n ? i + _b - 1 : _n - 1)) ;
      
      int tmpc ;
      size_t tmpo ;
      EncodeBits(bits, tmpc, tmpo) ;
       
      //printf("%d %llu. %llu\n", tmpc, tmpo, bits) ;  
      if (blocki % _pb == 0)
      {
        _P[blocki/_pb] = j ;
      }
      _C.Write(blocki, tmpc) ;
      if (_L[tmpc] > 0)
      {
        Utils::BitsWrite(_O, j, j + _L[tmpc] - 1, tmpo) ;
      }
      j += _L[tmpc] ;
      
      // _Process the information for rank operation
      if (blocki % _pb == 0)
        _R[blocki / _pb]  = onecntSum ;
      
      // _Process the information for select operation 
      if (_selectSpeed && (onecntSum / _sb != (onecntSum + tmpc) / _sb 
          || (!locateFirstOne && tmpc > 0)))
      {
        int localOneCnt = 0 ; 
        int l = 0 ;
        for (l = 0 ; l < _b ; ++l)
          if ((bits >> l)&1)
          {
            if ((onecntSum + localOneCnt) % _sb == 0)
            {
              _S[(onecntSum + localOneCnt) / _sb] = i + l ;
              break ;
            }
            ++localOneCnt ;
          }
        locateFirstOne = true ;
      }  
      onecntSum += tmpc ;
    }
  }

  void Free()
  {
    _C.Free() ;
    if (_n != 0)
    {
      int i ;
      for (i = 0 ; i <= _b ; ++i)
        free(_choose[i]) ;
      free(_choose) ;
      free(_L) ;
      
      free(_O) ;
      free(_P) ;
        
      free(_R) ;

      free(_S) ;
      _n = 0 ;
    }
  }

  // Return the ith bits (0-based)
  int Access(size_t i) const
  {
    // Get the partial sum from _P
    size_t bi = i / _b ;
    size_t pi = bi / _pb ;

    size_t j ; 
    int blockc = _C.Read(bi) ;
    if (blockc == 0)
      return 0 ;
    else if (blockc == _b)
      return 1 ;
    
    size_t blocko ;
    size_t os = _P[pi] ; // start position in o
    // j to index the block offsets
    for (j = pi * _pb ; j < bi ; ++j)
      os += _L[ _C.Read(j) ] ;
    blocko = Utils::BitsRead(_O, os, os + _L[blockc] - 1) ;

    WORD bits = DecodeBits(blockc, blocko) ;

    int residuali = i % _b ;
    return (bits >> residuali) & 1 ;
  }

  // Return the _number of 1s before i
  size_t Rank1(size_t i, int inclusive = 1) const
  {
    size_t j ;
    size_t bi = i / _b ; // index for block
    size_t ri = bi / _pb ; // index for R

    size_t ret = _R[ri] ;
    size_t os = _P[ri] ; 
    for (j = ri * _pb ; j < bi ; ++j)
    {
      int onecnt = _C.Read(j) ;
      ret += onecnt ;
      os += _L[onecnt] ;
    }
    int blockc = _C.Read(bi) ;
    if (blockc == 0)
      return ret ;
    else if (blockc == _b)
      return ret + i%_b + inclusive ;
    else
    {
      size_t blocko = Utils::BitsRead(_O, os, os + _L[blockc] - 1) ;

      WORD bits = DecodeBits(blockc, blocko) ;

      int residuali = i % _b ;
      return ret + Utils::Popcount( bits & MASK_WCHECK(residuali + inclusive) ) ;
    }
  }

  // Return the index of th i-th (1-based, so rank and select are inversible) 1
  size_t Select(size_t i) const
  {
    if (i == 0)
      return POSITIVE_INF ;
    // Unlike the uncompressed case, binary search might be less efficient
    // because we _need to sequentially find the appropriate _O in rank.
    size_t j ;
    size_t si = (i-1) / _sb ; 
    size_t bi = _S[si] / _b ; // it aligns to block bi
    size_t pi = bi / _pb ; // block bi belongs to the _P-block recording the offset in _O

    size_t os = _P[pi] ;
    // We rollback the index a little bit to align with the information of _O
    uint64_t onecntSum = _R[pi] ; // Another bless that R and _P are aligned
    // j index the block
    for (j = pi * _pb ; j * _b < _n ; ++j)
    {
      int blockc = _C.Read(j) ;
      if (onecntSum + blockc >= i)
      {
        // the desired 1 is in this block
        size_t blocko = Utils::BitsRead(_O, os, os + _L[blockc] - 1) ; 
        WORD bits = DecodeBits(blockc, blocko) ; 

        int l ;
        for (l = 0 ; l < _b ; ++l)
          if ((bits >> l) & 1)
          {
            ++onecntSum ;
            if (onecntSum == i)
              return j * _b + l ;
          }
        break ;
      }
      os += _L[blockc] ; 
      onecntSum += blockc ;
    }
    return 0 ;
  }

  size_t Select(int type, size_t i) const
  {
    return 0 ;
  }

  size_t GetSpace()
  {
    return _space + sizeof(*this) ;
  }
} ;
}

#endif

#ifndef _MOURISL_COMPACTDS_SEQUENCE_HOMOPOLYMER
#define _MOURISL_COMPACTDS_SEQUENCE_HOMOPOLYMER

#include <math.h>

#include "Sequence.hpp"
#include "Sequence_WaveletTree.hpp"

// Split the original sequence into fixed-length blocks,
//   compress the single-run block by reducing it to one character
namespace compactds {
class Sequence_RunBlock: public Sequence
{
private:
  size_t _b ; // block size
  size_t _blockCnt ;
  Bitvector_Plain _useRunBlock ; // 0-plain sequence, 1-homo polymer sequence 
  //size_t **_alphabetBlockPartialSum ;
  Sequence_WaveletTree<Bitvector_Plain> _waveletSeq ;
  Sequence_WaveletTree<Bitvector_Plain> _runBlockSeq ;

  // Variables and functions related to automatic block size estimation
  size_t _blockSizeInferLength ; // use this amount of numbers to infer block size
  size_t EstimateSpace(const FixedSizeElemArray &S, size_t n, size_t b, int alphabetBit)
  {
    size_t i, j ;
    size_t runBlockCnt = 0 ;
    size_t runBlockLen = 0 ;
    for (i = 0 ; i < n ; i += b)
    {
      uint64_t c = S.Read(i) ;
      bool runBlockFlag = true ;
      for (j = i + 1 ; j < i + b && j < n ; ++j)
      {
        if (S.Read(j) != c)
        {
          runBlockFlag = false ;
          break ;
        }
      }
      if (runBlockFlag)
      {
        ++runBlockCnt ;
        runBlockLen += (j - i) ;
      }
    }
    return DIV_CEIL(n, b) + alphabetBit * (runBlockCnt + n - runBlockLen) ; 
  }

  // Use the first m characters from S to determine block size
  size_t ComputeBlockSize(const FixedSizeElemArray &S, size_t n, size_t alphabetSize)
  {
    size_t i ;
    int alphabetBit = Utils::Log2Ceil(alphabetSize) ;

    size_t bestSpace = 0 ;
    size_t bestTag = 0 ;
    size_t m = (n < _blockSizeInferLength ? n : _blockSizeInferLength) ;
    for (i = 2 ; i <= m ; i *= 2)
    {
      size_t space = EstimateSpace(S, m, i, alphabetBit) ;
      if (bestSpace == 0 || space < bestSpace)
      {
        bestSpace = space ;
        bestTag = i ;
      }
    }

    if (bestTag <= m)
    {
      size_t space = EstimateSpace(S, m, bestTag / 2 * 3, alphabetBit) ;
      if (space < bestSpace)
      {
        bestSpace = space ;
        bestTag = bestTag / 2 * 3 ;
      }
      
      size_t r = 0 ;
      size_t c = S.Read(0) ;
      for (i = 1 ; i < m ; ++i)
      {
        size_t tmp = S.Read(i) ;
        if (tmp != c)
        {
          ++r ;
          c = tmp ;
        }
      }
      size_t testSize = CEIL(sqrt((double)m/(double)r)) ;
      if (testSize > 2)
      {
        space = EstimateSpace(S, m, testSize, alphabetBit) ;
        if (space < bestSpace)
        {
          bestSpace = space ;
          bestTag = testSize ;
        }
      }
    }
    return bestTag ;
  }

public:
  Sequence_RunBlock() 
  {
    _b = 0 ;
    _blockSizeInferLength = (1<<20) ;
  }

  ~Sequence_RunBlock() 
  {
    Free() ;
  }
  
  void Free()
  {
    if (_n > 0)
    {
      //size_t i ;
      //size_t alphabetSize = _alphabets.GetSize() ;
      //for (i = 0 ; i < alphabetSize ; ++i) 
      //  free(_alphabetBlockPartialSum[i]) ;
      //free(_alphabetBlockPartialSum) ;

      _useRunBlock.Free() ;
      _waveletSeq.Free() ;
      _runBlockSeq.Free() ;

      Sequence::Free() ;
    }
  }
 
  void SetBlockSize(size_t b)
  {
    _b = b ;
  }

  void SetExtraParameter(void *p)
  {
    SetBlockSize((size_t)p) ;
  }

  void SetBlockSizeInferLength(size_t l)
  {
    _blockSizeInferLength = l ;
  }
  
  size_t GetSpace() 
  {
    bool inclusive = true ;
    return _space + _alphabets.GetSpace() - sizeof(_alphabets)
      + (inclusive ? sizeof(*this) : 0) ;
  }

  void Init(const FixedSizeElemArray &S, size_t sequenceLength, const ALPHABET *alphabetMap) 
  {
    size_t i, j ;

    _n = sequenceLength ;
    size_t alphabetSize = _alphabets.GetSize() ;
    if (alphabetSize == 0)
    {
      _alphabets.InitFromList(alphabetMap, strlen(alphabetMap)) ;
      alphabetSize = _alphabets.GetSize() ;
    }
    size_t alphabetBits = Utils::Log2Ceil(alphabetSize) ;
    
    if (_b == 0)
      _b = ComputeBlockSize(S, sequenceLength, alphabetSize) ;
    if (_b == 1)
      _b = _n ;
    _blockCnt = DIV_CEIL(_n, _b) ;
    
    WORD *B = Utils::MallocByBits(_blockCnt) ; // block indicator 

    for (i = 0 ; i < _n ; i += _b)
    {
      int prevc = S.Read(i) ;
      size_t rcnt = 1 ;
      for (j = 1 ; j < _b && i + j < _n ; ++j)
      {
        int c = S.Read(i + j) ;
        if (c != prevc)
        {
          ++rcnt ;
          prevc = c ;
          break ;
        }
      }
      if (rcnt == 1)
      {
        Utils::BitSet(B, i / _b) ;
      }
    }
    _useRunBlock.SetSelectSpeed(DS_SELECT_SPEED_NO) ;
    _useRunBlock.Init(B, _blockCnt) ;
    
    // Split the sequence into two parts
    FixedSizeElemArray tmpS ;
    
    size_t rbCount = _useRunBlock.Rank(1, _blockCnt) ;
    if (rbCount > 0)
      tmpS.Malloc(S.GetElemLength(), MAX(rbCount, _n - _b * (rbCount-1)) + 1) ; // The minus 1 to the rbCount is to handle the case where the last block is a run block and will overcount 
    else
      tmpS.Malloc(S.GetElemLength(), _n) ;

    int k ; // use run block
    for ( k = 0 ; k <= 1 ; ++k)
    {
      size_t size = 0 ; 
      size_t elemPerWord = size_t(WORDBITS/alphabetBits) ; // each word can hold this number of elements  
      WORD w = 0 ; // w holding tempoary elements that will be write into tmpS in chunk
      size_t wElem = 0 ; // How many element w is holding now
      for (i = 0 ; i < _n ; i += _b)
      {
        if (Utils::BitRead(B, i / _b) != k)
          continue ;
        if (k == 0)
        {
          for (j = 0 ; j < _b && i + j < _n ; ++j)
          {
            w |= (S.Read(i + j) << (wElem * alphabetBits)) ;
            ++wElem ;
            if (wElem >= elemPerWord)
            {
              tmpS.PackWrite(size, w, wElem) ;
              size += wElem ;
              
              w = 0 ;
              wElem = 0 ;
            }
          }
        }
        else
        {
          w |= (S.Read(i) << (wElem * alphabetBits)) ;
          ++wElem ;
          if (wElem >= elemPerWord)
          {
            tmpS.PackWrite(size, w, wElem) ;
            size += wElem ;

            w = 0 ;
            wElem = 0 ;
          }
        }
      }

      if (wElem > 0)
      {
        tmpS.PackWrite(size, w, wElem) ;
        size += wElem ;

        w = 0 ;
        wElem = 0 ;
      }
      
      tmpS.SetSize(size) ;
      //printf("%d %d\n", _b, size) ;
      if (k == 0)
      {
        if (size > 0)
        {
          _waveletSeq.SetSelectSpeed( DS_SELECT_SPEED_NO ) ;
          _waveletSeq.Init(tmpS, size, alphabetMap) ;
        }
      }
      else
      {
        if (size > 0)
        {
          _runBlockSeq.SetSelectSpeed( DS_SELECT_SPEED_NO ) ;
          _runBlockSeq.Init(tmpS, size, alphabetMap) ;
        }
      }
    }
    _space += _useRunBlock.GetSpace() - sizeof(_useRunBlock) ;
    _space += _waveletSeq.GetSpace() - sizeof(_waveletSeq) ;
    _space += _runBlockSeq.GetSpace() - sizeof(_runBlockSeq) ;
    //printf("%d %d %d\n", sizeof(*this), sizeof(_waveletSeq), sizeof(_runBlockSeq)) ;

    free(B) ;
  }

  ALPHABET Access(size_t i) const 
  {
    size_t bi = i / _b ;
    int type = _useRunBlock.Access(bi) ;
    if (type == 0)
    {
      size_t r = _useRunBlock.Rank(1, bi) ;
      i -= _b * r ;
      return _waveletSeq.Access(i) ;
    }
    else
    {
      size_t r = _useRunBlock.Rank(0, bi) ;
      i -= _b * r ;
      return _runBlockSeq.Access(i/_b) ;
    }
  }

  size_t Rank(ALPHABET c, size_t i, int inclusive = 1) const
  {
    if (!inclusive)
    {
      if (i == 0)
        return 0 ;
      --i ;
    }

    size_t bi = i / _b ;
    int type = _useRunBlock.Access(bi) ;
    //size_t ranki = _useRunBlock.Rank(type, bi) ;
    size_t ranki = _b < _n ? _useRunBlock.Rank(type, bi) : 1 ;
    size_t otherRanki = (bi + 1) - ranki ;
     
    size_t ret = 0 ;
    if (type == 0)
      ret = _waveletSeq.Rank(c, (ranki - 1) * _b + i % _b) ; // ranki>=1 because bi is of type.
    else
    {
      bool inRun = true ;
      size_t rbRank = _runBlockSeq.RankAndTest(c, ranki - 1, inRun) ; // type==1 makes sure ranki >= 1
      if (inRun) // This makes sure rbRank>=1 at (ranki-1)
        ret = (rbRank - 1) * _b + i % _b + 1;
      else
        ret = rbRank * _b ;
    }

    if (otherRanki == 0)
    {
      return ret ;
    }
    if (type == 0)
      ret += _runBlockSeq.Rank(c, otherRanki - 1) * _b ;
    else
      ret += _waveletSeq.Rank(c, otherRanki * _b - 1) ;

    return ret ;
  }

  size_t Select(ALPHABET c, size_t i) const
  {
    return 0 ;
  }

  void Decompress(FixedSizeElemArray &S)
  {
    S.Free() ;
    
    size_t i, j, k ;
    size_t rbIdx = 0 ;
    size_t alphabetBit = Utils::Log2Ceil(_alphabets.GetSize()) ;
    S.Malloc(alphabetBit, _n) ;
    for (i = 0 ; i < _n ; i += _b)
    {
      size_t elemPerWord = size_t(WORDBITS/alphabetBit) ; // each word can hold this number of elements  
      WORD w = 0 ;
      if (_useRunBlock.Access(i / _b) == 1)
      {
        WORD c = _alphabets.Encode(_runBlockSeq.Access(rbIdx)) ;
        for (j = i ; j < i + _b && j < _n ; )
        {
          w = 0 ;
          for (k = j ; k < j + elemPerWord && k < i + _b && k < _n; ++k)
          {
            w |= (c << ((k - j) * alphabetBit)) ; 
          }
          S.PackWrite(j, w, k - j) ;
          j = k ;
        }
        ++rbIdx ;     
      }
      else
      {
        size_t l = i - _b * rbIdx ; 
        for (j = i ; j < i + _b && j < _n ;)
        {
          w = 0 ;
          for (k = j ; k < j + elemPerWord && k < i + _b && k < _n ; ++k, ++l)
          {
            WORD c = _alphabets.Encode(_waveletSeq.Access(l)) ;
            w |= (c << ((k - j) * alphabetBit)) ;
          }
          S.PackWrite(j, w, k - j) ;
          j = k ;
        }
      }
    }
  }

  void Save(FILE *fp)
  {
    Sequence::Save(fp) ;
    SAVE_VAR(fp, _b) ;
    SAVE_VAR(fp, _blockCnt) ;
    _useRunBlock.Save(fp) ;
    _waveletSeq.Save(fp) ;
    _runBlockSeq.Save(fp) ;
  }

  void Load(FILE *fp)
  {
    Free() ;
    
    Sequence::Load(fp) ;
    LOAD_VAR(fp, _b) ;
    LOAD_VAR(fp, _blockCnt) ;
    _useRunBlock.Load(fp) ;
    _waveletSeq.Load(fp) ;
    _runBlockSeq.Load(fp) ;
  }
  
  void PrintStats()
  {
    Utils::PrintLog("Sequence_RunBlock: total_length: %lu block_size: %lu runBlock_block: %lu",
        _n, _b, _useRunBlock.Rank(1, _blockCnt - 1)) ;
    _runBlockSeq.PrintStats() ;
    _waveletSeq.PrintStats() ;
  }
} ;
}

#endif

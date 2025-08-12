#ifndef _MOURISL_COMPACTDS_SEQUENCE_HYBRID
#define _MOURISL_COMPACTDS_SEQUENCE_HYBRID

#include "Sequence.hpp"
#include "Sequence_WaveletTree.hpp"
#include "Sequence_RunLength.hpp"

namespace compactds {
class Sequence_Hybrid: public Sequence
{
private:
  size_t _b ; // block size
  size_t _blockCnt ;
  size_t _minAvgRunLength ; // minimum average run length in a block
  Bitvector_Plain _useRunLength ; // 0-plain sequence, 1-run length sequence 
  //size_t **_alphabetBlockPartialSum ;
  Sequence_WaveletTree<Bitvector_Plain> _waveletSeq ;
  Sequence_RunLength _runlengthSeq ;
  
  size_t _blockSizeInferLength ; // use this amount of numbers to infer block size

  size_t EstimateSpace(const FixedSizeElemArray &S, size_t n, size_t b, size_t minRl, int alphabetBit)
  {
    size_t i, j ;
    size_t rlBlockCnt = 0 ; // the number of blocks for run-length representation
    size_t rlBlockLen = 0 ;
    size_t runCnt = 0 ; // run count in runlength-endcoed sequence.
    size_t lastRunChr = 0 ;
    for (i = 0 ; i < n ; i += b)
    {
      uint64_t c = S.Read(i) ;
      size_t localRunCnt = 1 ;
      for (j = i + 1 ; j < i + b && j < n ; ++j)
      {
        if (S.Read(j) != c)
        {
          ++localRunCnt ;
          c = S.Read(j) ;
        }
      }
      if ((j - i) / localRunCnt >= minRl)
      {
        size_t reduce = 0 ;
        if (S.Read(i) == lastRunChr) 
          reduce = 1 ;
        runCnt += localRunCnt - reduce ; 
        rlBlockLen += (j - i) ; 
        lastRunChr = c ;
        ++rlBlockCnt ;
      }
    }
    size_t ret = DIV_CEIL(n, b) + alphabetBit * (n - rlBlockLen) ;
    
    if (runCnt > 0)
      ret += runCnt * Utils::Log2Ceil(n / runCnt) + alphabetBit * runCnt + runCnt * Utils::Log2Ceil(n * 4 / runCnt) ;

    return ret ;
  }

  // Use the first m characters from S to determine the best block size
  //   the blocksize shall minimize the block bit overhead and 
  //   maximize the number of characters that are in the rl-block
  size_t ComputeBlockSize(const FixedSizeElemArray &S, size_t n, size_t alphabetSize)
  {
    size_t i ;
    int alphabetBit = Utils::Log2Ceil(alphabetSize) ;

    size_t bestSpace = 0 ;
    size_t bestTag = 0 ;
    size_t m = (n < _blockSizeInferLength ? n : _blockSizeInferLength) ;
    for (i = 4 ; i <= m ; i *= 2)
    {
      size_t space = EstimateSpace(S, m, i, _minAvgRunLength, alphabetBit) ;
      if (bestSpace == 0 || space < bestSpace)
      {
        bestSpace = space ;
        bestTag = i ;
      }
    }

    if (bestTag <= m)
    {
      size_t space = EstimateSpace(S, m, bestTag / 2 * 3, _minAvgRunLength, alphabetBit) ;
      if (space < bestSpace)
      {
        bestSpace = space ;
        bestTag = bestTag / 2 * 3 ;
      }
    }
    return bestTag ;
  }

public:
  Sequence_Hybrid() 
  {
    _b = 0 ;
    _minAvgRunLength = 6 ;
    _blockSizeInferLength = (1<<20) ;
  }

  ~Sequence_Hybrid() 
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
      _useRunLength.Free() ;
      _waveletSeq.Free() ;
      _runlengthSeq.Free() ;
      Sequence::Free() ;
    }
  }
 
  void SetBlockSize(size_t b)
  {
    _b = b ;
  }

  void SetblockSizeInferLength(size_t l)
  {
    _blockSizeInferLength = l ;
  }

  void SetMinAvgRunLength(size_t r)
  {
    _minAvgRunLength = r ;
  }

  size_t GetSpace() 
  {
    return _space + sizeof(*this) + _alphabets.GetSpace() - sizeof(_alphabets);
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
    //size_t *psums ; // use this to avoid access the rank in another type of array
    //psums = (size_t *)calloc(alphabetSize, sizeof(size_t)) ;
    
    if (_b == 0)
      _b = ComputeBlockSize(S, sequenceLength, alphabetSize) ;

    _blockCnt = DIV_CEIL(_n, _b) ;
    size_t runlengthBlockCnt = 0 ; 
    
    WORD *B = Utils::MallocByBits(_blockCnt) ; // block indicator 
    //_alphabetBlockPartialSum = (size_t **)malloc(sizeof(size_t *) * alphabetSize) ;
    //_space += sizeof(size_t *) * alphabetSize ;
    /*for (i = 0 ; i < alphabetSize ; ++i)
    {
      //_alphabetBlockPartialSum[i] = (size_t *)malloc(sizeof(size_t) * (_blockCnt + 1)) ;
      //_space += sizeof(size_t) * (_blockCnt + 1) ;
    }*/

    for (i = 0 ; i < _n ; i += _b)
    {
      //for (j = 0 ; j < alphabetSize ; ++j)
      //  _alphabetBlockPartialSum[j][i / _b] = psums[j] ;
      
      int prevc = S.Read(i) ;
      //++psums[prevc] ;
      size_t rcnt = 1 ;
      for (j = 1 ; j < _b && i + j < _n ; ++j)
      {
        int c = S.Read(i + j) ;
        //++psums[c] ;
        if (c != prevc)
        {
          ++rcnt ;
          prevc = c ;
        }
      }
      if (_b / rcnt >= _minAvgRunLength)
      {
        ++runlengthBlockCnt ;
        Utils::BitSet(B, i / _b) ;
      }
    }
    //for (j = 0 ; j < alphabetSize ; ++j)
    //  _alphabetBlockPartialSum[j][i / _b] = psums[j] ;
    _useRunLength.Init(B, _blockCnt) ;
    
    // Split the sequence into two parts
    FixedSizeElemArray tmpS ;
    tmpS.Malloc(S.GetElemLength(), _n) ;
    int k ; // use run length
    for ( k = 0 ; k <= 1 ; ++k)
    {
      size_t size = 0 ; 
      for (i = 0 ; i < _n ; i += _b)
      {
        if (Utils::BitRead(B, i / _b) != k)
          continue ;
        for (j = 0 ; j < _b && i + j < _n ; ++j)
        {
          tmpS.Write(size, S.Read(i + j)) ;
          ++size ;
        }
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
          _runlengthSeq.Init(tmpS, size, alphabetMap) ;
      }
    }
    _space += _useRunLength.GetSpace() - sizeof(_useRunLength) ;
    _space += _waveletSeq.GetSpace() - sizeof(_waveletSeq) ;
    _space += _runlengthSeq.GetSpace() - sizeof(_runlengthSeq) ;
    //printf("%d %d %d\n", sizeof(*this), sizeof(_waveletSeq), sizeof(_runlengthSeq)) ;

    //free(psums) ;
    free(B) ;
  }

  ALPHABET Access(size_t i) const 
  {
    size_t bi = i / _b ;
    int type = _useRunLength.Access(bi) ;
    if (type == 0)
    {
      size_t r = _useRunLength.Rank(1, bi) ;
      i -= _b * r ;
      return _waveletSeq.Access(i) ;
    }
    else
    {
      size_t r = _useRunLength.Rank(0, bi) ;
      i -= _b * r ;
      return _runlengthSeq.Access(i) ;
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
    int type = _useRunLength.Access(bi) ;
    size_t ranki = _useRunLength.Rank(type, bi) ;
    size_t otherRanki = (bi + 1) - ranki ;
     
    size_t ret = 0 ;
    size_t typei = (ranki - 1) * _b + i % _b ; // ranki>=1 because bi is of type. 
    if (type == 0)
      ret = _waveletSeq.Rank(c, typei) ;
    else
      ret = _runlengthSeq.Rank(c, typei) ;
    if (otherRanki == 0)
      return ret ;
    
    size_t otheri = otherRanki * _b - 1 ;
    if (type == 0)
      ret += _runlengthSeq.Rank(c, otheri) ;
    else
      ret += _waveletSeq.Rank(c, otheri) ;

    return ret ;
  }

  size_t Select(ALPHABET c, size_t i) const
  {
    return 0 ;
  }

  void Save(FILE *fp)
  {
    Sequence::Save(fp) ;
    SAVE_VAR(fp, _b) ;
    SAVE_VAR(fp, _blockCnt) ;
    SAVE_VAR(fp, _minAvgRunLength) ;
    _useRunLength.Save(fp) ;
    _waveletSeq.Save(fp) ;
    _runlengthSeq.Save(fp) ;
  }

  void Load(FILE *fp)
  {
    Free() ;
    
    Sequence::Load(fp) ;
    LOAD_VAR(fp, _b) ;
    LOAD_VAR(fp, _blockCnt) ;
    LOAD_VAR(fp, _minAvgRunLength) ;
    _useRunLength.Load(fp) ;
    _waveletSeq.Load(fp) ;
    _runlengthSeq.Load(fp) ;
  }
  
  void PrintStats()
  {
    Utils::PrintLog("Sequence_Hybrid: total_length: %lu block_size: %lu min_avg_runlength: %lu runlength_block: %lu",
        _n, _b, _minAvgRunLength, _useRunLength.Rank(1, _blockCnt - 1)) ;
    _runlengthSeq.PrintStats() ;
    _waveletSeq.PrintStats() ;
  }
} ;
}

#endif

#ifndef _MOURISL_COMPACTDS_SEQUENCE_RUNBLOCKONETREE
#define _MOURISL_COMPACTDS_SEQUENCE_RUNBLOCKONETREE

#include <math.h>

#include "Sequence.hpp"
#include "Sequence_WaveletTree.hpp"

// Split the original sequence into fixed-length blocks,
//   compress the single-run block by reducing it to one character
// This is implemented using a single-wavelet tree. It uses O(m) more bits than a two-wavelet tree implementation.
namespace compactds {
class Sequence_RunBlockOneTree: public Sequence
{
private:
  size_t _b ; // block size
  size_t _blockCnt ;
  Bitvector_Plain _useRunBlock ; // 0-plain sequence, 1-runblock
  Bitvector_Plain *_alphabetRB ; // 0-plain sequence, 1-runblock. for each alphabet
  //size_t **_alphabetBlockPartialSum ;
  Sequence_WaveletTree<Bitvector_Plain> _compressedSeq ;

  // Variables and functions related to automatic block size estimation
  size_t _blockSizeInferLength ; // use this amount of numbers to infer block size

  // Get the total run block compressible length for S[s..e] (inclusive) with block size b
  size_t GetRunBlockLength(const FixedSizeElemArray &S, size_t n, size_t s, size_t e, size_t b)
  {
    size_t i, j ;
    size_t runBlockLen = 0 ;
    for (i = s ; i <= e && i < n ; i += b)
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
        runBlockLen += (j - i) ;
      }
    }

    return runBlockLen ;
  }

  size_t EstimateSpace(const FixedSizeElemArray &S, size_t n, size_t b, int alphabetBit)
  {
    size_t runBlockLen = 0 ;
    size_t len = _blockSizeInferLength ; 
    size_t testCases = 1024 ; 
    size_t m = 0 ; // Actual number of bases used for estimate space
    
    if (len * testCases >= n)
    {
      runBlockLen = GetRunBlockLength(S, n, 0, n - 1, b) ;
      m = n ;
    }
    else
    {
      size_t i ;
      m = 0 ;
      for (i = 0 ; i < n ; i += DIV_CEIL(n, testCases))
      {
        size_t e = i + len - 1 ;
        if (e >= n)
          e = n - 1 ;
        runBlockLen += GetRunBlockLength(S, n, i, i + len - 1, b) ;
        m += e - i + 1 ;
      }
    }
    size_t runBlockCnt = DIV_CEIL(runBlockLen, b) ;
    if (b > 1)
      return DIV_CEIL(m, b) + (alphabetBit + 1) * (runBlockCnt + m - runBlockLen) ; 
    else
      return alphabetBit * m ; // b==1 will be handled specially
  }

  // Use the same chunks for estimating space to estimate the average run length 
  double EstimateAverageRunLength(const FixedSizeElemArray &S, size_t n)
  {
    size_t len = _blockSizeInferLength ; 
    size_t testCases = 1024 ; 
    size_t m = 0 ; // Actual number of bases used for estimate space

    size_t r = 0 ;
    if (len * testCases >= n)
    {
      size_t i ;
      size_t c = S.Read(0) ;
      for (i = 1 ; i < n ; ++i)
      {
        size_t tmp = S.Read(i) ;
        if (tmp != c)
        {
          ++r ;
          c = tmp ;
        }
      }
      ++r ;
      m = n ;
    }
    else
    {
      size_t i, j ;
      m = 0 ;
      for (i = 0 ; i < n ; i += DIV_CEIL(n, testCases))
      {
        size_t e = i + len - 1 ;
        if (e >= n)
          e = n - 1 ;
        
        size_t c = S.Read(i) ;
        for (j = i + 1 ; j <= e ; ++j)
        {
          size_t tmp = S.Read(j) ;
          if (tmp != c)
          {
            ++r ;
            c = tmp ;
          }
        }
        ++r ;
        m += e - i + 1 ;
      }
    }
    return (double)m / (double)r ;
  }

  // Use the characters from dispersed chunks of size m of S to determine block size
  size_t ComputeBlockSize(const FixedSizeElemArray &S, size_t n, size_t alphabetSize)
  {
    size_t i ;
    int alphabetBit = Utils::Log2Ceil(alphabetSize) ;

    size_t bestSpace = 0 ;
    size_t bestTag = 0 ;
    size_t m = _blockSizeInferLength ;
    for (i = 1 ; i <= m ; i *= 2)
    {
      size_t space = EstimateSpace(S, n, i, alphabetBit) ;
      if (bestSpace == 0 || space < bestSpace)
      {
        bestSpace = space ;
        bestTag = i ;
      }
    }

    if (bestTag <= m)
    {
      size_t space = 0 ;
      if (bestTag >= 2)
      {
        space = EstimateSpace(S, n, bestTag / 2 * 3, alphabetBit) ;
        if (space < bestSpace)
        {
          bestSpace = space ;
          bestTag = bestTag / 2 * 3 ;
        }
      }
      
      size_t testSize = CEIL(sqrt(EstimateAverageRunLength(S, n))) ;
      if (testSize > 2)
      {
        space = EstimateSpace(S, n, testSize, alphabetBit) ;
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
  Sequence_RunBlockOneTree() 
  {
    _b = 0 ;
    _blockSizeInferLength = (1<<10) ; // Up to 1024 by default
  }

  ~Sequence_RunBlockOneTree() 
  {
    Free() ;
  }
  
  void Free()
  {
    if (_n > 0)
    {
      delete[] _alphabetRB ;
      _useRunBlock.Free() ;
      _compressedSeq.Free() ;
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
    size_t *alphabetCnt = (size_t *)calloc(alphabetSize, sizeof(size_t)) ;

    if (_b != _n)
    {
      WORD *blockContent = (WORD *)calloc(_b, sizeof(WORD)) ;
      for (i = 0 ; i < _n ; i += _b)
      {
        uint64_t prevc = S.Read(i) ;
        blockContent[0] = prevc ;
        size_t rcnt = 1 ;
        for (j = 1 ; j < _b && i + j < _n ; ++j)
        {
          uint64_t c = S.Read(i + j) ;
          blockContent[j] = c ;
          if (c != prevc)
          {
            ++rcnt ;
            prevc = c ;
          }
        }
        if (rcnt == 1)
        {
          ++alphabetCnt[prevc] ;
          Utils::BitSet(B, i / _b) ;
        }
        else 
        {
          for (j = 0 ; j < _b && i + j < _n ; ++j)
            ++alphabetCnt[blockContent[j]] ;
        }
      }
      free(blockContent) ;
    }
    else if (_n > 0)//_b==_n
    {
      uint64_t c = S.Read(0) ;
      for (i = 1 ; i < _n ; ++i)
      {
        if (S.Read(i) != c)
          break ;
      }

      if (i == _n)
      {
        alphabetCnt[c] = _n ;
        Utils::BitSet(B, 0) ;
      } 
      // Else: when there is no run block. we can set alphabetRB as empty 
      //   and let rank9 handle the special case
    }

    _useRunBlock.SetSelectSpeed(DS_SELECT_SPEED_NO) ;
    _useRunBlock.Init(B, _blockCnt) ;
    _alphabetRB = new Bitvector_Plain[alphabetSize] ;
    for (i = 0 ; i < alphabetSize ; ++i)
    {
      _alphabetRB[i].SetSelectSpeed(DS_SELECT_SPEED_NO) ;
      _alphabetRB[i].Malloc(alphabetCnt[i]) ;
    }

    // Create a new sequence that is mixed of runblock comrpessed and uncompressed
    FixedSizeElemArray tmpS ;
    tmpS.Malloc(S.GetElemLength(), _n) ;
    size_t size = 0 ; 
    size_t elemPerWord = size_t(WORDBITS/alphabetBits) ; // each word can hold this number of elements  
    WORD w = 0 ; // w holding tempoary elements that will be write into tmpS in chunk
    size_t wElem = 0 ; // How many element w is holding now
    for (i = 0 ; i < alphabetSize ; ++i)
      alphabetCnt[i] = 0 ;

    for (i = 0 ; i < _n ; i += _b)
    {
      if (Utils::BitRead(B, i / _b) == 0) // a non-run-block 
      {
        for (j = 0 ; j < _b && i + j < _n ; ++j)
        {
          uint64_t c = S.Read(i+j) ;
          w |= (c << (wElem * alphabetBits)) ;
          ++alphabetCnt[c] ;
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
      else // This is a run-block
      {
        uint64_t c = S.Read(i) ;
        w |= (c << (wElem * alphabetBits)) ;
        _alphabetRB[c].BitSet(alphabetCnt[c]) ;
        ++alphabetCnt[c] ;
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
    _compressedSeq.SetSelectSpeed(DS_SELECT_SPEED_NO) ;
    _compressedSeq.Init(tmpS, size, alphabetMap) ;
    for (i = 0 ; i < alphabetSize ; ++i)
    {
      _alphabetRB[i].Init() ;
    }

    _space = _useRunBlock.GetSpace() - sizeof(_useRunBlock) ;
    for (i = 0 ; i < alphabetSize ; ++i)
      _space += _alphabetRB[i].GetSpace() - sizeof(_alphabetRB[i]) ;
    _space += _compressedSeq.GetSpace() - sizeof(_compressedSeq) ;
    //printf("%d %d %d\n", sizeof(*this), sizeof(_waveletSeq), sizeof(_runBlockSeq)) ;
    
    free(alphabetCnt) ;
    free(B) ;
  }

  ALPHABET Access(size_t i) const 
  {
    size_t bi = i / _b ;
    size_t r = _useRunBlock.Rank(1, bi, /*inclusive=*/0) ;
    int type = _useRunBlock.Access(bi) ;
    
    // remap i to the compressed string position
    if (type == 1) 
    {
      // if it is in a run block, then move to the "representative" position
      i -= (i % _b) ; 
    }
    i -= (_b - 1) * r ; // Each runblock before introduce _b-1 extra positions
    return _compressedSeq.Access(i) ;
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
    size_t r = _useRunBlock.Rank(1, bi, /*inclusive=*/0) ;
    int type = _useRunBlock.Access(bi) ; // is the block containing i a run block or not.
    
    // remap i to the compressed string position
    size_t remainder = 0 ;
    size_t ci = i ; // index on the compressed seq
    if (type == 1) 
    {
      // if it is in a run block, then move to the "representative" position
      remainder = i % _b; 
      ci -= (i % _b) ; 
    }
    ci -= (_b - 1) * r ; // Each runblock before introduce _b-1 extra positions
    
    size_t ret = 0 ;
    bool inRun = false ;
    if (type == 1)
      ret = _compressedSeq.RankAndTest(c, ci, inRun) ;
    else
      ret = _compressedSeq.Rank(c, ci) ;
    if (ret > 0 && (_b < _n || inRun) )// The second condition will handle the case when we _b==_n and everything is compressed.
      // ret is 1-based rank value, so need to -1 to get the index in the _alphabetRB
      ret += _alphabetRB[ _alphabets.Encode(c) ].Rank(1, ret - 1) * (_b - 1) ; 
    if (inRun) // inRun can only be true when type==1
      ret = ret - (_b - 1) + remainder ;
    //printf("%c %lu %lu %lu %lu. %lu %d %d\n", c, _alphabets.Encode(c), i, ci, ret, remainder, type, inRun) ;
    return ret ;
  }

  size_t Select(ALPHABET c, size_t i) const
  {
    return 0 ;
  }

  void Decompress(FixedSizeElemArray &S)
  {
    S.Free() ;
    
    /*size_t i, j, k ;
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
    }*/
  }

  void Save(FILE *fp)
  {
    Sequence::Save(fp) ;
    SAVE_VAR(fp, _b) ;
    SAVE_VAR(fp, _blockCnt) ;
    _useRunBlock.Save(fp) ;
    size_t alphabetSize = _alphabets.GetSize() ;
    size_t i ;
    for (i = 0 ; i < alphabetSize ; ++i)
      _alphabetRB[i].Save(fp) ;
    _compressedSeq.Save(fp) ;
  }

  void Load(FILE *fp)
  {
    Free() ;
    
    Sequence::Load(fp) ;
    LOAD_VAR(fp, _b) ;
    LOAD_VAR(fp, _blockCnt) ;
    _useRunBlock.Load(fp) ;
    size_t alphabetSize = _alphabets.GetSize() ;
    size_t i ;
    _alphabetRB = new Bitvector_Plain[alphabetSize] ;
    for (i = 0 ; i < alphabetSize ; ++i)
      _alphabetRB[i].Load(fp) ;
    _compressedSeq.Load(fp) ;
  }
  
  void PrintStats()
  {
    Utils::PrintLog("Sequence_RunBlockOneTree: total_length: %lu block_size: %lu runBlock_block: %lu",
        _n, _b, _useRunBlock.Rank(1, _blockCnt - 1)) ;
    _compressedSeq.PrintStats() ;
  }
} ;
}

#endif

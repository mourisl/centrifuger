#ifndef _MOURISL_COMPACTDS_SEQUENCE_RUNLENGTH
#define _MOURISL_COMPACTDS_SEQUENCE_RUNLENGTH

#include "Sequence.hpp"
#include "Sequence_WaveletTree.hpp"
#include "PartialSum.hpp"

// This sequence type assumes the alphabet coder is plain.
namespace compactds {
class Sequence_RunLength : public Sequence
{
private:
  Bitvector_Sparse _runs ; // mark the beginning of each runs, E in the manuscript
  Sequence_WaveletTree<Bitvector_Plain> _runChars ; // the character for each run, supporting ranking, L' in the manuscript
  PartialSum *_alphabetPartialSum ; // the partial length with respect to each alphabet, D in the manuscript
  size_t _rcnt ;
public:
  Sequence_RunLength() 
  {
    _n = _rcnt = 0 ;
  }

  ~Sequence_RunLength() 
  {
    Free() ;
  }

  void Free()
  {
    if (_n > 0)
    {
      _rcnt = 0 ;
      delete[] _alphabetPartialSum  ;
      _runChars.Free() ;
      _runs.Free() ;
      Sequence::Free() ;
    }
  }

  size_t GetSpace()
  {
    return _space + _alphabets.GetSize() - sizeof(_alphabets) + sizeof(*this) ;
  }
  
  void Init(const FixedSizeElemArray &S, size_t sequenceLength, const ALPHABET *alphabetMap) 
  {
    size_t i ;
    uint8_t c ; // character
  
    if (_alphabets.GetSize() == 0)
      _alphabets.InitFromList(alphabetMap, strlen(alphabetMap)) ;

    // Get the runs
    _n = sequenceLength ;
    _rcnt = 1 ;
    c = S.Read(0) ;
    for (i = 1 ; i < sequenceLength ; ++i)
    {
      if ( S.Read(i) != c)
      {
        ++_rcnt ;
        c = S.Read(i) ;
      }
    }
    
    FixedSizeElemArray chars ;
    WORD *W = Utils::MallocByBits(sequenceLength + 2) ;
    chars.Malloc(S.GetElemLength(), _rcnt) ;
    
    c = S.Read(0) ;
    chars.Write(0, c) ;
    Utils::BitSet(W, 0) ;
    _rcnt = 1 ;
    for (i = 1 ; i < sequenceLength ; ++i)
    {
      if (S.Read(i) != c)
      {
        c = S.Read(i) ;
        
        Utils::BitSet(W, i) ;
        chars.Write(_rcnt, c) ;
        ++_rcnt ;
      }
    }
    //Utils::BitSet(W, sequenceLength) ;
    _runs.Init(W, sequenceLength + 2) ;
    _runChars.SetSelectSpeed( DS_SELECT_SPEED_NO ) ;
    _runChars.Init(chars, _rcnt, alphabetMap) ;
    _space = _runs.GetSpace() - sizeof(_runs) + _runChars.GetSpace() - sizeof(_runChars) ;
    
    
    // Process the runs/partial sums for each alphabet
    int alphabetSize = _alphabets.GetSize() ;
    _alphabetPartialSum = new PartialSum[alphabetSize] ;
    for (c = 0 ; c < alphabetSize ; ++c)
    {
      memset(W, 0, Utils::BitsToWords(sequenceLength) * sizeof(WORD)) ;
      size_t psum = 0 ;
      Utils::BitSet(W, 0) ;
      for (i = 0 ; i < sequenceLength ; )
      {
        if (S.Read(i) != c)
        {
          ++i ;
          continue ;
        }
        size_t j ;
        for (j = i ; j < sequenceLength ; ++j)
          if (S.Read(j) != c)
            break ;
        psum += j - i ;
        Utils::BitSet(W, psum) ;
        i = j ;
      }
      _alphabetPartialSum[c].InitFromBitvector(W, psum + 1) ;
      _space += _alphabetPartialSum[c].GetSpace() - sizeof(_alphabetPartialSum[c]) ;
    }
    free(W) ;
  }

  ALPHABET Access(size_t i) const 
  {
    return _runChars.Access(_runs.Rank(1, i) - 1) ;
  }

  size_t Rank(ALPHABET c, size_t i, int inclusive = 1) const  
  {
    if (!inclusive)
    {
      if (i == 0)
        return 0 ;
      else
        --i ;
    }

    size_t cid = _alphabets.Encode(c) ;
    size_t rrank = _runs.Rank(1, i) ; // rank in runs
    
    bool inRun = true ;
    size_t crank = _runChars.RankAndTest(c, rrank - 1, inRun) ; // rank for this character
    //printf("%c %d: rrank=%d crank=%d\n", c, i, rrank, crank) ; 
    if (inRun) 
    {
      size_t psum = _alphabetPartialSum[cid].Sum(crank - 1) ;
      //printf("%d %d. ret=%d\n", psum, _runs.Select(1, rrank),
      //    psum + i - _runs.Select(1, rrank) + 1) ;
      return psum + i - _runs.Select(1, rrank) + 1 ;
    }
    else
    {
      //printf("other %d\n", _alphabetPartialSum[cid].Sum(crank)) ;
      return _alphabetPartialSum[cid].Sum(crank) ;
    }
  }

  // Not supported
  size_t Select(ALPHABET c, size_t i) const
  {
    return 0 ;
  }

  void Save(FILE *fp)
  {
    Sequence::Save(fp) ;  
    SAVE_VAR(fp, _rcnt) ;
    _runs.Save(fp) ;
    _runChars.Save(fp) ;
    int alphabetSize = _alphabets.GetSize() ;
    for (int i = 0 ; i < alphabetSize ; ++i)
      _alphabetPartialSum[i].Save(fp) ;
  }

  void Load(FILE *fp)
  {
    Free() ;
    Sequence::Load(fp) ;
    LOAD_VAR(fp, _rcnt) ;
    _runs.Load(fp) ;
    _runChars.Load(fp) ;
    int alphabetSize = _alphabets.GetSize() ;
    _alphabetPartialSum = new PartialSum[alphabetSize] ;
    for (int i = 0 ; i < alphabetSize ; ++i)
      _alphabetPartialSum[i].Load(fp) ;
  }

  void PrintStats()
  {
    Utils::PrintLog("Sequence_RunLength: total_length: %lu run count: %lu", _n, _rcnt) ;
  }
} ;
}

#endif 

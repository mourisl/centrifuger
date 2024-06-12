#ifndef _MOURISL_COMPACTDS_FM_INDEX
#define _MOURISL_COMPACTDS_FM_INDEX

#include <stdio.h>

#include "Alphabet.hpp"
#include "FixedSizeElemArray.hpp"
#include "FMBuilder.hpp"

// Auxiliary data, other than the BWT and F (alphabet partial sum), for FM index
// Should be directly initalized through FMBuilderParam, simplifies the parameter passing
namespace compactds {
struct _FMIndexAuxData 
{
  size_t n ; // the length of the text

  int sampleStrategy ;
  int sampleRate ;
  size_t sampleSize ;
  FixedSizeElemArray sampledSA ;

  // precomputedRange: the BWT range for a prefix of size param.precomputeWidth
  //                  The pair format is (the start position, and the length of the range).
  //                  The advantage is that we can easily tell whether a range is empty.
  size_t precomputeWidth ;
  size_t precomputeSize ;
  std::pair<size_t, size_t> *precomputedRange ;

  size_t maxLcp ; // only consider LCP up to this point
  WORD *semiLcpGreater ; // The LCP is between current suffix and its previous one
  WORD *semiLcpEqual ;

  size_t adjustedSA0 ;
  std::map<size_t, size_t> selectedSA ; // SAs for speical purposes: e.g. boundary of genomes 
  WORD *selectedSAFilter ; // Quick test whether a SA could be selectedSA 
  int selectedSAFilterSampleRate ;

  bool printLog ;

  _FMIndexAuxData()
  {
    sampleStrategy = 0 ;
    sampleRate = 0 ;
    sampleSize = 0 ;
    precomputeWidth = 0 ;
    precomputeSize = 0 ;
    precomputedRange = NULL ;
    
    maxLcp = 0 ;
    semiLcpGreater = NULL ;
    semiLcpEqual = NULL ;

    adjustedSA0 = 0 ;
    selectedSAFilter = NULL ;
    selectedSAFilterSampleRate = 1024 ;

    printLog = true ;
  }

  ~_FMIndexAuxData()
  {
    // NOTE: has to be explicitly called through Free to release the memory.
  } ;

  void Free()
  {
    sampledSA.Free() ;
    
    if (precomputedRange)
    {
      free(precomputedRange) ;
      precomputedRange = NULL ;
    }

    if (semiLcpGreater)
    {
      free(semiLcpGreater) ;
      free(semiLcpEqual) ;
      semiLcpGreater = NULL ;
      semiLcpEqual = NULL ;
    }

    if (selectedSA.size() > 0)
    {
      selectedSA.clear() ;
      free(selectedSAFilter) ;
    }
  }

  void Save(FILE *fp) 
  {
    SAVE_VAR(fp, n) ;
    SAVE_VAR(fp, sampleStrategy) ;
    SAVE_VAR(fp, sampleRate) ;
    SAVE_VAR(fp, sampleSize) ;
    SAVE_VAR(fp, precomputeWidth) ;
    SAVE_VAR(fp, precomputeSize) ;
    SAVE_VAR(fp, adjustedSA0) ;

    sampledSA.Save(fp) ;
    SAVE_ARR(fp, precomputedRange, precomputeSize) ;

    SAVE_VAR(fp, maxLcp) ;
    if (maxLcp > 0)
    {
      fwrite(semiLcpGreater, sizeof(*semiLcpGreater), Utils::BitsToWords(n), fp) ;
      fwrite(semiLcpEqual, sizeof(*semiLcpEqual), Utils::BitsToWords(n), fp) ;
    }

    // For speical SAs
    size_t tmpSize = selectedSA.size() ;
    SAVE_VAR(fp, tmpSize) ;
    SAVE_VAR(fp, selectedSAFilterSampleRate) ;
    for (std::map<size_t, size_t>::iterator iter = selectedSA.begin() ;
        iter != selectedSA.end() ; ++iter)
    {
      size_t pair[2] = {iter->first, iter->second} ;
      fwrite(pair, sizeof(size_t), 2, fp) ;
    }
  }

  void Load(FILE *fp)
  {
    Free() ;
    size_t i ;

    LOAD_VAR(fp, n) ;
    LOAD_VAR(fp, sampleStrategy) ;
    LOAD_VAR(fp, sampleRate) ;
    LOAD_VAR(fp, sampleSize) ;
    LOAD_VAR(fp, precomputeWidth) ;
    LOAD_VAR(fp, precomputeSize) ;
    LOAD_VAR(fp, adjustedSA0) ;

    sampledSA.Load(fp) ; 
    precomputedRange = (std::pair<size_t, size_t> *)malloc(
        sizeof(std::pair<size_t, size_t>) * precomputeSize) ;
    LOAD_ARR(fp, precomputedRange, precomputeSize) ;

    LOAD_VAR(fp, maxLcp) ;
    if (maxLcp > 0)
    {
      semiLcpGreater = Utils::MallocByBits(n) ;
      semiLcpEqual = Utils::MallocByBits(n) ;
      fread(semiLcpGreater, sizeof(*semiLcpGreater), Utils::BitsToWords(n), fp) ;
      fread(semiLcpEqual, sizeof(*semiLcpEqual), Utils::BitsToWords(n), fp) ;
    }

    size_t tmpSize = 0 ;
    LOAD_VAR(fp, tmpSize) ;
    LOAD_VAR(fp, selectedSAFilterSampleRate) ;
    if (tmpSize > 0)
    {
      selectedSAFilter = Utils::MallocByBits(DIV_CEIL(n, selectedSAFilterSampleRate)) ; 
      for (i = 0 ; i < tmpSize ; ++i)
      {
        size_t pair[2] ;
        fread(pair, sizeof(size_t), 2, fp) ;
        selectedSA[pair[0]] = pair[1] ;
        Utils::BitSet(selectedSAFilter, pair[0] / selectedSAFilterSampleRate) ;
      }
    }
  }
} ;

template <class SeqClass>
class FMIndex
{
private:
  SeqClass _BWT ;
  size_t _n ;
  Alphabet _alphabets ; // May handle more complex mapping, e.g. Huffman coding
  Alphabet _plainAlphabetCoder ; // for plain mapping, important for partial sum access
  size_t *_plainAlphabetPartialSum ;
  size_t _plainAlphabetBits ; // Needed for coding index accessing precomputedRange
  size_t _firstISA ; // ISA[0]
  ALPHABET _lastChr ; // last character in the original text 
  
  // @return: whether SA[i] information is stored
  // the SA information is returned through the reference sa 
  bool GetSampledSA(size_t i, size_t &sa)
  {
    if (i == _firstISA)
    {
      sa = _auxData.adjustedSA0 ;
      return true ;
    }
    else if (i % _auxData.sampleRate == 0)
    {
      sa = _auxData.sampledSA[i / _auxData.sampleRate] ;
      return true ;
    }
    else if (_auxData.selectedSAFilter)
    {
      if (Utils::BitRead(_auxData.selectedSAFilter,  i / _auxData.selectedSAFilterSampleRate)
          && (_auxData.selectedSA.find(i) != _auxData.selectedSA.end()))
      {
        sa = _auxData.selectedSA[i] ;
        return true ;
      }
    }

    return false ;
  }
public:
  struct _FMIndexAuxData _auxData ; // the data used for locate operation
  
  FMIndex() 
  {
    _n = 0 ;
  }
  
  ~FMIndex() 
  {
    Free() ;
  }
  
  void SetAlphabetCode(const Alphabet &a)
  {
    _alphabets = a ;  
  }

  void SetSequenceExtraParameter(void *p)
  {
    _BWT.SetExtraParameter(p) ;
  }

  void Free()
  {
    if (_n > 0)
    {
      _n = 0 ;
      free(_plainAlphabetPartialSum) ;
      _auxData.Free() ;
    }
  }

  void InitAuxData(struct _FMBuilderParam &builderParam)
  {
    _auxData.n = builderParam.n ;
  
    _auxData.sampleRate = builderParam.sampleRate ;
    _auxData.sampleSize = builderParam.sampleSize ;
    _auxData.sampleStrategy = builderParam.sampleStrategy ;
    //_auxData.sampledSA = builderParam.sampledSA ;
    _auxData.sampledSA.InitFromArray(0, builderParam.sampledSA, _auxData.sampleSize) ;
    free(builderParam.sampledSA) ;
    
    _auxData.precomputeWidth = builderParam.precomputeWidth ;
    _auxData.precomputeSize = builderParam.precomputeSize ;
    _auxData.precomputedRange = builderParam.precomputedRange ;

    _auxData.maxLcp = builderParam.maxLcp ;
    _auxData.semiLcpGreater = builderParam.semiLcpGreater ; 
    _auxData.semiLcpEqual = builderParam.semiLcpEqual ;

    _auxData.adjustedSA0 = builderParam.adjustedSA0 ;

    if (builderParam.selectedSA.size() > 0)
    {
      _auxData.selectedSAFilter = Utils::MallocByBits(DIV_CEIL(_auxData.n, 
            _auxData.selectedSAFilterSampleRate)) ; 
      
      _auxData.selectedSA = builderParam.selectedSA ;
      for (std::map<size_t, size_t>::iterator iter = _auxData.selectedSA.begin() ;
                    iter != _auxData.selectedSA.end(); ++iter)
      {
        Utils::BitSet(_auxData.selectedSAFilter, 
            iter->first / _auxData.selectedSAFilterSampleRate) ;
      }
    }
  }

  void Init(FixedSizeElemArray &BWT, size_t n,
    size_t firstISA, struct _FMBuilderParam& builderParam, 
    const ALPHABET *alphabetMapping, int alphabetSize) 
  {
    size_t i ;
    
    _plainAlphabetCoder.InitFromList(alphabetMapping, alphabetSize) ; // The input BWT string should be also plain coded in the same fashion
    _plainAlphabetBits = Utils::Log2Ceil(alphabetSize) ;
    
    if (_alphabets.GetSize() == 0)
      _alphabets.InitFromList(alphabetMapping, alphabetSize) ;
    
    // Auxiliary data structures
    _n = n ;
    _firstISA = firstISA ;
    _lastChr = alphabetMapping[ BWT.Read(firstISA) ] ;
    InitAuxData(builderParam) ;

    // L list
    _BWT.SetAlphabet(_alphabets) ;
    _BWT.Init(BWT, n, alphabetMapping) ;
    if (_auxData.printLog)
      _BWT.PrintStats() ;

    // F list
    _plainAlphabetPartialSum = (size_t *)calloc(alphabetSize + 1,
        sizeof(*_plainAlphabetPartialSum)) ;
    for (i = 0 ; i < n ; ++i)
    {
      ++_plainAlphabetPartialSum[BWT.Read(i)] ; 
    }
    for (i = 1 ; (int)i < alphabetSize ; ++i)
      _plainAlphabetPartialSum[i] += _plainAlphabetPartialSum[i - 1] ;
    for (i = alphabetSize ; i >= 1 ; --i)
      _plainAlphabetPartialSum[i] = _plainAlphabetPartialSum[i - 1] ;
    _plainAlphabetPartialSum[0] = 0 ;
  }

  size_t Rank(ALPHABET c, size_t p, int inclusive = 1)
  {
    size_t ret = _BWT.Rank(c, p, inclusive) ;
    // Since we do not use $, the last character in the original string 
    //   will be moved to the _firstISA instead of the first position
    //   We need to move this back
    // Potential future refactoring: appending an A to the end of the string
    if (c == _lastChr && (p < _firstISA || (!inclusive && p == _firstISA)))
      ++ret ;
    return ret ;
  }

  void BackwardExtend(ALPHABET c, size_t sp, size_t ep, 
      size_t &nextSp, size_t &nextEp)
  {
    size_t offset = _plainAlphabetPartialSum[ _plainAlphabetCoder.Encode(c) ] ; 
    //printf("%c: %d %d %d. %d %d\n", c, offset, sp, ep, _BWT.Rank(c, sp, 0),
    //    _BWT.Rank(c, ep)) ;
    // Need minus 1 here because the return of Rank is 1-based.
    nextSp = offset + Rank(c, sp, /*inclusive=*/0) + 1 - 1 ;
    
    // TODO: Fix a potential issue of underflow.
    //       Now it is handled by out side
    if (sp != ep)
      nextEp = offset + Rank(c, ep) - 1 ;
    else
      nextEp = nextSp + ((_BWT.Access(ep) == c) ? 0 : -1) ;
  }

  // This one is essentially LF mapping 
  size_t BackwardExtend(ALPHABET c, size_t p)
  {
    size_t offset = _plainAlphabetPartialSum[ _plainAlphabetCoder.Encode(c) ] ;
    return offset + Rank(c, p) - 1 ;
  }

  // m - length of s
  // Return the [sp, ep] through the option, and the length of matched prefix in size_t
  size_t BackwardSearch(char *s, size_t m, size_t &sp, size_t &ep)
  {
    size_t i ;
    if (m < _auxData.precomputeWidth)
      return 0 ;

    if (_auxData.precomputeWidth > 0)
    {
      WORD initW = 0 ;
      for (i = 0 ; i < _auxData.precomputeWidth ; ++i)
      {
        if (!_alphabets.IsIn(s[m - 1 - i]))
        {
          sp = 1 ;
          ep = 0 ;
          return i ;
        }
        initW = (initW << _plainAlphabetBits) | (_plainAlphabetCoder.Encode(s[m - 1 - i])) ;
      }
      
      if (_auxData.precomputedRange[initW].second == 0)
      {
        sp = 1 ;
        ep = 0 ;
        return _auxData.precomputeWidth - 1 ;
      }
      sp = _auxData.precomputedRange[initW].first ;
      ep = sp + _auxData.precomputedRange[initW].second - 1 ;
    }
    else
    {
      sp = 0 ;
      ep = _n - 1 ;
    }

    size_t l = _auxData.precomputeWidth ; 
    size_t nextSp = sp ;
    size_t nextEp = ep ;
    while (l < m)
    {
      if (!_alphabets.IsIn(s[m - 1 - l]))
        break ;
      BackwardExtend(s[m - 1 - l], sp, ep, nextSp, nextEp) ;
      if ( nextSp > nextEp || nextEp > _n)
        break ;
      sp = nextSp ;
      ep = nextEp ;
      ++l ;
    }
    return l ;
  }

  // @return: the value of the sampled SA for BWT[i]
  //          l is the offset between 
  size_t BackwardToSampledSA(size_t i, size_t &l)
  {
    l = 0 ;
    size_t ret = 0 ;
    while (!GetSampledSA(i, ret))
    {
      i = BackwardExtend( _BWT.Access(i), i) ;
      ++l ;
    }
    return ret ;
  }

  // return ISA[n - 1]
  size_t GetLastISA()
  {
    return _plainAlphabetPartialSum[ _plainAlphabetCoder.Encode(_lastChr) ] ;
  }

  // Calculate the values for SA[sp..ep]
  void LocateRange(size_t sp, size_t ep, bool withOffset, std::vector<size_t> &locatedSA)
  {
    size_t i ;
    locatedSA.clear() ;
    for (i = sp ; i <= ep ; ++i)
    {
      size_t l ;
      size_t sa = BackwardToSampledSA(i, l) ;
      if (withOffset)
        locatedSA.push_back(sa + l) ;
      else
        locatedSA.push_back(sa) ;
    }
  }

  size_t GetSize()
  {
    return _n ;
  }

  size_t GetAlphabetSize()
  {
    return _alphabets.GetSize() ;
  }

  void PrintSpace()
  {
    Utils::PrintLog("FM-index space usage (bytes):") ;
    Utils::PrintLog("BWT: %llu", _BWT.GetSpace()) ;
    Utils::PrintLog("sampledSA: %llu", _auxData.sampledSA.GetSpace()) ;
    Utils::PrintLog("precomputedRange: %llu", _auxData.precomputeSize * sizeof(*_auxData.precomputedRange)) ;
  }

  void Save(FILE *fp)
  {
    SAVE_VAR(fp, _n) ;
    SAVE_VAR(fp, _plainAlphabetBits) ;
    SAVE_VAR(fp, _firstISA) ;
    SAVE_VAR(fp, _lastChr) ;

    _BWT.Save(fp) ;

    _alphabets.Save(fp) ;
    _plainAlphabetCoder.Save(fp) ;
    size_t alphabetSize = _plainAlphabetCoder.GetSize() ;
    SAVE_ARR(fp, _plainAlphabetPartialSum, alphabetSize + 1) ;
    
    _auxData.Save(fp) ;
  }

  void Load(FILE *fp)
  {
    Free() ;

    LOAD_VAR(fp, _n) ;
    LOAD_VAR(fp, _plainAlphabetBits) ;
    LOAD_VAR(fp, _firstISA) ;
    LOAD_VAR(fp, _lastChr) ;

    _BWT.Load(fp) ;

    _alphabets.Load(fp) ;
    _plainAlphabetCoder.Load(fp) ;
    size_t alphabetSize = _plainAlphabetCoder.GetSize() ;
    _plainAlphabetPartialSum = (size_t *)calloc(alphabetSize + 1,
        sizeof(*_plainAlphabetPartialSum)) ;
    LOAD_ARR(fp, _plainAlphabetPartialSum, alphabetSize + 1) ;

    _auxData.Load(fp) ; 
  }
} ;
}

#endif

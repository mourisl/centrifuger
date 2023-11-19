#ifndef _MOURISL_COMPACTDS_COMPRESSED_SUFFIX_ARRAY
#define _MOURISL_COMPACTDS_COMPRESSED_SUFFIX_ARRAY

#include "Bitvector_Sparse.hpp"
#include "Sequence_WaveletTree.hpp"

namespace compactds {
class CompressedSuffixArray
{
private:
  size_t _space ;
  Bitvector_Sparse *_Psi ; // Psi for each alphabet
  Bitvector_Sparse _D ; // mark the positions of starting alphabet in suffix
  Alphabet _alphabets ; // use plain alphabet set here.
  size_t _n ;
  size_t *_alphabetPartialSum ;
  size_t _firstISA ;
  ALPHABET _lastChr ;  
  WORD **_psiB ; // bits for encoding Psis
  
  size_t Rank(Sequence_WaveletTree<Bitvector_Plain> &BWT, ALPHABET c, size_t p, int inclusive = 1)
  {
    size_t ret = BWT.Rank(c, p, inclusive) ;
    // Since we do not use $, the last character in the original string 
    //   will be moved to the _firstISA instead of the first position
    //   We need to move this back
    // Potential future refactoring: appending an A to the end of the string
    if (c == _lastChr && (p < _firstISA || (!inclusive && p == _firstISA)))
      ++ret ;
    return ret ;
  }

public:
  CompressedSuffixArray() 
  {
    _n = _space = 0 ;
  }
  ~CompressedSuffixArray() {}
  
  void Free()
  {
    if (_n > 0)
    {
      delete[] _Psi ;
    }
    
    if (_psiB != NULL)
    {
      int alphabetSize = _alphabets.GetSize() ;
      int i ;
      for (i = 0 ; i < alphabetSize ; ++i)
        free(_psiB[i]) ;   
      free(_psiB) ;
    }
  }
  
  // Allocate necessary memories for reading SAs 
  void Prepare(size_t n, ALPHABET *alphabetList)
  {
    int i ;

    _n = n ;
    _alphabets.InitFromList(alphabetList, strlen(alphabetList)) ;
    
    int alphabetSize = _alphabets.GetSize() ;
    _Psi = new Bitvector_Sparse[alphabetSize] ;
    _psiB = (WORD **)malloc(sizeof(_psiB[0]) * alphabetSize) ;
    for (i = 0 ; i < alphabetSize ; ++i)
      _psiB[i] = Utils::MallocByBits(n) ;
    _alphabetPartialSum = (size_t *)calloc(alphabetSize + 1, sizeof(size_t)) ;
  }
  
  
  // sa corresponding to SA[from..to], inclusive
  void ReadSaChunk(FixedSizeElemArray &T, size_t n, size_t *sa, size_t from, size_t to)
  {
    size_t i ;

    if (to >= n)
      to = n - 1 ;
    int alphabetSize = _alphabets.GetSize() ;
    size_t *alphabetCount = (size_t *)calloc(alphabetSize + 1, sizeof(size_t)) ;
    for (i = from ; i <= to ; ++i)
    {
      size_t s = sa[i - from] ;
      if (s == 0)
        continue ;
      int c = T.Read(s - 1) ;
      Utils::BitSet(_psiB[c], s) ;
      ++alphabetCount[c] ;   
    }
    _lastChr = T.Read(n - 1) ;
    for (i = 1 ; i < alphabetSize ; ++i)
      alphabetCount[i] += alphabetCount[i - 1] ;
    for (i = 0 ; i < alphabetSize ; ++i)
      _alphabetPartialSum[i + 1] += alphabetCount[i] ;
  }

  // Compress everything
  void Init()
  {
    size_t i ;
    int alphabetSize = _alphabets.GetSize() ;
    
    for (i = 0 ; i < alphabetSize ; ++i)
      _Psi[i].Init(_psiB[i], _n) ;

    for (i = 0 ; i < alphabetSize ; ++i)
      free(_psiB[i]) ;   
    free(_psiB) ;
    _psiB = NULL ;
  }

  void InitFromBWT(FixedSizeElemArray &BWT, size_t n, size_t firstISA, ALPHABET *alphabetList)
  {
    size_t i ;
    _n = n ;
    _firstISA = firstISA ;

    Prepare(n, alphabetList) ;

    // Prepare the auxiliary data, e.g. rank, for BWT
    int alphabetSize = _alphabets.GetSize() ;
    _alphabetPartialSum = (size_t *)calloc(alphabetSize, sizeof(size_t)) ;
    
    for (i = 0 ; i < n ; ++i)
      ++_alphabetPartialSum[ BWT.Read(i)] ;
    for (i = 1 ; i < alphabetSize ; ++i)
      _alphabetPartialSum[i] += _alphabetPartialSum[i - 1] ;
    for (i = alphabetSize ; i > 0 ; --i)
      _alphabetPartialSum[i] = _alphabetPartialSum[i - 1] ;
    _alphabetPartialSum[0] = 0 ;
    _lastChr = alphabetList[ BWT.Read(firstISA) ] ;
    
    Sequence_WaveletTree<Bitvector_Plain> seqBWT ;
    seqBWT.Init(BWT, n, alphabetList) ;
  
    // Compute Psi from BWT. The Psi for each alphabet is marked on the bit array psiB
    size_t lastISA = _alphabetPartialSum[ _alphabets.Encode(_lastChr) ] ; 
    size_t p = lastISA;
    for (i = 0 ; i < n ; ++i)
    {
      int c = BWT.Read(p) ;
      size_t lf = _alphabetPartialSum[c] + Rank(seqBWT, alphabetList[c], p) ; 
      Utils::BitSet(_psiB[c], p) ; // Psi[lf] = p. Psi is the inverse function of LF mapping  
      p = lf ;
    }

    Init() ;
  }    
} ;
}

#endif 

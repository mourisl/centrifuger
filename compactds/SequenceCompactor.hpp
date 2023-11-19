#ifndef _MOURISL_COMPACTDS_SEQUENCECOMPACTOR
#define _MOURISL_COMPACTDS_SEQUENCECOMPACTOR

// The class that handles convert the raw sequence to FixedSizeElemArray
// I put this class in compactds because FM and Sequence classes assumes 
//  the input is from compact representation
#include "FixedSizeElemArray.hpp"
#include "Alphabet.hpp"

namespace compactds {
class SequenceCompactor
{
private:
  bool _capitalize ;
  ALPHABET _missingReplace ; 
  Alphabet _alphabets ;
public: 
  SequenceCompactor() 
  {
    _capitalize = false ;
    _missingReplace = '\0' ;  
  };

  ~SequenceCompactor() {} ;
  
  void Init(const char *alphabetList)
  {
    _alphabets.InitFromList(alphabetList, strlen(alphabetList)) ;
  }

  void Init(const char *alphabetList, FixedSizeElemArray &compactSeq, size_t reserveLength)
  {
    int alphabetCodeLen = _alphabets.InitFromList(alphabetList, strlen(alphabetList)) ;
    compactSeq.Malloc(alphabetCodeLen, reserveLength) ; 
    compactSeq.SetSize(0) ;
  }

  void SetCapitalize(bool c)
  {
    _capitalize = c ;
  }

  void SetMissingReplace(ALPHABET c)
  {
    _missingReplace = c ;
  }

  // @return: number of chars added to seq
  size_t Compact(const char *rawseq, FixedSizeElemArray &seq) 
  {
    size_t i ;
    size_t origLen = seq.GetSize() ;
    for (i = 0 ; rawseq[i] ; ++i)
    {
      char c = rawseq[i] ;
      if (_capitalize)
      {
        if (c >= 'a' && c <= 'z')
          c = c - 'a' + 'A' ;
      }

      if (!_alphabets.IsIn(c))
      {
        if (_missingReplace == '\0')
          continue ;
        else
          c = _missingReplace ;
      }
      seq.PushBack( _alphabets.Encode(c) ) ;
    }

    return seq.GetSize() - origLen ;
  }
} ;
} 
#endif

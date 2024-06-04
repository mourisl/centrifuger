#ifndef _MOURISL_COMPACTDS_DS_ALPHABET
#define _MOURISL_COMPACTDS_DS_ALPHABET

#include "Utils.hpp"
#include "HuffmanCode.hpp"
#include "FixedSizeElemArray.hpp"

typedef char ALPHABET ;

#define ALPHABET_CODE_NOCODE 0
#define ALPHABET_CODE_PLAIN 1
#define ALPHABET_CODE_HUFFMAN 2

// The data structe for mapping alphabet
// Conceptually, all the other data structure regard the alphabet as {0,...,|sigma|-1},
// This function serves to map these numeric alphabet to actually alphabet(char by default).
namespace compactds {
class Alphabet
{
private:
  size_t _space ;
  int _method ;  
  ALPHABET *_alphabetList ;    
  int _alphabetCode[1<<(sizeof(ALPHABET) * 8)] ;
  short _alphabetCodeLen[1<<(sizeof(ALPHABET) * 8)] ; // the length of encoded bits.
  size_t _n ;

  HuffmanCode _huffmanCode ; 
public:
  Alphabet() 
  {
    _n = _space = 0 ;
    _method = ALPHABET_CODE_NOCODE ;
  }

  ~Alphabet() { Free() ; }

  void Free()
  {
    if (_n != 0)
    {
      free(_alphabetList) ;
      _n = 0 ;
    }
  }
  
  size_t GetSpace() { return _space + sizeof(*this); } 

  // Use plain binary number sequentially for the characters in s.
  // @return: code length
  int InitFromList(const ALPHABET *s, size_t n)
  {
    size_t i ;
    this->_n = n ;
    _alphabetList = (ALPHABET *)malloc(sizeof(ALPHABET) * n) ;
    _space = sizeof(ALPHABET) * n ;  
    memset(_alphabetCode, 0, sizeof(_alphabetCode)) ;
    memset(_alphabetCodeLen, 0, sizeof(_alphabetCodeLen)) ;
    
    int codeLen = Utils::Log2Ceil(n) ;
    for (i = 0 ; i < n ; ++i)
    {
      _alphabetList[i] = s[i] ;
      _alphabetCode[ (int)s[i] ]= i ;
      _alphabetCodeLen[ (int)s[i] ] = codeLen ; 
    }
    _method = ALPHABET_CODE_PLAIN ;
    return codeLen ;
  }

  // s: list of the characters
  // freq: list of the frequencies for each character
  // n: number of character
  void InitHuffman(const ALPHABET *s, const uint64_t *freq, size_t n)
  {
    size_t i ;
    this->_n = n ;
    _alphabetList = (ALPHABET *)malloc(sizeof(ALPHABET) * n) ;
    for (i = 0 ; i < n ; ++i)
      _alphabetList[i] = s[i] ;
   
    _huffmanCode.InitFromFrequency(freq, n) ;

    for (i = 0 ; i < n ; ++i)
    {
      int l ;
      _alphabetCode[i] = _huffmanCode.Encode(i, l) ;
      _alphabetCodeLen[i] = l ;
    }
    _method = ALPHABET_CODE_HUFFMAN ;
  }

  size_t GetAlphabetCapacity() const
  {
    if (ALPHABET_CODE_NOCODE)
      return 0 ;
    else if (ALPHABET_CODE_PLAIN)
      return 1<<(Utils::Log2Ceil(_n)) ;
    else if (ALPHABET_CODE_HUFFMAN)
      return _n ;
    return 0 ;
  }

  size_t GetLongestCodeLength() const
  {
    if (ALPHABET_CODE_NOCODE)
      return 0 ;
    else
    {
      size_t i ;
      size_t ret = 0 ;
      for (i = 0 ; i < _n ; ++i)
      {
        size_t tmp = _alphabetCodeLen[(int)_alphabetList[i]] ;
        if (tmp > ret)
          ret = tmp ;
      }
      return ret ;
    }
    return 0 ;
  }

  size_t GetSize() const
  {
    return _n ;
  }
  
  // l: how many bits used in the coding
  ALPHABET Decode(WORD c, int l) const
  {
    //l = _alphabetCodeLen[ (int)_alphabetList[i] ] ;
    size_t i ;
    if (_method == ALPHABET_CODE_NOCODE)
    {
      return c ;
    }

    if (_method == ALPHABET_CODE_PLAIN)
      i = c ;
    else
      i = _huffmanCode.Decode(c, l) ;
    return _alphabetList[i] ;
  }

  WORD Encode(ALPHABET c, int &l) const
  {
    if (_method == ALPHABET_CODE_NOCODE)
    {
      //l = Utils::CountBits(c) ;
      l = 0 ;
      return c ;
    }
    else
    {
      l = _alphabetCodeLen[(int)c] ;
      return _alphabetCode[(int)c] ;
    }
  }
  
  WORD Encode(ALPHABET c) const
  {
    if (_method == ALPHABET_CODE_NOCODE)
      return c ;
    else
      return _alphabetCode[(int)c] ;
  }

  // test whether the alphabet c is in the list
  bool IsIn(ALPHABET c) const
  {
    size_t i ;
    for (i = 0 ; i < _n ; ++i)
      if (_alphabetList[i] == c)
        return true ;
    return false ;
  }

  Alphabet& operator=(const Alphabet &in)
  {
    Free() ;
    _n = in._n ;
    _space = in._space ;
    _method = in._method ;

    _alphabetList = (ALPHABET *)malloc(sizeof(ALPHABET) * _n) ;
    _space = sizeof(ALPHABET) * _n ;  
    memcpy(_alphabetList, in._alphabetList, sizeof(ALPHABET) * _n ) ;
    memcpy(_alphabetCode, in._alphabetCode, sizeof(_alphabetCode)) ;
    memcpy(_alphabetCodeLen, in._alphabetCodeLen, sizeof(_alphabetCodeLen)) ;
    _huffmanCode = in._huffmanCode ;
    return *this ;
  }

  void Save(FILE *fp)
  {
    SAVE_VAR(fp, _space) ;
    SAVE_VAR(fp, _method) ;
    SAVE_VAR(fp, _n) ;
    if (_n != 0)
    {
      fwrite(_alphabetList, sizeof(ALPHABET), _n, fp) ;    
      fwrite(_alphabetCode, sizeof(_alphabetCode[0]), 1<<(sizeof(ALPHABET) * 8), fp) ;    
      fwrite(_alphabetCodeLen, sizeof(_alphabetCodeLen[0]), 1<<(sizeof(ALPHABET) * 8), fp) ; 
    }
  }

  void Load(FILE *fp)
  {
    Free() ;
    LOAD_VAR(fp, _space) ;
    LOAD_VAR(fp, _method) ;
    LOAD_VAR(fp, _n) ;
    
    if (_n != 0) // no list for empty sequence
    {
      _alphabetList = (ALPHABET *)malloc(sizeof(ALPHABET) * _n) ;
      fread(_alphabetList, sizeof(ALPHABET), _n, fp) ;   
      fread(_alphabetCode, sizeof(_alphabetCode[0]), 1<<(sizeof(ALPHABET) * 8), fp) ;    
      fread(_alphabetCodeLen, sizeof(_alphabetCodeLen[0]), 1<<(sizeof(ALPHABET) * 8), fp) ;   
    }
  }
} ;
}
#endif 

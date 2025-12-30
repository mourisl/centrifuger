#ifndef _MOURISL_COMPACTDS_DS_RANK
#define _MOURISL_COMPACTDS_DS_RANK

#include "Utils.hpp"
#include "FixedSizeElemArray.hpp"

// The standalone data structe for rank query on a plain bitvector
// Time complexity: constant time 
// Extra space complexity (bits): n/b + n/w*log(bw)
namespace compactds {
class DS_Rank
{
private:
  uint64_t *_R ; // the partial sum of 1s for blocks of the bit vector (right exclusive)
  FixedSizeElemArray _subR ; // the partial sum within each block for constant access
  int _b ; // block size, with respective to the word.
  int _bshift ; // the number of bits-1 of b
  size_t _wordCnt ;
  size_t _space ;
public:
  DS_Rank() 
  {
    _R = NULL ;
    _b = _space = 0 ;
  }

  DS_Rank(int blockSize, const WORD *B, const int &n) 
  {
    Init(blockSize, B, n) ;
  }

  ~DS_Rank() { Free() ; }

  void Free()
  {
    if (_R != NULL)
    {
      free(_R) ;
      _R = NULL ;
    }
    _b = 0 ;
  }
  
  size_t GetSpace() { return _space + sizeof(*this); } 

  // blockSize is the number of WORDs for each R 
  void Init(int blockSize, const WORD *B, const size_t &n)
  {
    size_t i ;
    _b = blockSize ;
    if (_b <= 0)
      _b = WORDBITS ;

    i = _b >> 1 ;
    for (_bshift = 0 ; i != 0 ; i >>= 1, ++_bshift) 
      ;
      
    _wordCnt = Utils::BitsToWords(n) ;
    size_t blockCnt = DIV_CEIL(_wordCnt, _b) ;
    _R = (uint64_t *)malloc(sizeof(uint64_t) * blockCnt) ;
    _space += sizeof(uint64_t) * blockCnt ;
    _subR.Malloc(Utils::Log2Ceil((_b-1)*WORDBITS), _wordCnt - blockCnt) ; // we don't need to store the first sub-block in each block
    _space += _subR.GetSpace() - sizeof(_subR) ;
    uint64_t onecntSum = 0 ;
    size_t localOneCntSum = 0 ;
    for (i = 0 ; i < _wordCnt ; ++i)
    {
      if (i % _b == 0)
      {
        _R[i/_b] = onecntSum ;
        localOneCntSum = 0 ;
      }
      else
      {
        _subR.Write(i - i / _b - 1, localOneCntSum) ;
      }
      int onecnt = Utils::Popcount(B[i]) ;
      onecntSum += onecnt ;
      localOneCntSum += onecnt ;
    }
  }
  
  int GetBlockSize() const // unit in word
  {
    return _b ;
  }

  int GetSubBlockSize() const
  {
    return WORDBITS ;  
  }

  uint64_t *GetR() const 
  {
    return _R ;
  }

  const FixedSizeElemArray *GetSubR() const
  {
    return &_subR ;
  }


  size_t Query(size_t i, const WORD *B, const size_t &n, int inclusive = 1) const
  {
    if (i >= n)
      return Query(n - 1, B, n, inclusive) ;

    size_t wi = i >> WORDBITS_WIDTH ;
    return _R[wi >> _bshift] + ((wi&(_b - 1)) ? _subR.Read(wi - (wi >> _bshift) - 1) : 0) 
      + Utils::Popcount(B[wi] & ((MASK(i&(WORDBITS - 1))<<inclusive) + inclusive)) ;
      //+ Utils::Popcount(B[wi] & MASK_WCHECK((i&(WORDBITS - 1)) + inclusive)) ;
    /*else
    {
      // The implemenation without _subR
      size_t j ;
      uint64_t onecntSum = _R[wi / b] ;
      for (j = (wi / b) * _b ; (j + 1) * WORDBITS <= i ; ++j)
        onecntSum += Utils::Popcount(B[j]) ;     
      return onecntSum + Utils::Popcount(B[j] & MASK(i % WORDBITS + inclusive)) ;
    }*/
  }
  
  void Save(FILE *fp)
  {
    fwrite(this, sizeof(*this), 1, fp) ;
    size_t blockCnt = DIV_CEIL(_wordCnt, _b) ;
    fwrite(_R, sizeof(_R[0]), blockCnt, fp) ;
    _subR.Save(fp) ;
  }

  void Load(FILE *fp)
  {
    fread(this, sizeof(*this), 1, fp) ;
    size_t blockCnt = DIV_CEIL(_wordCnt, _b) ;
    if (_R != NULL)
      free(_R) ;
    _R = (uint64_t *)malloc(sizeof(uint64_t) * blockCnt * 2) ;
    fread(_R, sizeof(_R[0]), blockCnt, fp) ;
    _subR.Load(fp) ;
  }

} ;

// Optimized for recording offset every 8 words.
class DS_Rank9
{
private:
  uint64_t *_R ; // the partial sum of 1s for blocks of the bit vector (right exclusive), even position for large block, odd positins for compact subblock
  size_t _wordCnt ;
  size_t _space ;
public:
  DS_Rank9() 
  {
    _R = NULL ;
    _space = 0 ;
  }

  DS_Rank9(const WORD *B, const int &n) 
  {
    Init(B, n) ;
  }

  ~DS_Rank9() { Free() ; }

  void Free()
  {
    if (_R != NULL)
    {
      free(_R) ;
      _R = NULL ;
    }
  }
  
  size_t GetSpace() { return _space + sizeof(*this); }

  int GetBlockSize() const
  {
    return 8 ;
  }

  int GetSubBlockSize() const
  {
    return WORDBITS ;
  }
  
  size_t GetWordCnt() const
  {
    return _wordCnt ;
  }
 
  // wi is the word index
  size_t GetRValue(size_t wi) const
  {
    const size_t ri = (wi >> 3) * 2 ; // region/block id
    const size_t t = (wi & 7) - 1 ; // the offset in the subblock
    return _R[ri] + ((_R[ri + 1]>> ((t + ((t>>60)&8))*9)) & 0x1ff) ; 
  }

  uint64_t *GetR() const
  {
    return _R ;
  }

  // blockSize is the number of WORDs for each R 
  void Init(const WORD *B, const size_t &n)
  {
    size_t i ;
    const int b = 8 ; // number of word in each block 
    const int subrWidth = Utils::Log2Ceil((b-1) * 64) ; // should equal to 9 
    _wordCnt = Utils::BitsToWords(n) ;
    size_t blockCnt = DIV_CEIL(_wordCnt, b) ;
    _R = (uint64_t *)calloc(blockCnt * 2, sizeof(uint64_t)) ;
    _space = sizeof(uint64_t) * blockCnt * 2 ;
    uint64_t onecntSum = 0 ;
    size_t localOneCntSum = 0 ;
    for (i = 0 ; i < _wordCnt ; ++i)
    {
      size_t bi = i/b * 2 ; // block index
      int br = i % b ; //remainder
      if (br == 0)
      {
        _R[bi] = onecntSum ;
        _R[bi + 1] = 0 ;
        localOneCntSum = 0 ;
      }
      else
      {
        _R[bi + 1] |= (localOneCntSum << ((br - 1)*subrWidth))  ;
      }
      int onecnt = Utils::Popcount(B[i]) ;
      onecntSum += onecnt ;
      localOneCntSum += onecnt ;
    }
    // Fill in the remaining subr blocks
    //   so other module don't need to worry
    //   too much about boundary case
    if ((i-1) % b > 0)
    {
      size_t bi = i/b * 2 ; // block index
      for ( ; i % b ; ++i)
      {
        int br = i % b ;
        _R[bi + 1] |= (localOneCntSum << ((br - 1)*subrWidth))  ;
      }
    }
  }
  
  // read the si-th subblock in block bi
  int DecodeSubR(size_t bi, size_t si) const
  {
    return (_R[2 * bi + 1] >> (si * 9)) & 0x1ff;
  }
  
  size_t Query(size_t i, const WORD *B, const size_t &n, int inclusive = 1) const
  {
    //if (n == 0). Should be handle externally to make sure no empty array is accessed (for now)
    //  return 0 ;
    if (i >= n)
      return Query(n - 1, B, n, inclusive) ;
    
    const size_t wi = (i>>WORDBITS_WIDTH) ; // word id
    const size_t ri = (wi >> 3) * 2 ; // region/block id
    const size_t t = (wi & 7) - 1 ; // the offset in the subblock
    // 0x1ff is the mask for 9 bit
    // The ((t>>60)&8))*9) portion is to avoid branching when wi%8 == 0
    //  In this case, t=0xffff.., and (t + ((t>>60)&8))*9) == 63
    //  and the top bit of _R[ri+1] is 0, which makes the whole portion == 0
    // The implementation of inclusive also avoids branching using property that 
    //  inclusive variable is binary.
    return _R[ri] + ((_R[ri + 1]>> ((t + ((t>>60)&8))*9)) & 0x1ff) 
      + Utils::Popcount(B[wi] & ((MASK(i&(WORDBITS - 1))<<inclusive) + inclusive)) ;
  }

  void Save(FILE *fp)
  {
    SAVE_VAR(fp, _space) ;
    SAVE_VAR(fp, _wordCnt) ;
    const int b = 8 ;
    size_t blockCnt = DIV_CEIL(_wordCnt, b) ;
    fwrite(_R, sizeof(_R[0]), blockCnt * 2, fp) ; 
  }

  void Load(FILE *fp)
  {
    Free() ;

    LOAD_VAR(fp, _space) ;
    LOAD_VAR(fp, _wordCnt) ;
    const int b = 8 ;
    size_t blockCnt = DIV_CEIL(_wordCnt, b) ;
    if (_R != NULL)
      free(_R) ;
    _R = (uint64_t *)malloc(sizeof(uint64_t) * blockCnt * 2) ;
    fread(_R, sizeof(_R[0]), blockCnt * 2, fp) ; 
  }
} ;
}

#endif 

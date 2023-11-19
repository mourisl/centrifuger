#ifndef _MOURISL_COMPACTDS_DS_SELECT
#define _MOURISL_COMPACTDS_DS_SELECT

#include "Utils.hpp"
#include "DS_Rank.hpp"

// The standalone data structe for select query on a plain bitvector with precomputed rank information 
// _n - bitvector length. m - number of 1s (or 0s for select0)
// Speed 1: Time complexity: O(log n/m) [_space: O(n/w)]
// Speed 2: Time complexity: O(log log n) [_space: O(n/log n)] 
//          Most of the _space are reuse the rank structure though, 
//          so in practice the extra _space is stil only O(n/w)
// Speed 3: Time complexity: O(log log n) [_space: O(n/log n)]
//          inspired by the implementation in SDSL 
// Seems speed 4 does not work properly...
// Speed 4: Time complexity: O(1)  [_space: O(n loglog _n / sqrt(log n) + sqrt(n))]
//          The textbook O(n/log log n)-_space algorithm has too large factor
//          for precomputed short _miniblocks. Not very pratical.
//

#define DS_SELECT_SPEED_NO 0
#define DS_SELECT_SPEED_SAMPLED 1
#define DS_SELECT_SPEED_RANKBINARY 2
#define DS_SELECT_SPEED_DENSESAMPLE 3
#define DS_SELECT_SPEED_CONSTANT 4

namespace compactds {
class DS_Select
{
private:
  size_t *_S[2] ; // sampled position for 0's and 1's

  // Data structures for long blocks
  size_t _longBlockLength ;
  WORD *_V[2] ; // indicator whether a S block is long (1) or short
  DS_Rank9 _rankV[2] ;
  FixedSizeElemArray _I[2] ; // precomputed index within long block

  int _precomputeb ; // the precomputed offsets within a word of size b
  int _precomputebElem ; // how many 1s we should consider for such word.

  int _minib ; // mini block size (the number of 1's) 
  size_t _longMiniBlockLength ; // long mini block length, for _speed==3,4
  WORD *_Vmini[2] ; // indicator whether a S block is long mini or not
  size_t _VminiSize[2] ;
  DS_Rank9 _rankVmini[2] ;
  FixedSizeElemArray _Imini[2] ; // offset for the beginning of mini block 
  FixedSizeElemArray _Ilongmini[2] ; // offset for each element in long-mini block
    
  // Concatenated precomputed short mini block's S. We need concatenation, otherwise too 
  // much overhead in the FixedSizeElemArray structure. 
  // Even without FixedSizeElemArray, the pointers will take too much _space.
  FixedSizeElemArray _precomputedShortMiniBlock[2] ; 

  int _b ; // block size (the number of 1's in a block) or sampling rate
  size_t _n ;
  size_t _totalOneCnt ; 
   
  size_t _space ;
  
  int _speed ; // 0: do not allocate; 1: slow, 2: medium, 3: medium-fast 4: fastest, constant time

  // The select that handles both types (0 or 1).
  size_t GeneralQuery(size_t i, const DS_Rank &rank, const WORD *B, const size_t &n, int type) const
  {
    if (i < 1 || i > _n)
      return POSITIVE_INF ;
    if (_n == 1) // i must == 1 here
      return 1 ;

    size_t si = (i - 1) / _b ;

    size_t l, m, r ; // variables for binary search. They are coordinates on B 
    l = _S[type][si] ;
    r = _S[type][si + 1] ;
    if ((i - 1) % _b == 0) 
      return l ;
    if (_speed == 1 || r - l < _longBlockLength) // r-l is more efficient than V.access.
    {
      if (_speed == 3)
      {
        // Adjust l, r using _miniblocks
        size_t oldl = l ;
        l = oldl + _Imini[type].Read((i - 1) / _minib) ;
        if ((i - 1) % _minib + 1 < (unsigned int)(_b / _minib)) // only adjust r if it is not in the last _miniblock in current block
          r = oldl + _Imini[type].Read((i - 1) / _minib + 1) ;
      }
      if (_speed <= 3)
      {
        --r ;

        // Locate the R
        uint64_t *rankR = rank.GetR() ; // rankR is right open
        int rankBlockSize = rank.GetBlockSize() ; // this block size is with respect to WORD 
        size_t rl, rr ; 
        size_t tmp ;
        rl = l / (rankBlockSize * WORDBITS) ;
        rr = r / (rankBlockSize * WORDBITS) ;
        while (rl <= rr)
        {
          m = (rl + rr) / 2 ;
          tmp = rankR[m] ;
          if (type == 0)
            tmp = m * rankBlockSize * WORDBITS - tmp ; 
          
          if (tmp < i)
            rl = m + 1 ;
          else
            rr = m - 1 ; // rankR[0]==0 makes sure m>=1 in the process
        }

        // Locate the subR
        size_t remaining ;
        if (type == 1)
          remaining = i - rankR[rr] ;
        else
          remaining = i - (rr * (rankBlockSize * WORDBITS) - rankR[rr]) ;
        size_t subrl, subrr, fixedSubrl ;
        int rankSubBlockSize = rank.GetSubBlockSize() ; // This block size is with respecto to bits
        const FixedSizeElemArray &rankSubR = *(rank.GetSubR()) ;
       
        subrl = rr * (rankBlockSize * WORDBITS / rankSubBlockSize - 1) ; // the first subr has offset 0, so we don't store them.
        subrr = subrl + rankBlockSize * WORDBITS / rankSubBlockSize - 2 ;
        if (subrr >= rankSubR.GetSize())
          subrr = rankSubR.GetSize() - 1 ;
        bool inFirstSubBlock = false ;
        if (rankSubR.GetSize() == 0 || (type == 1 && (uint32_t)rankSubR.Read(subrl) >= remaining)
            || (type == 0 && rankSubBlockSize - (uint32_t)rankSubR.Read(subrl) >= remaining)
            || subrl >= rankSubR.GetSize()) // The case that the last block has only one subblock, which will not be allocated
          inFirstSubBlock = true ;
        
        fixedSubrl = subrl ;
        if ( !inFirstSubBlock )
        {
          while (subrl <= subrr)
          {
            m = (subrl + subrr) / 2 ;
            tmp = rankSubR.Read(m) ;
            if (type == 0)
              tmp = (m - fixedSubrl + 1) * rankSubBlockSize - tmp ; // plus 1 here to incorporate the first sub block
            if (tmp < remaining)
              subrl = m + 1 ;
            else
              subrr = m - 1 ; // the in firstsubblock test makes sure this part won't under-flow 
          }
          
          if (type == 1)
            remaining -= rankSubR.Read(subrr) ;
          else
            remaining -= ((subrr - fixedSubrl + 1) * rankSubBlockSize - rankSubR.Read(subrr)) ;
        }

        // Processing the last WORD
        size_t lastWi = 0 ; // index of the last word
        WORD lastW = 0 ;
        if (inFirstSubBlock)
          lastWi = rr * rankBlockSize ;
        else
          lastWi = rr + (subrr + 1) * (rankSubBlockSize / WORDBITS) ; // here the rr is to compensate for the first subblock missed in every sub block
        lastW = B[lastWi] ;
        size_t j ;
        
        int sum = 0 ;
        for (j = 0 ; j < WORDBITS ; j += _precomputeb)
        {
          WORD x = (lastW >> j) & MASK(_precomputeb) ;
          int tmp = Utils::Popcount(x) ;
          if (type == 0)
            tmp = _precomputeb - tmp ;
          if (sum + tmp >= (int)remaining)
          {
            return lastWi * WORDBITS + j + _precomputedShortMiniBlock[type].Read(x * _precomputebElem + remaining - sum - 1) ;
          }
          sum += tmp ;
        }
        return POSITIVE_INF ; // should not reach here.
      }
      else // _speed >= 4
      {
        size_t skippedMiniBlocksInLong = _rankV[type].Query(si, _V[type], n, 0) * (_b/_minib) ;
        size_t iMini = (i - 1) / _minib - skippedMiniBlocksInLong ; 
        if ( Utils::BitRead(_Vmini[type], iMini) )
        {
          // long mini block
          if ((i-1) % _minib == 0)
          {
            return l + _Imini[type].Read(iMini) ;
          }
          else
          {
            size_t iLongMini = _rankVmini[type].Query(iMini - skippedMiniBlocksInLong,
                  _Vmini[type], _VminiSize[type], 0) ;
            //printf("%d\n", iLongMini * (_minib - 1) + (i - 1)%_minib - 1) ;
            /*printf("b=%d _minib=%d. i=%d iMini=%d skippedMini=%d iLongMini=%d l=%d _Imini[iMini]=%d x=%d _Ilongmini[x]=%d. ret=%d\n", 
                b, _minib, i, iMini, skippedMiniBlocksInLong, 
                iLongMini, l, _Imini[type].Read(iMini), 
                iLongMini * (_minib - 1) + (i-1)%_minib - 1, _Ilongmini[type].Read(iLongMini * (_minib - 1) + (i-1)%_minib - 1), 
                l + _Ilongmini[type].Read(iLongMini * (_minib - 1) + (i-1) % _minib - 1)) ;*/ 
            return l + _Ilongmini[type].Read(iLongMini * (_minib - 1) + (i-1) % _minib - 1) ; 
          }
        }
        else
        {
          // short mini block
          size_t offset = l + _Imini[type].Read(iMini) ;
          WORD localw = Utils::BitsRead(B, offset, offset + _longMiniBlockLength - 2) ;
          return offset + _precomputedShortMiniBlock[type].Read(localw * _minib + (i-1)%_minib) ;
        }
      }
    }
    else
    {
      // long block with sparse 1's
      size_t iI = (_rankV[type].Query(si, _V[type], n) - 1) * (_b - 1); // block index in I
      //printf("long block %d %d %d %d %d. %d\n", i, b, si, Utils::BitRead(V[type], si), iI, _I[type].Read(iI + (i - 1)%b - 1)) ;
      return _I[type].Read(iI + (i - 1)%_b - 1) ;
    }
  }
  
  // The select that handles both types (0 or 1). Optimized for rank9
  size_t GeneralQuery(size_t i, const DS_Rank9 &rank, const WORD *B, const size_t &n, int type) const
  {
    if (i < 1 || i > _n)
      return POSITIVE_INF ;
    if (_n == 1) // i must == 1 here
      return 1 ;
    
    size_t si = (i - 1) / _b ;

    size_t l, m, r ; // variables for binary search. They are coordinates on B 
    l = _S[type][si] ;
    r = _S[type][si + 1] ;
    if ((i - 1) % _b == 0) 
      return l ;
    if (_speed == 1 || r - l < _longBlockLength) // r-l is more efficient than V.access.
    {
      if (_speed == 3)
      {
        // Adjust l, r using _miniblocks
        size_t oldl = l ;
        const int factor = _b / _minib ; // there are factor _miniblocks in each block
        l = oldl + _Imini[type].Read((i - 1) / _minib) ;
        if ( (int)((i - 1)/_minib) % factor + 1 < factor // only adjust r if it is not in the last _miniblock in current block
           && r < _n  // and the _miniblock is not the last _miniblock in the all the array 
            ) 
          r = oldl + _Imini[type].Read((i - 1) / _minib + 1) ;
      }
      if (_speed <= 3)
      {
        --r ;

        // Locate the R
        uint64_t *rankR = rank.GetR() ;
        const int rankBlockSize = rank.GetBlockSize() ; // this block size is with respect to WORD 
        size_t rl, rr ; 
        size_t tmp ;
        rl = l / (rankBlockSize * WORDBITS) ;
        rr = r / (rankBlockSize * WORDBITS) ;
        while (rl <= rr)
        {
          m = (rl + rr) / 2 ;
          tmp = rankR[m << 1] ;
          if (type == 0)
            tmp = m * rankBlockSize * WORDBITS - tmp ; 
          
          if (tmp < i)
            rl = m + 1 ;
          else
            rr = m - 1 ; // rankR[0]==0 makes sure m>=1 in the process
        }

        // Locate the subR
        size_t remaining ;
        if (type == 1)
          remaining = i - rankR[rr<<1] ;
        else
          remaining = i - (rr * (rankBlockSize * WORDBITS) - rankR[rr<<1]) ;
        if (remaining == 512) // Happens only when the block is all 1 and we are query the last element. Number 512 requires more than 9 bits to represent
        {
          return rr * rankBlockSize * WORDBITS + 511 ;
        }
        // Mark the lowest bit for every 9-bit block
        const uint64_t l9 = 0x40201008040201ull ;
        const uint64_t h9 = l9 << 8 ; // mark the highest bit for every 9-bit block
        const uint64_t expandRem = remaining * l9 ;
        uint64_t subrWord = rankR[rr * 2 + 1] ;
        if (type == 0) // need to take corresponding complement to get the accumulate counts for 0
        {
          //64 + ((64*2)<<(9*1)) + ((64*3)<<(9*2)) + ((64*4)<<(9*3)) + ((64*5)<<(9*4)) + ((64*6)<<(9*5)) + ((64*7)<<(9*6)) = 0x7030140803000040ull
          subrWord = 0x7030140803010040ull - subrWord ;
        }
        uint64_t bitblockComp = BITBLOCK_LT(subrWord, expandRem, h9) ;
        size_t subrr = (((bitblockComp >> 8) * l9) >> 54ull) & 7ull ;
        // Processing the last WORD
        size_t lastWi ; // index of the last word
        lastWi = rr * rankBlockSize + subrr ;
        if (lastWi >= rank.GetWordCnt())
        {
          lastWi = rank.GetWordCnt() - 1 ;
          subrr = lastWi - rr * rankBlockSize ;
        }
        WORD lastW = B[lastWi] ;
        if (subrr > 0)
        {
          remaining -= ((subrWord >> ((subrr-1) * 9)) & 0x1ff) ;
        }
        
        if (type == 0)
          lastW = ~lastW ;
        
        return lastWi * WORDBITS + Utils::SelectInWord(lastW, remaining) ;
        //return POSITIVE_INF ; // should not reach here.
      }
      else // _speed >= 4
      {
        size_t skippedMiniBlocksInLong = _rankV[type].Query(si, _V[type], n, 0) * (_b/_minib) ;
        size_t iMini = (i - 1) / _minib - skippedMiniBlocksInLong ; 
        if ( Utils::BitRead(_Vmini[type], iMini) )
        {
          // long mini block
          if ((i-1) % _minib == 0)
          {
            return l + _Imini[type].Read(iMini) ;
          }
          else
          {
            size_t iLongMini = _rankVmini[type].Query(iMini - skippedMiniBlocksInLong,
                  _Vmini[type], _VminiSize[type], 0) ;
            //printf("%d\n", iLongMini * (_minib - 1) + (i - 1)%_minib - 1) ;
            /*printf("b=%d _minib=%d. i=%d iMini=%d skippedMini=%d iLongMini=%d l=%d _Imini[iMini]=%d x=%d _Ilongmini[x]=%d. ret=%d\n", 
                b, _minib, i, iMini, skippedMiniBlocksInLong, 
                iLongMini, l, _Imini[type].Read(iMini), 
                iLongMini * (_minib - 1) + (i-1)%_minib - 1, _Ilongmini[type].Read(iLongMini * (_minib - 1) + (i-1)%_minib - 1), 
                l + _Ilongmini[type].Read(iLongMini * (_minib - 1) + (i-1) % _minib - 1)) ;*/ 
            return l + _Ilongmini[type].Read(iLongMini * (_minib - 1) + (i-1) % _minib - 1) ; 
          }
        }
        else
        {
          // short mini block
          size_t offset = l + _Imini[type].Read(iMini) ;
          WORD localw = Utils::BitsRead(B, offset, offset + _longMiniBlockLength - 2) ;
          return offset + _precomputedShortMiniBlock[type].Read(localw * _minib + (i-1)%_minib) ;
        }
      }
    }
    else
    {
      // long block with sparse 1's
      size_t iI = (_rankV[type].Query(si, _V[type], n) - 1) * (_b - 1); // block index in I
      //printf("long block %d %d %d %d %d. %d\n", i, b, si, Utils::BitRead(V[type], si), iI, _I[type].Read(iI + (i - 1)%b - 1)) ;
      return _I[type].Read(iI + (i - 1)%_b - 1) ;
    }
  }
public:
  DS_Select() 
  {
    _S[0] = _S[1] = NULL ;
    _V[0] = _V[1] = NULL ;
    _Vmini[0] = _Vmini[1] = NULL ;
    _n = _totalOneCnt = _b = _space = 0 ;
  }

  DS_Select(int blockSize, const WORD *B, const int &n, int selectSpeed, int selectTypeSupport) 
  {
    Init(blockSize, B, n, selectSpeed, selectTypeSupport) ;
  }

  ~DS_Select() { Free() ; }

  void Free()
  {
    int i ;
    for (i = 0 ; i <= 1 ; ++i)
    {
      if (_S[i] != NULL)
      {
        free(_S[i]) ;
        _S[i] = NULL ;
      }
      
      if (_V[i] != NULL)
      {
        free(_V[i]) ;
        _V[i] = NULL ; 
      }
      _rankV[i].Free() ;
      _I[i].Free() ;
    
      if (_Vmini[i] != NULL)
      {
        free(_Vmini[i]) ;
        _Vmini[i] = NULL ;
      }
      _rankVmini[i].Free() ;
      _Imini[i].Free() ;
      _Ilongmini[i].Free() ;
      _precomputedShortMiniBlock[i].Free() ;
    } 
    _n = _b = 0 ;
  }
  
  size_t GetSpace() { return _space + sizeof(*this); } 

  // blockSize is the number of WORDs for each R 
  // selectTypeSupport: bit coding for whether allocate memory to support select0 and select1
  //  0-bit: select 0, 1-bit: selct1; so 3 means support both
  void Init(int blockSize, const WORD *B, const size_t &n, int selectSpeed, int selectTypeSupport)
  {
    _speed = selectSpeed ;
    this->_n = n ;
    if (selectSpeed == 0 || selectTypeSupport == 0 || n <= 1)
      return ;
    size_t i, j ;
    size_t wordCnt = Utils::BitsToWordBytes(n) / sizeof(WORD) ;
    size_t *posBuffer = NULL;
    _space = 0 ;
    _b = blockSize ;

    // Set the parameters based the desired _speed
    if (_b <= (int)WORDBITS)
    {
      _b = WORDBITS * WORDBITS;
      if (_speed == 2)
        _b = WORDBITS * Utils::Log2Ceil(n) ; //* Utils::Log2Ceil( Utils::Log2Ceil(n) ) ;
      if (_speed == 3)
        _b = WORDBITS * WORDBITS ;
      if (_speed == 4)
        _b = WORDBITS * Utils::Log2Ceil(n) ;
    }
    
    int logn = Utils::Log2Ceil(n) ;
    //_longBlockLength = _b * Utils::Log2Ceil(n) * Utils::Log2Ceil(n) ; // Two sampled 1's are too far apart. It should be b*log^2 n
    _longBlockLength = logn * logn * logn * logn ; // Two sampled 1's are too far apart. It should be log^4 n
    if (_longBlockLength < (unsigned int)_b)
      _longBlockLength = _b ;

    _longMiniBlockLength = 0 ;
    if (_speed == 2 || _speed == 3)
    {
      if (n >= (1<<30))
        _precomputeb = 16 ;  // relate to precomputed select 
      else
        _precomputeb = 8 ;
      _precomputebElem = _precomputeb ;
      if (_speed == 3)
      {
        _minib = 2 * WORDBITS ;//logn * logn ;//CEIL(sqrt((double)b)) ;
        _minib -= _b % _minib ;  
        if (_minib < 3)
        {
          _minib = 3 ;
          if (_b % 3)
            _minib = 3 + _b%3 ;
        }
      }
    }
    else if (_speed == 4)
    {
      //_minib = sqrt(log n)
      _minib = CEIL(pow((double)_b, 0.25)) ; // We make _minib depends on the choice of _b so it is easier to control the block size.
      _minib -= _b % _minib ;  
      if (_minib < 3)
      {
        _minib = 3 ;
        if (_b % 3)
          _minib = 3 + _b%3 ;
      }
      _longMiniBlockLength = DIV_CEIL(_minib * _minib, 2) ;
      posBuffer = (size_t*)malloc(sizeof(*posBuffer) * (_b+1)) ;
      _precomputeb = _longMiniBlockLength - 1 ;
      _precomputebElem = _minib ;
    }

    _totalOneCnt = 0 ;
    for (i = 0 ; i < wordCnt ; ++i)
      _totalOneCnt += Utils::Popcount(B[i]) ;
    
    // Sample every other _b 1's (or 0's)  
    size_t blockCnt[2] ;
    blockCnt[0] = DIV_CEIL((n - _totalOneCnt), _b) + 1 ;
    blockCnt[1] = DIV_CEIL(_totalOneCnt, _b) + 1 ;
    for (i = 0 ; i <= 1 ; ++i)
    {
      if (!(selectTypeSupport & (1<<i)))
        continue ;
      _S[i] = (size_t *)calloc(blockCnt[i], sizeof(size_t)) ;
      _space += sizeof(size_t) * blockCnt[i] ;
      
      _S[i][blockCnt[i] - 1] = _n ; 
    }
    
    uint64_t sum[2] = {0, 0} ;
    // _S[bit][x] mark a block of "bit" starting at index x
    for (i = 0 ; i < wordCnt ; ++i)
    {
      WORD tmp = B[i] ;
      for (j = 0 ; j < WORDBITS && (i * WORDBITS + j < n); ++j)
      {
        int bit = (tmp >>j) & 1ull ;
        if (!(selectTypeSupport & (1<<bit)))
          continue ;
        if (sum[bit] % _b == 0)
          _S[bit][sum[bit] / _b] = i * WORDBITS + j ;
        ++sum[bit] ;
      }
    }
    
    // Sample more index for faster query
    if (_speed >= 2)
    {
      int k ;
      for (k = 0 ; k <= 1 ; ++k)
      {
        if (!(selectTypeSupport & (1<<k)))
          continue ;
        _V[k] = Utils::MallocByBits(blockCnt[k]) ;
        _space += Utils::BitsToWords(blockCnt[k]) ;
        _I[k].Malloc(Utils::Log2Ceil(n), DIV_CEIL(n, _longBlockLength)*_b) ;
        
        if (_speed >= 3)
          _Imini[k].Malloc(Utils::Log2Ceil(_longBlockLength), blockCnt[k] * (_b/_minib)) ;
        
        if (_speed >= 4)
        {
          _Vmini[k] = Utils::MallocByBits(blockCnt[k] * (_b / _minib)) ;
          // The long mini block can be almost as larage as long block length - 1
          _Ilongmini[k].Malloc(Utils::Log2Ceil(_longBlockLength), DIV_CEIL(n, _longMiniBlockLength) * _minib) ;
        }
        
        size_t newISize = 0 ;
        size_t new_IminiSize = 0 ;
        size_t new_IlongminiSize = 0 ;
        for (i = 0 ; i < blockCnt[k] - 1 ; ++i)
        {
          if (_S[k][i + 1] - _S[k][i] >= _longBlockLength) 
          {
            Utils::BitSet(_V[k], i) ;
            // The first element is already stored in S, so no need to store it
            for (j = _S[k][i] + 1 ; j < _S[k][i + 1] ; ++j)
            {
              if (Utils::BitRead(B, j) == k)
              {
                _I[k].Write(newISize, j) ;
                ++newISize ;
              }
            }
            
            if (_speed == 3) // For _speed 3, we still fill up I mini 
                            // so we don't need to acces _rankV for efficiency.
                            // Maybe I should do this to _speed 4 as well.
            {
              for (j = 0 ; j < (size_t)(_b / _minib) ; ++j)
              {
                _Imini[k].Write(new_IminiSize, 0) ;
                ++new_IminiSize ;
              }
            }
          }
          else if (_speed >= 3) // short block case, we only need to process them when _speed==3
          {
            int minicnt = 1; 
            size_t prevj = _S[k][i] ;
            if (_speed >= 4)
              posBuffer[0] = _S[k][i] ;
            // j reaches the beginning of the next block so we can wrap up any unadded 
            //   k's to the _miniblock. This handles both case that the last _miniblock in 
            //   a block or the last _miniblock in the whole bit vector.
            for (j = _S[k][i] + 1; j <= _S[k][i + 1] ; ++j)
            {
              int bit = 0 ;
              if (j < n)
                bit = Utils::BitRead(B, j) ;
              if (bit == k || (j == _S[k][i + 1] && minicnt > 0))
              {
                if (minicnt == _minib || (j == _S[k][i + 1] && minicnt > 0))
                {
                  _Imini[k].Write(new_IminiSize, prevj - _S[k][i]) ;
                  ++new_IminiSize ;

                  if (_speed >= 4 && j - prevj >= _longMiniBlockLength)
                  {
                    int l ;
                    Utils::BitSet(_Vmini[k], new_IminiSize - 1) ;
                    for (l = 1 ; l < minicnt ; ++l) // we don't need to store the first element
                    {
                      _Ilongmini[k].Write(new_IlongminiSize, posBuffer[l] - _S[k][i]) ;
                      ++new_IlongminiSize ;
                    }
                  }
                  
                  prevj = j ;
                  minicnt = 0 ;
                }
                
                if (bit == k)
                {
                  if (_speed >= 4)
                    posBuffer[minicnt] = j ;
                  ++minicnt ;
                }
              }
            }
          }
        }
        _I[k].Resize(newISize) ;
        _space += _I[k].GetSpace() - sizeof(_I[k]) ;
        
        _rankV[k].Init(_V[k], blockCnt[k]) ;
        _space += _rankV[k].GetSpace() - sizeof(_rankV[k]) ;
        
        if (_speed >= 3)
        {
          _Imini[k].Resize(new_IminiSize) ;
          _space += _Imini[k].GetSpace() - sizeof(_Imini[k]) ;
          //printf("%d %d. %d. %d %d\n", _Imini[k].GetSpace(), new_IminiSize, Utils::Log2Ceil(_longBlockLength), _minib, n/_minib) ;  
        }

        if (_speed >= 4)
        {
          _Vmini[k] = (WORD *)realloc(_Vmini[k], 
              Utils::BitsToWordBytes(new_IminiSize)) ;
          _VminiSize[k] = new_IminiSize ;
          _space += Utils::BitsToWordBytes(new_IminiSize) ;
          _rankVmini[k].Init(_Vmini[k], new_IminiSize) ;
          _space += _rankVmini[k].GetSpace() - sizeof(_rankVmini[k]) ;

          _Ilongmini[k].Resize(new_IlongminiSize) ;
          _space += _Ilongmini[k].GetSpace() - sizeof(_Ilongmini[k]) ;
        }
      }
    }
    
    if (0 && _speed >= 2) // Now we are using Rank9 and bit operator, so no need for the precomputed element
    {
      // The precomputed short _miniblocks 
      unsigned int k ;
      for (k = 0 ; k <= 1 ; ++k) 
      {
        if (!(selectTypeSupport & (1<<k)))
          continue ;
        size_t size = 1<<_precomputeb ;
        _precomputedShortMiniBlock[k].Malloc(Utils::Log2Ceil(_precomputeb), size * _precomputebElem) ;
        for (i = 0 ; i < size ; ++i) 
        {
          // Only consider the first _minib 1s 
          j = 0 ;
          size_t l ;
          for (l = 0 ; l < (unsigned int)_precomputeb ; ++l)
          {
            if (((i >> l) & 1ull)==k)
            {
              _precomputedShortMiniBlock[k].Write(i * _precomputebElem + j, l) ;
              ++j ;
              if ((int)j >= _precomputebElem)
                break ;
            }
          }
        }
        _space += _precomputedShortMiniBlock[k].GetSpace() - sizeof(_precomputedShortMiniBlock[k]) ;
      }
      if (_speed >= 4)
        free(posBuffer) ;
    }
  }

  // Return the index of the ith (1-index ith) 1.
  size_t Query(size_t i, const DS_Rank9 &rank, const WORD *B, const size_t &n) const
  {
    return GeneralQuery(i, rank, B, n, 1) ; 
  }

  // Return the index of the ith (1-index ith) 0.
  size_t Query0(size_t i, const DS_Rank9 &rank, const WORD *B, const size_t &n) const
  {
    return GeneralQuery(i, rank, B, n, 0) ;
  }

  void Save(FILE *fp)
  {
    SAVE_VAR(fp, _space) ; 
    SAVE_VAR(fp, _n) ;
    SAVE_VAR(fp, _speed) ;
    
    if (_speed == DS_SELECT_SPEED_NO || _n == 0)
      return ;
    
    SAVE_VAR(fp, _longBlockLength) ; 
    SAVE_VAR(fp, _minib);
    SAVE_VAR(fp, _longMiniBlockLength) ;
    SAVE_VAR(fp, _b) ;
    SAVE_VAR(fp, _totalOneCnt) ;

    size_t blockCnt[2] ;
    blockCnt[0] = DIV_CEIL((_n - _totalOneCnt), _b) + 1 ;
    blockCnt[1] = DIV_CEIL(_totalOneCnt, _b) + 1 ;
    for (int i = 0 ; i <= 1 ; ++i)
    {
      size_t size = Utils::BitsToWords(blockCnt[i]) ;
      fwrite(_S[i], sizeof(_S[i][0]), blockCnt[i], fp) ;
      if (_speed >= 2)
      {
        fwrite(_V[i], sizeof(_V[i][0]), size, fp) ;
        _rankV[i].Save(fp) ;
        _I[i].Save(fp) ;
      }
      if (_speed >= 3)
        _Imini[i].Save(fp) ;
    }
  }

  void Load(FILE *fp)
  {
    Free() ;

    LOAD_VAR(fp, _space) ; 
    LOAD_VAR(fp, _n) ;
    LOAD_VAR(fp, _speed) ;
    
    if (_speed == DS_SELECT_SPEED_NO || _n == 0)
      return ;
    
    LOAD_VAR(fp, _longBlockLength) ; 
    LOAD_VAR(fp, _minib);
    LOAD_VAR(fp, _longMiniBlockLength) ;
    LOAD_VAR(fp, _b) ;
    LOAD_VAR(fp, _totalOneCnt) ;

    size_t blockCnt[2] ;
    blockCnt[0] = DIV_CEIL((_n - _totalOneCnt), _b) + 1 ;
    blockCnt[1] = DIV_CEIL(_totalOneCnt, _b) + 1 ;
    for (int i = 0 ; i <= 1 ; ++i)
    {
      size_t size = Utils::BitsToWords(blockCnt[i]) ;
      _S[i] = (size_t *)malloc(sizeof(_S[i][0]) * blockCnt[i]) ;
      fread(_S[i], sizeof(_S[i][0]), blockCnt[i], fp) ;
      
      if (_speed >= 2)
      {
        _V[i] = Utils::MallocByBits(blockCnt[i]) ;
        fread(_V[i], sizeof(_V[i][0]), size, fp) ;
        _rankV[i].Load(fp) ;
        _I[i].Load(fp) ;
      }
      
      if (_speed >= 3)
        _Imini[i].Load(fp) ;
    }
  }
} ;
}

#endif 

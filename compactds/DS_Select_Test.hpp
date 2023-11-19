#ifndef _MOURISL_COMPACTDS_DS_SELECT_TEST
#define _MOURISL_COMPACTDS_DS_SELECT_TEST

#include "Utils.hpp"
#include "DS_Rank.hpp"

// The standalone data structe for select query on a plain bitvector with precomputed rank information 
// n - bitvector length. m - number of 1s (or 0s for select0)
// Speed 1: Time complexity: O(log n/m) [space: O(n/w)]
// Speed 2: Time complexity: O(log log n) [space: O(n/log n)] 
//          Most of the space are reuse the rank structure though, 
//          so in practice the extra space is stil only O(n/w)
// Seems speed 3, 4 does not work properly...
// Speed 3: Time complexity: O(log log n) [space: O(n/log n)]
//          inspired by the implementation in SDSL 
// Speed 4: Time complexity: O(1)  [space: O(n loglog n / sqrt(log n) + sqrt(n))]
//          The textbook O(n/log log n)-space algorithm has too large factor
//          for precomputed short miniblocks. Not very pratical.
//

#define DS_SELECT_SPEED_NO 0
#define DS_SELECT_SPEED_SAMPLED 1
#define DS_SELECT_SPEED_RANKBINARY 2
#define DS_SELECT_SPEED_DENSESAMPLE 3
#define DS_SELECT_SPEED_CONSTANT 4

namespace compactds {
class DS_Select_Test
{
private:
  size_t *S[2] ; // sampled position for 0's and 1's

  // Data structures for long blocks
  size_t longBlockLength ;
  WORD *V[2] ; // indicator whether a S block is long (1) or short
  DS_Rank rankV[2] ;
  FixedSizeElemArray I[2] ; // precomputed index within long block

  int precomputeb ; // the precomputed offsets within a word of size b
  int precomputebElem ; // how many 1s we should consider for such word.

  int minib ; // mini block size (the number of 1's) 
  size_t longMiniBlockLength ; // long mini block length, for speed==3,4
  WORD *Vmini[2] ; // indicator whether a S block is long mini or not
  size_t VminiSize[2] ;
  DS_Rank rankVmini[2] ;
  FixedSizeElemArray Imini[2] ; // offset for the beginning of mini block 
  FixedSizeElemArray Ilongmini[2] ; // offset for each element in long-mini block
    
  // Concatenated precomputed short mini block's S. We need concatenation, otherwise too 
  // much overhead in the FixedSizeElemArray structure. 
  // Even without FixedSizeElemArray, the pointers will take too much space.
  FixedSizeElemArray precomputedShortMiniBlock[2] ; 

  int b ; // block size (the number of 1's in a block) or sampling rate
  size_t n ;
  size_t totalOneCnt ; 
   
  size_t space ;
  
  int speed ; // 0: do not allocate; 1: slow, 2: medium, 3: medium-fast 4: fastest, constant time

  // The select that handles both types (0 or 1).
  size_t GeneralQuery(size_t i, const DS_Rank9 &rank, const WORD *B, const size_t &n, int type) const
  {
    if (i < 1 || i > n)
      return POSITIVE_INF ;
    if (n == 1) // i must == 1 here
      return 1 ;

    size_t si = (i - 1) / b ;

    size_t l, m, r ; // variables for binary search. They are coordinates on B 
    l = S[type][si] ;
    r = S[type][si + 1] ;
    if ((i - 1) % b == 0) 
      return l ;

    if (speed == 1 || r - l < longBlockLength) // r-l is more efficient than V.access.
    {
      if (speed == 3)
      {
        // Adjust l, r using miniblocks
        size_t oldl = l ;
        l = oldl + Imini[type].Read((i - 1) / minib) ;
        if ((i - 1) / minib + 1 < (unsigned int)(b / minib)) // only adjust r if it is not in the last miniblock in current block
          r = oldl + Imini[type].Read((i - 1) / minib + 1) ;
      }
      if (speed <= 3)
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
        
        int subrl, subrr ;
        subrl = 0 ; // the first subr has offset 0, so we don't store them.
        subrr = 6 ;
        if (rr * rankBlockSize + 1 + subrr >= rank.GetWordCnt())
          subrr = rank.GetWordCnt() - 2 - rr * rankBlockSize;
        bool inFirstSubBlock = false ;
        if (rank.GetWordCnt() <= 1 
            || (type == 1 && rank.DecodeSubR(rr, 0) >= remaining)
            || (type == 0 && WORDBITS - rank.DecodeSubR(rr, 0) >= remaining)
            || subrr < subrl) // The case that the last block has only one subblock, which will not be allocated
          inFirstSubBlock = true ;
        
        if ( !inFirstSubBlock )
        {
          size_t rword = rankR[2 * rr + 1] ;
          while (subrl <= subrr)
          {
            m = (subrl + subrr) / 2 ;
            tmp = (rword >> (m * 9ull)) & 0x1ff ;
            //printf("%d %d. %llu\n", m, tmp, rword) ;
            if (type == 0)
              tmp = (m + 1) * WORDBITS - tmp ; // plus 1 here to incorporate the first sub block
            if (tmp < remaining)
              subrl = m + 1 ;
            else
              subrr = m - 1 ; // the in firstsubblock test makes sure this part won't under-flow 
          }
          
          if (type == 1)
            remaining -= (rword >> (subrr * 9)) & 0x1ff ;
          else
            remaining -= ((subrr + 1) * WORDBITS - ((rword >> (subrr * 9)) & 0x1ff)) ;
        }
        
        // Processing the last WORD
        size_t lastWi = 0 ; // index of the last word
        WORD lastW = 0 ;
        if (inFirstSubBlock)
          lastWi = rr * rankBlockSize ;
        else
          lastWi = rr * rankBlockSize + subrr + 1 ; // here the rr is to compensate for the first subblock missed in every sub block
        lastW = B[lastWi] ;
        size_t j ;
          
        int sum = 0 ;
        for (j = 0 ; j < WORDBITS ; j += precomputeb)
        {
          WORD x = (lastW >> j) & MASK(precomputeb) ;
          int tmp = Utils::Popcount(x) ;
          if (type == 0)
            tmp = precomputeb - tmp ;
          if (sum + tmp >= (int)remaining)
          {
            return lastWi * WORDBITS + j + precomputedShortMiniBlock[type].Read(x * precomputebElem + remaining - sum - 1) ;
          }
          sum += tmp ;
        }
        return POSITIVE_INF ; // should not reach here.
      }
      else // speed >= 4
      {
        size_t skippedMiniBlocksInLong = rankV[type].Query(si, V[type], n, 0) * (b/minib) ;
        size_t iMini = (i - 1) / minib - skippedMiniBlocksInLong ; 
        if ( Utils::BitRead(Vmini[type], iMini) )
        {
          // long mini block
          if ((i-1) % minib == 0)
          {
            return l + Imini[type].Read(iMini) ;
          }
          else
          {
            size_t iLongMini = rankVmini[type].Query(iMini - skippedMiniBlocksInLong,
                  Vmini[type], VminiSize[type], 0) ;
            //printf("%d\n", iLongMini * (minib - 1) + (i - 1)%minib - 1) ;
            /*printf("b=%d minib=%d. i=%d iMini=%d skippedMini=%d iLongMini=%d l=%d Imini[iMini]=%d x=%d Ilongmini[x]=%d. ret=%d\n", 
                b, minib, i, iMini, skippedMiniBlocksInLong, 
                iLongMini, l, Imini[type].Read(iMini), 
                iLongMini * (minib - 1) + (i-1)%minib - 1, Ilongmini[type].Read(iLongMini * (minib - 1) + (i-1)%minib - 1), 
                l + Ilongmini[type].Read(iLongMini * (minib - 1) + (i-1) % minib - 1)) ;*/ 
            return l + Ilongmini[type].Read(iLongMini * (minib - 1) + (i-1) % minib - 1) ; 
          }
        }
        else
        {
          // short mini block
          size_t offset = l + Imini[type].Read(iMini) ;
          WORD localw = Utils::BitsRead(B, offset, offset + longMiniBlockLength - 2) ;
          return offset + precomputedShortMiniBlock[type].Read(localw * minib + (i-1)%minib) ;
        }
      }
    }
    else
    {
      // long block with sparse 1's
      size_t iI = (rankV[type].Query(si, V[type], n) - 1) * (b - 1); // block index in I
      //printf("long block %d %d %d %d %d. %d\n", i, b, si, Utils::BitRead(V[type], si), iI, I[type].Read(iI + (i - 1)%b - 1)) ;
      return I[type].Read(iI + (i - 1)%b - 1) ;
    }
  }
public:
  DS_Select_Test() 
  {
    S[0] = S[1] = NULL ;
    V[0] = V[1] = NULL ;
    Vmini[0] = Vmini[1] = NULL ;
    n = totalOneCnt = b = space = 0 ;
  }

  DS_Select_Test(int blockSize, const WORD *B, const int &n, int selectSpeed, int selectTypeSupport) 
  {
    Init(blockSize, B, n, selectSpeed, selectTypeSupport) ;
  }

  ~DS_Select_Test() { Free() ; }

  void Free()
  {
    int i ;
    for (i = 0 ; i <= 1 ; ++i)
    {
      if (S[i] != NULL)
      {
        free(S[i]) ;
        S[i] = NULL ;
      }
      
      if (V[i] != NULL)
      {
        free(V[i]) ;
        V[i] = NULL ; 
      }
      rankV[i].Free() ;
      I[i].Free() ;
    
      if (Vmini[i] != NULL)
      {
        free(Vmini[i]) ;
        Vmini[i] = NULL ;
      }
      rankVmini[i].Free() ;
      Imini[i].Free() ;
      Ilongmini[i].Free() ;
      precomputedShortMiniBlock[i].Free() ;
    } 
    n = b = 0 ;
  }
  
  size_t GetSpace() { return space + sizeof(*this); } 

  // blockSize is the number of WORDs for each R 
  // selectTypeSupport: bit coding for whether allocate memory to support select0 and select1
  //  0-bit: select 0, 1-bit: selct1; so 3 means support both
  void Init(int blockSize, const WORD *B, const size_t &n, int selectSpeed, int selectTypeSupport)
  {
    if (selectSpeed == 0 || selectTypeSupport == 0 || n <= 1)
      return ;
    size_t i, j ;
    size_t wordCnt = Utils::BitsToWordBytes(n) / sizeof(WORD) ;
    size_t *posBuffer = NULL;
    this->n = n ;
    speed = selectSpeed ;
    space = 0 ;
    b = blockSize ;

    // Set the parameters based the desired speed
    if (b <= (int)WORDBITS)
    {
      b = WORDBITS * WORDBITS;
      if (speed >= 2)
        b = WORDBITS * Utils::Log2Ceil(n) ; //* Utils::Log2Ceil( Utils::Log2Ceil(n) ) ;
      if (speed == 4)
        b = WORDBITS * Utils::Log2Ceil(n) ;
    }
    
    longBlockLength = b * Utils::Log2Ceil(n) * Utils::Log2Ceil(n) ; // Two sampled 1's are too far apart. It should be b*log^2 n
    if (speed == 2 || speed == 3)
    {
      if (n >= (1<<30))
        precomputeb = 16 ;  // relate to precomputed select 
      else
        precomputeb = 8 ;
      precomputebElem = precomputeb ;
      if (speed == 3)
      {
        minib = CEIL(sqrt((double)b)) ;
        minib -= b % minib ;  
        if (minib < 3)
        {
          minib = 3 ;
          if (b % 3)
            minib = 3 + b%3 ;
        }
      }
    }
    else if (speed == 4)
    {
      //minib = sqrt(log n)
      minib = CEIL(pow((double)b, 0.25)) ; // We make minib depends on the choice of b so it is easier to control the block size.
      minib -= b % minib ;  
      if (minib < 3)
      {
        minib = 3 ;
        if (b % 3)
          minib = 3 + b%3 ;
      }
      longMiniBlockLength = DIV_CEIL(minib * minib, 2) ;
      posBuffer = (size_t*)malloc(sizeof(*posBuffer) * (b+1)) ;
      precomputeb = longMiniBlockLength - 1 ;
      precomputebElem = minib ;
    }

    totalOneCnt = 0 ;
    for (i = 0 ; i < wordCnt ; ++i)
      totalOneCnt += Utils::Popcount(B[i]) ;
    
    // Sample every other b 1's (or 0's)  
    size_t blockCnt[2] ;
    blockCnt[0] = DIV_CEIL((n - totalOneCnt), b) + 1 ;
    blockCnt[1] = DIV_CEIL(totalOneCnt, b) + 1 ;
    for (i = 0 ; i <= 1 ; ++i)
    {
      if (!(selectTypeSupport & (1<<i)))
        continue ;
      S[i] = (size_t *)malloc(sizeof(size_t) * blockCnt[i]) ;
      space += sizeof(size_t) * blockCnt[i] ;
      
      S[i][blockCnt[i] - 1] = n ; 
    }
    
    uint64_t sum[2] = {0, 0} ;
    // S[bit][x] mark a block of "bit" starting at index x
    for (i = 0 ; i < wordCnt ; ++i)
    {
      WORD tmp = B[i] ;
      for (j = 0 ; j < WORDBITS && (i * WORDBITS + j < n); ++j)
      {
        int bit = (tmp >>j) & 1ull ;
        if (!(selectTypeSupport & (1<<bit)))
          continue ;
        if (sum[bit] % b == 0)
          S[bit][sum[bit] / b] = i * WORDBITS + j ;
        ++sum[bit] ;
      }
    }
    
    // Sample more index for faster query
    if (speed >= 2)
    {
      int k ;
      for (k = 0 ; k <= 1 ; ++k)
      {
        if (!(selectTypeSupport & (1<<k)))
          continue ;
        V[k] = Utils::MallocByBits(blockCnt[k]) ;
        space += Utils::BitsToWords(blockCnt[k]) ;
        I[k].Malloc(Utils::Log2Ceil(n), DIV_CEIL(n, longBlockLength)*b) ;
        
        if (speed >= 3)
          Imini[k].Malloc(Utils::Log2Ceil(longBlockLength), blockCnt[k] * (b/minib)) ;

        if (speed >= 4)
        {
          Vmini[k] = Utils::MallocByBits(blockCnt[k] * (b / minib)) ;
          // The long mini block can be almost as larage as long block length - 1
          Ilongmini[k].Malloc(Utils::Log2Ceil(longBlockLength), DIV_CEIL(n, longMiniBlockLength) * minib) ;
        }
        
        size_t newISize = 0 ;
        size_t newIminiSize = 0 ;
        size_t newIlongminiSize = 0 ;
        for (i = 0 ; i < blockCnt[k] - 1 ; ++i)
        {
          if (S[k][i + 1] - S[k][i] >= longBlockLength) 
          {
            Utils::BitSet(V[k], i) ;
            // The first element is already stored in S, so no need to store it
            for (j = S[k][i] + 1 ; j < S[k][i + 1] ; ++j)
            {
              if (Utils::BitRead(B, j) == k)
              {
                I[k].Write(newISize, j) ;
                ++newISize ;
              }
            }
            
            if (speed == 3) // For speed 3, we still fill up I mini 
                            // so we don't need to acces rankV for efficiency.
                            // Maybe I should do this to speed 4 as well.
            {
              for (i = 0 ; i < (size_t)(b / minib) ; ++i)
              {
                Imini[k].Write(newIminiSize, 0) ;
                ++newIminiSize ;
              }
            }
          }
          else if (speed >= 3) // short block case, we only need to process them when speed==3
          {
            int minicnt = 1; 
            size_t prevj = S[k][i] ;
            if (speed >= 4)
              posBuffer[0] = S[k][i] ;
            // j reaches the beginning of the next block so we can wrap up any unadded 
            //   k's to the miniblock. This handles both case that the last miniblock in 
            //   a block or the last miniblock in the whole bit vector.
            for (j = S[k][i] + 1; j <= S[k][i + 1] ; ++j)
            {
              int bit = 0 ;
              if (j < n)
                bit = Utils::BitRead(B, j) ;
              if (bit == k || (j == S[k][i + 1] && minicnt > 0))
              {
                if (minicnt == minib || (j == S[k][i + 1] && minicnt > 0))
                {
                  Imini[k].Write(newIminiSize, prevj - S[k][i]) ;
                  ++newIminiSize ;

                  if (speed >= 4 && j - prevj >= longMiniBlockLength)
                  {
                    int l ;
                    Utils::BitSet(Vmini[k], newIminiSize - 1) ;
                    for (l = 1 ; l < minicnt ; ++l) // we don't need to store the first element
                    {
                      Ilongmini[k].Write(newIlongminiSize, posBuffer[l] - S[k][i]) ;
                      ++newIlongminiSize ;
                    }
                  }
                  
                  prevj = j ;
                  minicnt = 0 ;
                }
                
                if (bit == k)
                {
                  if (speed >= 4)
                    posBuffer[minicnt] = j ;
                  ++minicnt ;
                }
              }
            }
          }
        }
        I[k].Resize(newISize) ;
        space += I[k].GetSpace() - sizeof(I[k]) ;
        
        rankV[k].Init(-1, V[k], blockCnt[k]) ;
        space += rankV[k].GetSpace() - sizeof(rankV[k]) ;
        
        if (speed >= 3)
        {
          Imini[k].Resize(newIminiSize) ;
          space += Imini[k].GetSpace() - sizeof(Imini[k]) ;
          //printf("%d %d. %d. %d %d\n", Imini[k].GetSpace(), newIminiSize, Utils::Log2Ceil(longBlockLength), minib, n/minib) ;  
        }

        if (speed >= 4)
        {
          Vmini[k] = (WORD *)realloc(Vmini[k], 
              Utils::BitsToWordBytes(newIminiSize)) ;
          VminiSize[k] = newIminiSize ;
          space += Utils::BitsToWordBytes(newIminiSize) ;
          rankVmini[k].Init(-1, Vmini[k], newIminiSize) ;
          space += rankVmini[k].GetSpace() - sizeof(rankVmini[k]) ;

          Ilongmini[k].Resize(newIlongminiSize) ;
          space += Ilongmini[k].GetSpace() - sizeof(Ilongmini[k]) ;
        }
      }
    }
    
    if (speed >= 2)
    {
      // The precomputed short miniblocks 
      unsigned int k ;
      for (k = 0 ; k <= 1 ; ++k) 
      {
        if (!(selectTypeSupport & (1<<k)))
          continue ;
        size_t size = 1<<precomputeb ;
        precomputedShortMiniBlock[k].Malloc(Utils::Log2Ceil(precomputeb), size * precomputebElem) ;
        for (i = 0 ; i < size ; ++i) 
        {
          // Only consider the first minib 1s 
          j = 0 ;
          size_t l ;
          for (l = 0 ; l < (unsigned int)precomputeb ; ++l)
          {
            if (((i >> l) & 1ull)==k)
            {
              precomputedShortMiniBlock[k].Write(i * precomputebElem + j, l) ;
              ++j ;
              if ((int)j >= precomputebElem)
                break ;
            }
          }
        }
        space += precomputedShortMiniBlock[k].GetSpace() - sizeof(precomputedShortMiniBlock[k]) ;
      }
      if (speed >= 4)
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
} ;
}
#endif 

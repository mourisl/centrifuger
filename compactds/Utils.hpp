#ifndef _MOURISL_COMPACTDS_UTILS
#define _MOURISL_COMPACTDS_UTILS


#include <stdint.h>
#include <stdlib.h>
#include <stdarg.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

namespace compactds {
#define WORD_64 // comment this out if word size is 32

#ifdef WORD_64
  typedef uint64_t WORD ;
  #define WORDBITS 64
  #define WORDBYTES 8
  #define WORDBITS_WIDTH 6 
#else
  typedef uint32_t WORD ;
  #define WORDBITS 32
  #define WORDBYTES 4
  #define WORDBITS_WIDTH 6 
  #define WORDBITS_WIDTH 5
#endif

#define BYTEBITS 8

#define DIV_CEIL(x,y) (((x)%(y))?((x)/(y)+1):((x)/(y)))
#define CEIL(x) (((int)(x) == (x))?((int)(x)):((int)(x) + 1))
#define MIN(x,y) ((x)<=(y)?(x):(y))
#define MAX(x,y) ((x)<=(y)?(y):(x))

// Create a mask of l 1s
#define MASK_WCHECK(l) (((l)>=(int)WORDBITS)?(0xffffffffffffffff):((1ull<<(l))-1ull))
#define MASK(l) ((1ull<<((uint64_t)(l)))-1ull)

const uint64_t MASK_ARRAY[65] =
{
  0ull, 
  0x1ull, 0x3ull, 0x7ull, 0xfull, 0x1full, 0x3full, 0x7full, 0xffull, 
  0x1ffull, 0x3ffull, 0x7ffull, 0xfffull, 0x1fffull, 0x3fffull, 0x7fffull, 0xffffull, 
  0x1ffffull, 0x3ffffull, 0x7ffffull, 0xfffffull, 0x1fffffull, 0x3fffffull, 0x7fffffull, 0xffffffull, 
  0x1ffffffull, 0x3ffffffull, 0x7ffffffull, 0xfffffffull, 0x1fffffffull, 0x3fffffffull, 0x7fffffffull, 0xffffffffull, 
  0x1ffffffffull, 0x3ffffffffull, 0x7ffffffffull, 0xfffffffffull, 0x1fffffffffull, 0x3fffffffffull, 0x7fffffffffull, 0xffffffffffull, 
  0x1ffffffffffull, 0x3ffffffffffull, 0x7ffffffffffull, 0xfffffffffffull, 0x1fffffffffffull, 0x3fffffffffffull, 0x7fffffffffffull, 0xffffffffffffull, 
  0x1ffffffffffffull, 0x3ffffffffffffull, 0x7ffffffffffffull, 0xfffffffffffffull, 0x1fffffffffffffull, 0x3fffffffffffffull, 0x7fffffffffffffull, 0xffffffffffffffull, 
  0x1ffffffffffffffull, 0x3ffffffffffffffull, 0x7ffffffffffffffull, 0xfffffffffffffffull, 0x1fffffffffffffffull, 0x3fffffffffffffffull, 0x7fffffffffffffffull, 0xffffffffffffffffull
} ;

// positive infinity
#define POSITIVE_INF ((uint64_t)-1)

// x-y modules by k-bit block , which are wide words has k bits subblocks 
// h masks/controls the block size
// Sebastiano Vigna, Broadword implementation of rank/select queries, 2008 
#define BITBLOCK_MODDIFF(x,y,h) (((x)|(h)) - ((y)&(~(h)))^(((x)^(~(y))&(h))))
// Test x<y in a subblock fashion
#define BITBLOCK_LT(x,y,h) ((((((x)|(h)) - ((y)&(~(h)))) | ((x)^(y)))^((x)|(~(y))))&(h))
// Test x<= y in a subblock fashion, kind of neg(y<x)
#define BITBLOCK_LEQ(x,y,h) ((((((y)|(h)) - ((x)&(~(h)))) | ((x)^(y)))^((x)&(~(y))))&(h))
// Test x>0 in a subblock fashion
#define BITBLOCK_GZERO(x, l, h) (((((x)|(h))-(l)) | (x)) & (h))

#define SAVE_VAR(fp, x) (fwrite(&(x), sizeof(x), 1, (fp)))
#define LOAD_VAR(fp, x) (fread(&(x), sizeof(x), 1, (fp)))

#define SAVE_ARR(fp, x, n) (fwrite((x), sizeof(*(x)), (n), (fp)))
#define LOAD_ARR(fp, x, n) (fread((x), sizeof(*(x)), (n), (fp)))

#ifdef __GNUC__
  #define CACHE_PREFETCH(x) __builtin_prefetch(x)
#else
  #define CACHE_PREFETCH(x)
#endif

class Utils
{
public:
  // How many bits in the input x 
  static int CountBits(WORD x)
  {
    int ret = 0 ;
    for (; x ; x >>= 1)
      ++ret ;
    return ret ;
  }
  
  // Count the number of 1's in x.
  static int Popcount(WORD x)
  {
#ifdef __GNUC__
      return __builtin_popcountll(x);
#else
#ifdef WORD_64
      x = x - ((x >> 1) & 0x5555555555555555ull) ;
      x = (x&0x3333333333333333ull) + ((x>>2)&0x3333333333333333ull) ;
      return (((x + (x >> 4)) & 0x0f0f0f0f0f0f0f0full) * 0x0101010101010101ull) >> 56 ; 
#else
      x = x - ((x >> 1) & 0x55555555) ;
      x = (x&0x33333333) + ((x>>2)&0x33333333) ;
      return (((x + (x >> 4)) & 0x0f0f0f0f) * 0x01010101) >> 24 ; 
#endif
#endif
    /*else
    {
      int ret = 0 ;
      for (; x ; x &= (x-1))
        ++ret ;
      return ret ;
    }*/
  }

  static int CountTrailingZeros(WORD x)
  {
    if (x == 0)
      return WORDBITS ;
#ifdef __GNUC__
    return __builtin_ctzll(x) ;
#else
    int ret = 0 ;
    for (; !(x & 1ull) ; x >>= 1ull, ++ret)
      ;
    return ret ;
#endif
  }

  // Select the r-th (1-index) 1 in word x
  static int SelectInWord(WORD x, int r)
  {
    const uint64_t l8 = 0x0101010101010101ull ;
    const uint64_t h8 = l8 << 7ull ;
    --r ;

    uint64_t s, b, l ;  
    // Calculate the byte-wise partial sums
    s = x - ((x & 0xAAAAAAAAAAAAAAAAull) >> 1) ;
    s = (s & 0x3333333333333333ull) + ((s>>2) & 0x3333333333333333ull) ;
    s = ((s + (s>>4))&0x0f0f0f0f0f0f0f0full) * l8 ;
    // Locate the byte
    // >> 53 is kind of make the byte unit to bit unit (>>56 << 3), makes later shift easier.
    b = (((BITBLOCK_LEQ(s, r * l8, h8)>>7) * l8) >> 53) & (~7ull) ;
    l = r - (((s<<8) >> b) & 0xff) ; // update remainder
    // Seems the 0x804..01ull trick is to expand the bit information into each byte
    //  each bit in a byte will be in its own byte of a 64bit integer
    s = (BITBLOCK_GZERO(((x >> b & 0xff) * l8 & 0x8040201008040201ull), l8, h8) >> 7) * l8 ;
    return b + (((BITBLOCK_LEQ(s, l * l8, h8) >> 7) * l8) >> 56) ;

  }

  // Compute ceil(log2(x)) without float computation
  static int Log2Ceil(WORD x)
  {
    int bcnt = CountBits(x) ;
    if (x == (1ull<<(bcnt - 1)))
      return bcnt - 1 ;
    else
      return bcnt ;
  }

  // The power function in the integer space
  static uint64_t PowerInt(int x, int y)
  {
    uint64_t ret = 1 ;
    uint64_t powerx = x ;
    while (y)
    {
      if (y & 1)
        ret *= powerx ;
      powerx *= powerx ;
      y >>= 1 ;
    }
    return ret ;
  }

  // The multiplication then take mode: (a*b)%m 
  // that make sure a*b not overflow
  static uint64_t SafeMultiMod(uint64_t a, uint64_t b, uint64_t m)
  {
    uint64_t ret = 0 ;
    a %= m ;
    while (b)
    {
      if (b & 1)
        ret += a%m ;
      a = (a * 2)%m ;
      b >>= 1 ;
    }
    return ret ;
  }
  
  // Assuming only span two words at most.
  //  s is j', e is j
  // Get B[s..e]
  static WORD BitsRead(const WORD *W, const size_t s, const size_t e) 
  {
    // In practice we should let other part be correct about this
    //if (s > e)   
    //  return 0 ;
    const size_t is = s >> WORDBITS_WIDTH ; // index for s
    const int rs = s & (WORDBITS - 1) ;
    const size_t ie = e >> WORDBITS_WIDTH ;
    
    //if (rs + (e - s) <= MASK(WORDBITS_WIDTH))
    //if ((e^s) <= MASK(WORDBITS_WIDTH))
    if (is == ie)
    {
      // in the same block
      return (W[is] >> rs) & MASK_WCHECK(e-s+1) ;
    }
    else
    {
      const int re = e & (WORDBITS - 1) ;// e%w, the residual offset within a word
      // Since ie!=is, re must be less than 63, so we don't need to check the MASK.
      return (W[is] >> rs) | ((W[ie] & MASK(re + 1)) << (WORDBITS - rs)) ;
    }
  }

  // Write B[s..e]=x. 
  static void BitsWrite(WORD *W, size_t s, size_t e, WORD x) 
  {
    if (s > e)
      return ;
    const int w = sizeof(WORD) * 8 ;
    int re = e & (w - 1) ;// e%w, the residual offset within a word
    int rs = s & (w - 1) ;

    size_t ie = e/w ; // index for e
    size_t is = s/w ;

    if (ie == is)
    {
      W[ie] = (W[ie] & ~(MASK_WCHECK(e-s+1) << rs)) | ((WORD)x<<rs); // masking | adding
    }
    else
    {
      W[is] = (W[is] & MASK(rs)) | ((WORD)x << rs) ;
      W[ie] = (W[ie] & ~MASK(re + 1)) | (x >> (w-rs)) ;
    }
  }

  static int BitRead(const WORD *W, size_t i) 
  {
    return (W[i>>WORDBITS_WIDTH] >> (i&(WORDBITS-1)))&1ull ;
  }

  static void BitSet(WORD *W, size_t i)
  {
    W[i>>WORDBITS_WIDTH] |= (1ull << (i&(WORDBITS-1))) ; 
  }
  
  static void BitFlip(WORD *W, size_t i)
  {
    W[i>>WORDBITS_WIDTH] ^= (1ull << (i&(WORDBITS-1))) ; 
  }
  
  static void BitClear(WORD *W, size_t i)
  {
    if (BitRead(W, i))
      W[i>>WORDBITS_WIDTH] -= (1ull << (i&(WORDBITS-1))) ; 
  }
  
  static size_t BitsToWordBytes(size_t l)
  {
    return sizeof(WORD) * DIV_CEIL(l, sizeof(WORD)*8) ;
  }

  static size_t BitsToWords(size_t l)
  {
    return DIV_CEIL(l, sizeof(WORD)*8) ;
  }

  static WORD *MallocByBits(size_t l)
  {
    return (WORD *)calloc(BitsToWords(l), sizeof(WORD)) ;
  }
  
  // Translate the space usage description (TB, GB, MB, KB) to bytes
  static size_t SpaceStringToBytes(const char *s) 
  {
    int i ;
    size_t ret = 0 ;
    for (i = 0 ; s[i] >= '0' && s[i] <= '9' ; ++i)
      ret = ret * 10 + s[i] - '0' ;

    switch (s[i])
    {
      case 'T':
      case 't':
        ret *= 1000000000000ull ; break ;
      case 'G':
      case 'g':
        ret *= 1000000000ull ; break ;
      case 'M':
      case 'm':
        ret *= 1000000ull ; break ;
      case 'K':
      case 'k':
        ret *= 1000ull ; break ;
    }

    return ret ;
  }

  // Deprive the path and extension from the file name
  // Further go one more extension if there is an extra extension 
  //  matching the given extraExtension string
  //  Can use "|" to add multiple potential extra extensions, like "fa|fna|faa"
  // return: length of the base name
  static int GetFileBaseName(const char *s, const char *extraExtension, char *result)
  {
    int i, j, k ;
    int len = strlen(s) ;
    int start, end ;

    start = 0 ;
    for (i = 0 ; i < len ; ++i)
      if (s[i] == '/')
        start = i + 1 ;
    
    for (j = len - 1 ; j >= start ; --j )
    {
      if (s[j] == '.')
        break ;
    }

    if (j < start)
      j = len - 1 ;
    else
      --j ;
    end = j ;
    if (extraExtension != NULL)
    {
      for (i = 0 ; i < extraExtension[i] ; )  
      {
        for (j = i + 1 ; extraExtension[j] && extraExtension[j] != '|' ; ++j )
          ;
        // [i..j) corresponds to one extension in the |-separated extensions
        bool flag = false ;
        for (k = j - 1 ; k >= i ; --k)
          if (end - (j - 1 - k ) < start || extraExtension[k] != s[end - (j - 1 - k)])
          {
            flag = true ;
            break ;
          }
        if (end - (j - 1 - k ) < start || s[end - (j - 1 - k)] != '.')
          flag = true ;

        if (!flag)
        {
          end = end - (j - 1 - k) - 1 ; // extra minus 1 to get rid of the ".".
          break ;
        }

        if (extraExtension[j] == '\0')
          break ;
        i = j + 1 ;
      }
    }

    for (k = start ; k <= end ; ++k)
      result[k - start] = s[k] ;
    result[k - start] = '\0' ;
    return end - start + 1 ;
  }

  static void PrintLog( const char *fmt, ... )
  {
    va_list args ;
    va_start( args, fmt ) ;
    char buffer[500] ;
    vsprintf( buffer, fmt, args ) ;

    time_t mytime = time(NULL) ;
    struct tm *localT = localtime( &mytime ) ;
    char stime[500] ;
    strftime( stime, sizeof( stime ), "%c", localT ) ;
    fprintf( stderr, "[%s] %s\n", stime, buffer ) ;
  }
} ;
}
#endif

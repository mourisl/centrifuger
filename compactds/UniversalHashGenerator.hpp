#ifndef _MOURISL_COMPACTDS_UNIVERSALHASHGENERATOR
#define _MOURISL_COMPACTDS_UNIVERSALHASHGENERATOR

#include "Utils.hpp"

// Universal hash family of ((a*x+b)%p)%m

namespace compactds {
class UniversalHashGenerator
{
private:
  const uint64_t p ; // The largest prime in 63bit, so 2*p can be in 64bit
  uint32_t state ;
  uint64_t m ;

  // Lehmer random generator
  // https://en.wikipedia.org/wiki/Lehmer_random_number_generator
  uint32_t Random()
  {
    return state = (state * 279470273ull) % 0xfffffffb;
  }

  // Not really 64 bit due to 0xfffffffb, but close enough
  uint64_t Random64()
  {
    uint32_t lower32 = Random() ;
    uint32_t upper32 = Random() ;
    return (upper32 * 0xfffffffbull) + lower32 ; 
  }
public:
  UniversalHashGenerator():p(9223372036854775783ull) {}
  ~UniversalHashGenerator() {}
  
  size_t GetSpace() {return sizeof(*this);}

  // map [0..n] to the range of [0,..,m-1]
  // @return: the big prime p. 0 if failed
  uint64_t Init(uint64_t m, uint32_t seed)  
  {
    if (seed == 0) 
      seed = 17 ;
    state = seed ;
    this->m = m ;

    return p ;
  }

  size_t GetP()
  {
    return p;
  }
  
  void SetSeed(uint32_t seed)
  {
    state = seed ;
  }

  // Generate a pair of (a, b)
  void Generate(uint64_t &a, uint64_t &b)
  {
    a = Random64() ;
    if (a == 0)
    {
      state = 17 ;
      a = Random64() ;
    }
    b = Random64() ;
  }

  // Though the outside program should have enough information
  // to do the mapping on its own, we provide the function here
  // for convenience.
  // The function should handle the overflow of a*x
  uint64_t Map(uint64_t a, uint64_t b, uint64_t x)
  {
    return (Utils::SafeMultiMod(x, a, p) + b)%p % m ;
  }
} ;
}

#endif

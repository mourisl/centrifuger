#ifndef _MOURISL_COMPACTDS_PERFECTHASH
#define _MOURISL_COMPACTDS_PERFECTHASH

// Generate a perfect hash function given the set of keys
#include "UniversalHashGenerator.hpp"
#include "FractionBitElemArray.hpp"
#include "Bitvector_Plain.hpp"
#include "SimpleVector.hpp"

#define PERFECT_MAP_KEY_TRIES 3

namespace compactds {
class PerfectHash
{
private:
  UniversalHashGenerator uh ;
  uint64_t a[PERFECT_MAP_KEY_TRIES], b[PERFECT_MAP_KEY_TRIES] ; // the parameters from the universal hash function
  FractionBitElemArray G ;
  size_t m ;

  // Map with hash, this include the shift
  uint64_t MapWithHashI(uint64_t key, int i)
  {
    return uh.Map(a[i], b[i], key) + i * (m/PERFECT_MAP_KEY_TRIES) ;
  }

  // The method is to give each key three potential slots,
  // the goal is to find a map that each slot is assigned by a unique key (one of the three).
  // So we process all the keys first, and start from the slot with unique key already
  // and release the assignment from the other two slots of the key.
  // This may release more slots with unique keys, and we repeat this process
  // If there are still ambiguous keys, we return FAIL(0)
  //
  // I think my implementation is better than the one suggested in the textbook,
  //  as it does not need to store the tuple and nodes/link map and also use queue instead 
  //  of priority_queue which also gives linear time speed
  int InitTry(uint64_t *keys, size_t n, SimpleVector<int> *L, size_t *nL,
      size_t *uniqueSlotQueue, WORD *keyIdxProcessed, size_t *S)
  {
    size_t i ;
    int j ; // in this function, j is to iterate hash function tries
    size_t Scnt = 0 ;
    size_t uniqueSlotQueueS, uniqueSlotQueueE ;

    // Initialize some parametrs.
    for (i = 0 ; i < m ; ++i)
      L[i].Clear() ;
    memset(keyIdxProcessed, 0, Utils::BitsToWordBytes(m)) ;
    uniqueSlotQueueS = 0 ; uniqueSlotQueueE = 0 ; //[S..E)

    for (j = 0 ; j < PERFECT_MAP_KEY_TRIES ; ++j)
      uh.Generate(a[j], b[j]) ;

    // Put all the keys to their slots
    for (i = 0 ; i < n ; ++i)
    {
      for (j = 0 ; j < PERFECT_MAP_KEY_TRIES ; ++j)
      {
        uint64_t target = MapWithHashI(keys[i], j);
        //printf("%llu %llu %d: %d %d: %d\n", a[j], b[j], m, i, j, target) ;
        L[target].PushBack(i) ;
      }
    }
    
    // Initialize the unique slot queue
    for (i = 0 ; i < m ; ++i)
    {
      if (L[i].Size() == 1)
      {
        uniqueSlotQueue[uniqueSlotQueueE] = i ;
        ++uniqueSlotQueueE ;
      }
      nL[i] = L[i].Size() ;
    }

    // main part, identify which slot is unique for a key until now
    while (uniqueSlotQueueS < uniqueSlotQueueE)
    {
      size_t slot = uniqueSlotQueue[uniqueSlotQueueS] ;
      ++uniqueSlotQueueS ;
      // Since each slot will be removed once
      // and the total length of the list PERFECT_MAP_KEY_TRIES*n 
      // , the overall time is still O(n)
      size_t size = L[slot].Size() ;
      size_t keyIdx = -1;
      for (i = 0 ; i < size ; ++i)
      {
        if (!Utils::BitRead(keyIdxProcessed, L[slot][i]))
        {
          keyIdx = L[slot][i] ;
          break ;
        }
      }
      if (i >= size)
      {
        // The l becomes empty, this could happen when a key
        // creates more than one unique-mapped slots
        continue ;
      }
      Utils::BitSet(keyIdxProcessed, keyIdx) ;
      S[Scnt] = keys[keyIdx] ;
      ++Scnt ;
      for (j = 0 ; j < PERFECT_MAP_KEY_TRIES ; ++j)
      {
        uint64_t target = MapWithHashI(keys[keyIdx], j) ; 
        --nL[target] ;
        if (nL[target] == 1) // it could be 0, so we should not use <=1
        {
          uniqueSlotQueue[uniqueSlotQueueE] = target ;
          ++uniqueSlotQueueE ;
        }
      }
    }
    if (Scnt < n)
      return 0 ;
    G.Malloc(3, m) ; // The value of G is {0, 1, 2}
    WORD *V = Utils::MallocByBits(m) ;
    for (i = 1 ; i <= Scnt ; ++i)
    {
      size_t key = keys[S[Scnt - i]] ;  
      uint64_t targets[PERFECT_MAP_KEY_TRIES] ;
      int gSumMod = 0 ;
      for (j = 0 ; j < PERFECT_MAP_KEY_TRIES ; ++j)
      {
        targets[j] = MapWithHashI(key, j) ;
        gSumMod += G.Read(targets[j]) % PERFECT_MAP_KEY_TRIES ;
      }
      for (j = 0 ; j < PERFECT_MAP_KEY_TRIES ; ++j)
      {
        if (!Utils::BitRead(V, targets[j]))
        {
          int tmp = (j - gSumMod) % PERFECT_MAP_KEY_TRIES ;
          if (tmp < 0)
            tmp += PERFECT_MAP_KEY_TRIES ;
          G.Write(targets[j], tmp) ;
          break ;
        }
      }

      for (j = 0 ; j < PERFECT_MAP_KEY_TRIES ; ++j)
        Utils::BitSet(V, targets[j]) ;
    }
    free(V) ;
    return 1 ;
  }
public:
  PerfectHash() {}
  ~PerfectHash() {}
  
  size_t GetSpace() 
  {
    return G.GetSpace() - sizeof(G) + sizeof(*this) ;
  }

  void Init(uint64_t *keys, size_t n, size_t m)  
  {
    if (m == 0)
      m = CEIL(1.25 * n / PERFECT_MAP_KEY_TRIES) * PERFECT_MAP_KEY_TRIES; 
    this->m = m ;
    SimpleVector<int> *L ; // the key list associated with each slot 
    size_t *nL ; // number of element in each L
    size_t *uniqueSlotQueue ; // the queue for slot with unique keys     
    size_t *S ; // the stack used to store keys 
    WORD *keyIdxProcessed ; // bit vector represent whether a key has been processed

    L = new SimpleVector<int>[m] ;
    nL = (size_t *)malloc(sizeof(size_t) * m) ;
    uniqueSlotQueue = (size_t *)malloc(sizeof(size_t) * m) ;
    S = (size_t *)malloc(sizeof(size_t) * n) ;
    keyIdxProcessed = Utils::MallocByBits(n) ;

    uh.Init(m/PERFECT_MAP_KEY_TRIES, 0) ;

    while (!InitTry(keys, n, L, nL, uniqueSlotQueue, keyIdxProcessed, S))
      ;

    delete[] L ;
    free(nL) ;
    free(uniqueSlotQueue) ;
    free(S) ;
    free(keyIdxProcessed) ;
  }

  uint64_t Map(uint64_t x)
  {
    size_t i ;
    uint64_t hs[PERFECT_MAP_KEY_TRIES] ;
    int gsum = 0 ;
    for (i = 0 ; i < PERFECT_MAP_KEY_TRIES ; ++i)
    {
      hs[i] = MapWithHashI(x, i) ; 
      gsum += G.Read(hs[i]) ;
    }
    return hs[gsum %PERFECT_MAP_KEY_TRIES] ;
  }
} ;
}

#endif

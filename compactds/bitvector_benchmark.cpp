#include <iostream>
#include <random>
#include <vector>
#include <chrono>
#include <functional>
#include "Bitvector_Plain.hpp"
#include "Bitvector_Sparse.hpp"
#include "DS_Select_Test.hpp"

using namespace std ;
using namespace std::chrono ;
using timer = std::chrono::high_resolution_clock;
using namespace compactds ; 

const int n = 800000000 ;
const int reps = 10000000 ;

void set_random_bits(std::vector<size_t> &v, int seed)
{
  std::mt19937_64 rng;
  if (0 == seed) {
    rng.seed(std::chrono::system_clock::now().time_since_epoch().count());
  } else
    rng.seed(seed);
  
  size_t *data = v.data() ;
  size_t size = v.size() ;
  *data = rng();
  for (size_t i=1; i < size; ++i) {
    *(++data) = rng();
  }
}

std::vector<size_t> rnd_positions(uint8_t log_s, uint64_t& mask, uint64_t mod=0, uint64_t seed=17)
{
  mask = (1<<log_s)-1;
  std::vector<size_t> rands(1<<log_s ,0);
  set_random_bits(rands, seed);
  if (mod > 0) {
    size_t i ;
    size_t size = rands.size() ;
    for (i = 0 ; i < size ; ++i)
      rands[i] %= mod ;
  }
  return rands;
}


int main(int argc, char *argv[])
{
  size_t i ;
  auto start = timer::now();
  Bitvector_Plain bv ;
  WORD *b = Utils::MallocByBits(n) ;

  std::mt19937_64 rng;
  std::uniform_int_distribution<uint64_t> distribution(0, n-1);
  auto dice = bind(distribution, rng);
  
  // populate vectors with some other bits
  for (i=0; i < n/25; ++i) {
    uint64_t x = dice();
    Utils::BitSet(b, x) ;
  }
  auto stop = timer::now();
  cout << "initialization in (ms): " << duration_cast<milliseconds>(stop-start).count() << endl;

  cout << "size in byptes: " << Utils::BitsToWordBytes(n) << endl;
  
  start = timer::now();
  bv.SetSelectTypeSupport(3) ;
  if (argc == 1)
    bv.SetSelectSpeed(3) ;
  else
    bv.SetSelectSpeed(atoi(argv[1])) ;
  DS_Rank9 ranktst ;
  bv.Init(b, n) ;
  ranktst.Init(b, n) ;
  DS_Select_Test selectTst ;
  selectTst.Init(0, b, n, 2, 3) ;
  stop = timer::now() ;
  cout << "construction in (ms): " << duration_cast<milliseconds>(stop-start).count() << endl;
  cout << "size in bytes: " << bv.GetSpace() << endl;

  auto ones = bv.Rank(1, n) ; 
  auto zeros = n-ones;
  if (0)
  {
    uint64_t mask = 0;
    
    auto rands = rnd_positions(20, mask, zeros, 17);
    for (uint64_t i=0; i<rands.size(); ++i) rands[i] = rands[i]+1;
    start = timer::now();
    uint64_t check = 0 ;
    for (i = 0 ; i < reps ; ++i)
      check += bv.Select(0, rands[i&mask]) ;
    stop = timer::now();

    cout << "# select0_time = " << duration_cast<nanoseconds>(stop-start).count()/(double)reps << endl;
    cout << "# select0_check = " << check << endl;
  }
  if (1)
  {
    uint64_t mask = 0;
    
    auto rands = rnd_positions(20, mask, ones, 17);
    for (uint64_t i=0; i<rands.size(); ++i) rands[i] = rands[i]+1;
    start = timer::now();
    uint64_t check = 0 ;
    for (i = 0 ; i < reps ; ++i)
      check += bv.Select(1, rands[i&mask]) ;
    stop = timer::now();

    cout << "# select1_time = " << duration_cast<nanoseconds>(stop-start).count()/(double)reps << endl;
    cout << "# select1_check = " << check << endl;
  }
  if (1)
  {
    uint64_t mask = 0;
    
    auto rands = rnd_positions(20, mask, ones, 17);
    for (uint64_t i=0; i<rands.size(); ++i) rands[i] = rands[i]+1;
    start = timer::now();
    uint64_t check = 0 ;
    for (i = 0 ; i < reps ; ++i)
      check += selectTst.Query(rands[i&mask], ranktst, b, n) ;
    stop = timer::now();

    cout << "# select1_time = " << duration_cast<nanoseconds>(stop-start).count()/(double)reps << endl;
    cout << "# select1_check = " << check << endl;
    cout << "# select1_space = " << selectTst.GetSpace() << endl;
  }
  {
    uint64_t mask = 0;
    
    auto rands = rnd_positions(20, mask, n, 17);
    start = timer::now();
    uint64_t check = 0 ;
    for (i = 0 ; i < reps ; ++i)
      check += bv.Rank(1, rands[i&mask]) ;
    stop = timer::now();

    cout << "# rank1_time = " << duration_cast<nanoseconds>(stop-start).count()/(double)reps << endl;
    cout << "# rank1_check = " << check << endl;
  }
  
  {
    uint64_t mask = 0;
    
    auto rands = rnd_positions(20, mask, n, 17);
    start = timer::now();
    uint64_t check = 0 ;
    for (i = 0 ; i < reps ; ++i)
      check += ranktst.Query(rands[i&mask], b, n) ;
    stop = timer::now();

    cout << "# rank1_time = " << duration_cast<nanoseconds>(stop-start).count()/(double)reps << endl;
    cout << "# rank1_check = " << check << endl;
    cout << "# rank1_test_space = " << ranktst.GetSpace() << endl;
  }
  return 0 ;
}

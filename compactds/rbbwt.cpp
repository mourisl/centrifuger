#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdarg.h>

#include <string>
#include <fstream>
#include <iostream>
#include <chrono>

#include "SequenceCompactor.hpp"
#include "Sequence_WaveletTree.hpp"
#include "Sequence_RunLength.hpp"
#include "Sequence_Hybrid.hpp"
#include "Sequence_RunBlock.hpp"
#include "FMBuilder.hpp"
#include "FMIndex.hpp"

// Usage: ./a.out (sequence_file|bwt) [load]
using namespace std::chrono ;
using timer = std::chrono::high_resolution_clock;

using namespace compactds ;

int main(int argc, char *argv[])
{
  std::string seq ;
  FixedSizeElemArray s ;
  
  char abList[] = "ACGT" ;
  FixedSizeElemArray BWT ;
  size_t n = 0 ;
  const size_t maxTestCnt = 10000000 ;

  if (atoi(argv[2]) == 0 || argc <= 2)
  {
    std::ifstream ifs(argv[1], std::ifstream::in) ;
    std::getline(ifs, seq) ;
    SequenceCompactor seqCompactor ;
    seqCompactor.Init(abList, s, 1000000) ;
    seqCompactor.Compact(seq.c_str(), s) ;

    n = s.GetSize() ;
    struct _FMBuilderParam param ;
    struct _FMIndexAuxData fmAuxData ;
    param.threadCnt = 4 ;
    param.saBlockSize = n / param.threadCnt ;
    
    FMBuilder::InferParametersGivenMemory(n, strlen(abList), Utils::SpaceStringToBytes("24G"), param) ;
    size_t firstISA = 0 ;
    FMBuilder::Build(s, n, strlen(abList),
        BWT, firstISA, param) ;
    param.Free() ;
    FILE *fp = fopen("tmp.idx", "w") ;
    BWT.Save(fp) ;
    fclose(fp) ;
  }
  else
  {
    FILE *fp = fopen(argv[1], "r") ;
    BWT.Load(fp) ;
    fclose(fp) ;
    
    n = BWT.GetSize() ;
  }
  printf("Total size: %lu\n", n) ;
  
  {
    Sequence_WaveletTree<> plbwt ; // plain bwt
    plbwt.SetSelectSpeed(0) ;
    plbwt.Init(BWT, n, abList) ;
    printf("Plain bwt space (bytes): %lu\n", plbwt.GetSpace()) ;
  
    auto start = timer::now();
    size_t check = 0 ;
    size_t i ;
    for (i = 0 ; i < n && i < maxTestCnt ; ++i)
    {
      size_t x = plbwt.Rank('A', i) ;
      check += x ;
    }
    auto stop = timer::now();
    std::cout << "# rank time (ns) from " << i << "  = " << duration_cast<nanoseconds>(stop-start).count()/(double)i << std::endl;
    std::cout << "# rank sum = " << check << std::endl;
  }
  
  {
    Sequence_RunLength rlbwt ;
    rlbwt.Init(BWT, n, abList) ;
    rlbwt.PrintStats() ;
    printf("Runlength bwt space (bytes): %lu\n", rlbwt.GetSpace()) ;
    
    auto start = timer::now();
    size_t check = 0 ;
    size_t i ;
    for (i = 0 ; i < n && i < maxTestCnt ; ++i)
    {
      size_t x = rlbwt.Rank('A', i) ;
      check += x ;
    }
    auto stop = timer::now();
    std::cout << "# rank time (ns) from " << i << "  = " << duration_cast<nanoseconds>(stop-start).count()/(double)i << std::endl;
    std::cout << "# rank sum = " << check << std::endl;
  }
  
  if (1) 
  {
    Sequence_Hybrid hybbwt ;
    //hybbwt.SetBlockSize(8) ;
    hybbwt.Init(BWT, n, abList) ;
    hybbwt.PrintStats() ;
    printf("Hybrid bwt space (bytes): %lu\n", hybbwt.GetSpace()) ;

    auto start = timer::now();
    size_t check = 0 ;
    size_t i ;
    for (i = 0 ; i < n && i < maxTestCnt ; ++i)
    {
      size_t x = hybbwt.Rank('A', i) ;
      check += x ;
    }
    auto stop = timer::now();
    std::cout << "# rank time (ns) from " << i << "  = " << duration_cast<nanoseconds>(stop-start).count()/(double)i << std::endl;
    std::cout << "# rank sum = " << check << std::endl;
  }

  {
    Sequence_RunBlock rbbwt ;
    //rbbwt.SetBlockSize(5) ;
    rbbwt.Init(BWT, n, abList) ;
    rbbwt.PrintStats() ;
    printf("RunBlock bwt space (bytes): %lu\n", rbbwt.GetSpace()) ;

    auto start = timer::now();
    size_t check = 0 ;
    size_t i ;
    for (i = 0 ; i < n && i < maxTestCnt ; ++i)
    {
      size_t x = rbbwt.Rank('A', i) ;
      check += x ;
    }
    auto stop = timer::now();
    std::cout << "# rank time (ns) from " << i << "  = " << duration_cast<nanoseconds>(stop-start).count()/(double)i << std::endl;
    std::cout << "# rank sum = " << check << std::endl;
  } 

  return 0 ;
}

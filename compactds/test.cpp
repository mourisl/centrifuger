#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdarg.h>

#include "FixedSizeElemArray.hpp"
#include "FractionBitElemArray.hpp"
#include "VariableSizeElemArray_SampledPointers.hpp"
#include "VariableSizeElemArray_DensePointers.hpp"
#include "InterleavedFixedSizeElemArray.hpp"

#include "Bitvector_Plain.hpp"
#include "Bitvector_Compressed.hpp"
#include "Bitvector_Sparse.hpp"
#include "Bitvector_RunLength.hpp"

#include "Sequence_Plain.hpp"
#include "Sequence_WaveletTree.hpp"
#include "Sequence_RunLength.hpp"
#include "Sequence_Hybrid.hpp"
#include "Sequence_RunBlock.hpp"

#include "PerfectHash.hpp"
#include "PartialSum.hpp"

#include "SuffixArrayGenerator.hpp"
#include "FMBuilder.hpp"
#include "FMIndex.hpp"

#include "DS_InvPermutation.hpp"
#include "Permutation.hpp"
#include "InvertedIndex.hpp"

#include "DS_Parenthesis.hpp"
#include "DS_PatternRankSelect.hpp"

#include "Tree_Plain.hpp"
#include "Tree_LOUDS.hpp"
#include "Tree_BP.hpp"
#include "Tree_DFUDS.hpp"

#include "Tree_Cardinal_Plain.hpp"
#include "Tree_Cardinal_LOUDS.hpp"
#include "Tree_Cardinal_Ordinal.hpp"

#include "Tree_Labeled.hpp"

using namespace compactds ; 

void PrintLog( const char *fmt, ... )
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

int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    fprintf(stderr, "Usage: ./test test_case\n") ;
    exit(1) ;
  }

  size_t i ;
  unsigned int mismatchCnt = 0 ;
  //int array[] = {20, 18, 22, 22, 16, 21, 11, 22, 21, 21, 5, 7, 31, 0, 3} ;
  //int array[] = {0xfffffff} ;
  if (!strcmp(argv[1], "array"))
  {
    /*int array[] = {0, 0xfff, 0, 1, 2, 3, 4, 5, 6, 8, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16} ;*/
    //int array[] = {0, 1, 2} ;
    const int n = 1000 ;
    unsigned int array[n] ;
    for (i = 0 ; i < n ; ++i)
      array[i] = (i * 7 + 3)%3 ; // trits
    unsigned int len = sizeof(array) / sizeof(array[0]) ;
    printf("Raw size: %d\n\n", (int)sizeof(array)) ;

    {
      FixedSizeElemArray fsea ;
      //B.Malloc(5, len) ;
      //for (i = 0 ; i < len ; ++i)
      //  B.Write(i, array[i]) ;
      fsea.InitFromArray(-1, array, len) ;
      mismatchCnt = 0 ;
      printf("Fixed-size element array:\n") ;
      for (i = 0 ; i < len ; ++i)
      {
        if (fsea.Read(i) != array[i])
        {
          ++mismatchCnt ;
        }
      }
      printf("mismatch count: %d\n", mismatchCnt) ;
      printf("Space usage (bytes): %d\n", (int)fsea.GetSpace());
    }
    
    {
      FixedSizeElemArray fsea ;
      //B.Malloc(5, len) ;
      //for (i = 0 ; i < len ; ++i)
      //  B.Write(i, array[i]) ;
      fsea.InitFromArray(-1, array, len) ;
      
      FILE *fp = fopen("tmp.out", "w") ;
      fsea.Save(fp) ;
      fclose(fp) ;
      
      fp = fopen("tmp.out", "r") ;
      fsea.Load(fp) ;
      fclose(fp) ;
      
      mismatchCnt = 0 ;
      printf("\nFixed-size element array load/save:\n") ;
      for (i = 0 ; i < len ; ++i)
      {
        if (fsea.Read(i) != array[i])
        {
          ++mismatchCnt ;
        }
      }
      printf("mismatch count: %d\n", mismatchCnt) ;
      printf("Space usage (bytes): %d\n", (int)fsea.GetSpace());
    }
    
    {
      FixedSizeElemArray fsea ;
      //B.Malloc(5, len) ;
      //for (i = 0 ; i < len ; ++i)
      //  B.Write(i, array[i]) ;
      fsea.Malloc(2, 0);
      fsea.Reserve(5);
      for (i = 0 ; i < len ; ++i)
        fsea.PushBack(array[i]);
      mismatchCnt = 0 ;
      printf("\nFixed-size element array with push back:\n") ;
      for (i = 0 ; i < len ; ++i)
      {
        if (fsea.Read(i) != array[i])
        {
          ++mismatchCnt ;
        }
      }
      printf("mismatch count: %d\n", mismatchCnt) ;
      printf("Space usage (bytes): %d\n", (int)fsea.GetSpace());
    }

    FractionBitElemArray fbea ;
    fbea.InitFromArray(0, array, len) ;
    printf("\nFraction bits element array:\n") ;
    mismatchCnt = 0 ;
    printf("Fixed-size element array:\n") ;
    for (i = 0 ; i < len ; ++i)
    {
      if (fbea.Read(i) != array[i])
        ++mismatchCnt ;
    }
    printf("mismatch count: %d\n", mismatchCnt) ;
    printf("Space usage (bytes): %d\n", (int)fbea.GetSpace());

    int blockSize = -1 ;
    VariableSizeElemArray_SampledPointers vseasp ;
    vseasp.InitFromArray(blockSize, array, len) ;
    printf("\nSampled pointers:\n") ;
    mismatchCnt = 0 ;
    for (i = 0 ; i < len ; ++i)
    {
      if (vseasp.Read(i) != array[i])
        ++mismatchCnt ;
    }
    printf("mismatch count: %d\n", mismatchCnt) ;
    printf("Space usage (bytes): %d\n", (int)vseasp.GetSpace());

    VariableSizeElemArray_DensePointers vseadp ;
    vseadp.InitFromArray(blockSize, array, len) ;
    printf("\nDense pointers:\n") ;
    mismatchCnt = 0 ;
    for (i = 0 ; i < len ; ++i)
    {
      if (vseadp.Read(i) != array[i])
        ++mismatchCnt ;
    }
    printf("mismatch count: %d\n", mismatchCnt) ;
    printf("Space usage (bytes): %d\n", (int)vseadp.GetSpace());

    {
      printf("\nInterleaved array:\n") ;
      ILArray il ;
      int block = 3 ;
      il.Malloc(2, DIV_CEIL(n, block), 2, block - 1) ;
      for (i = 0 ; i < len ; ++i)
      {
        if (i%block == 0)
          il.Write(0, i / block, array[i]) ;
        else
          il.Write(1, i - i / block, array[i]) ;
      }
      mismatchCnt = 0 ;
      for (i = 0 ; i < len ; ++i)
      {
        if (i%block == 0)
        {
          if (il.Read(0, i / block) != array[i])
            ++mismatchCnt ;
        }
        else
        {
          if (il.Read(1, i - i / block) != array[i])
            ++mismatchCnt ;
        }
      }
      printf("mismatch count: %d\n", mismatchCnt) ;
      printf("Space usage (bytes): %d\n", (int)il.GetSpace());
    }
    
    {
      printf("\nInterleaved64 array:\n") ;
      IL64Array il64 ;
      int block = 3 ;
      il64.Malloc(DIV_CEIL(n, block), 2, block - 1) ;
      for (i = 0 ; i < len ; ++i)
      {
        if (i%block == 0)
          il64.Write0(i / block, array[i]) ;
        else
          il64.Write1(i - i / block, array[i]) ;
      }
      mismatchCnt = 0 ;
      for (i = 0 ; i < len ; ++i)
      {
        if (i%block == 0)
        {
          if (il64.Read0(i / block) != array[i])
            ++mismatchCnt ;
        }
        else
        {
          if (il64.Read1(i - i / block) != array[i])
            ++mismatchCnt ;
        }
      }
      printf("mismatch count: %d\n", mismatchCnt) ;
      printf("Space usage (bytes): %d\n", (int)il64.GetSpace());
    }
  }
#if 0 // comment out large chunk of the code for compile efficiency. Remove this in future
  else if (!strcmp(argv[1], "bitvector"))
  {
    int k = 0 ;
    unsigned int sum ;
    WORD *B ;
    size_t n = 1000000 ;

    B = Utils::MallocByBits(n) ;
    
    //for (i = 0 ; i*1 < n ; ++i )
    //  Utils::BitSet(B, i*1) ;
    for (i = 0 ; i < n ; )
    {
      int rlen = 5 ;
      if (argc > 2)
        rlen = atoi(argv[2]) ;
      for (int j = 0 ; j < rlen ; ++j)
        Utils::BitSet(B, i + j) ;
      i += 4 * rlen ;
    }
    /*for (i = 0 ; i < n ; ++i)
    {
      if (rand() & 1)
        Utils::BitSet(B, i) ;
    }*/

    printf("Raw size: %d\n", (int)DIV_CEIL(n, 8)) ;
    
    //------
    {
      PrintLog("Plain bitvector:") ;
      Bitvector_Plain bvp ;
      bvp.SetSelectSpeed(1) ;
      bvp.Init(B, n) ;
      for (i = 0 ; i < n ; ++i)
      {
        if (bvp.Access(i) != Utils::BitRead(B, i))
          ++mismatchCnt ;
        //printf("%d %d\n", bvc.Access(i), array.Read(i)) ;
      }
      printf("Access mismatch count: %d\n", mismatchCnt) ;

      mismatchCnt = 0 ;
      sum = 0 ;
      for (i = 0 ; i < n ; ++i)
      {
        if (Utils::BitRead(B, i) == 1)
          ++sum ;
        if (bvp.Rank(1, i) != sum)
          ++mismatchCnt ;
      }
      printf("Rank mismatch count: %d\n", mismatchCnt) ;

      mismatchCnt = 0 ;
      for (int type = 1 ; type >= 1 ; --type)
      {
        k = 0 ;
        for (i = 0 ; i < n ; ++i)
        {
          if (Utils::BitRead(B, i) == type)
          {
            size_t s = bvp.Select(type, k + 1) ;
            if (s != i)
            {
              ++mismatchCnt ;
              //printf("mismatch %d: %d %d\n", k + 1, s, i) ;
            }
            ++k ;
          }
        }
      }
      printf("Select mismatch count: %d\n", mismatchCnt) ;

      mismatchCnt = 0 ;
      k = 0 ;
      for (i = 0 ; i < n ; ++i)
      {
        for (k = i ; k >= 0 ; --k)
          if (bvp.Access(k) == 1)
            break ;
        if (k < 0)
          break ;
        if ((int)bvp.Pred(i) != k)
          ++mismatchCnt ;

      }
      printf("Pred mismatch count: %d\n", mismatchCnt) ;

      mismatchCnt = 0 ;
      k = 0 ;
      for (i = 0 ; i < n ; ++i)
      {
        for (k = i ; k < (int)n ; ++k)
          if (bvp.Access(k) == 1)
            break ;
        if (k >= (int)n)
          break ;
        if ((int)bvp.Succ(i) != k)
          ++mismatchCnt ;
      }
      printf("Succ mismatch count: %d\n", mismatchCnt) ;

      printf("Space usage (byptes): %d\n\n", (int)bvp.GetSpace()) ;
    }
    

    // ------
    {
      PrintLog("Compressed bitvector:") ;
      Bitvector_Compressed bvc ;
      bvc.Init(B, n) ;  

      mismatchCnt = 0 ;
      for (i = 0 ; i < n ; ++i)
      {
        if (bvc.Access(i) != Utils::BitRead(B, i))
          ++mismatchCnt ;
        //printf("%d %d\n", bvc.Access(i), array.Read(i)) ;
      }
      printf("Access mismatch count: %d\n", mismatchCnt) ;

      mismatchCnt = 0 ;
      sum = 0 ;
      for (i = 0 ; i < n ; ++i)
      {
        if (Utils::BitRead(B, i) == 1)
          ++sum ;
        if (bvc.Rank(1, i/*, inclusive=1*/) != sum)
          ++mismatchCnt ;
        //printf("%d %d\n", bvc.Rank(i), sum) ;
      }
      printf("Rank mismatch count: %d\n", mismatchCnt) ;

      /*mismatchCnt = 0 ;
        k = 0 ;
        for (i = 0 ; i < n ; ++i)
        {
        if (Utils::BitRead(B, i) == 1)
        {
        if (bvc.Select(k + 1) != i)
        ++mismatchCnt ;
        ++k ;
        }
        }
        printf("Select mismatch count: %d\n", mismatchCnt) ;*/
      printf("Space usage (byptes): %d\n\n", (int)bvc.GetSpace()) ;
    }
    //-----
    {
      PrintLog("Sparse bitvector:") ;
      Bitvector_Sparse bvs ;
      bvs.Init(B, n) ;  

      mismatchCnt = 0 ;
      for (i = 0 ; i < n ; ++i)
      {
        if (bvs.Access(i) != Utils::BitRead(B, i))
        {
          ++mismatchCnt ;
          //printf("%d: %d %d\n", i, bvs.Access(i), Utils::BitRead(B, i)) ;
        }
      }
      printf("Access mismatch count: %d\n", mismatchCnt) ;

      mismatchCnt = 0 ;
      sum = 0 ;
      for (i = 0 ; i < n ; ++i)
      {
        if (Utils::BitRead(B, i) == 1)
          ++sum ;
        if (bvs.Rank(1, i) != sum)
        {
          ++mismatchCnt ;
          //printf("compare %d: %d %d\n", i, bvs.Rank(1, i), sum) ;
        }
      }
      printf("Rank mismatch count: %d\n", mismatchCnt) ;

      mismatchCnt = 0 ;
      k = 0 ;
      for (i = 0 ; i < n ; ++i)
      {
        if (Utils::BitRead(B, i) == 1)
        {
          if (bvs.Select(k + 1) != i)
            ++mismatchCnt ;
          ++k ;
        }
      }
      printf("Select mismatch count: %d\n", mismatchCnt) ;
      printf("Space usage (byptes): %d\n\n", (int)bvs.GetSpace()) ;
    }
   
    if (1)
    {
      PrintLog("Sparse bitvector load/save:") ;
      Bitvector_Sparse bvs ;
      bvs.Init(B, n) ;  

      FILE *fp = fopen("tmp.out", "w") ;
      bvs.Save(fp) ;
      fclose(fp) ;
      
      fp = fopen("tmp.out", "r") ;
      bvs.Load(fp) ;
      fclose(fp) ;

      mismatchCnt = 0 ;
      for (i = 0 ; i < n ; ++i)
      {
        if (bvs.Access(i) != Utils::BitRead(B, i))
        {
          ++mismatchCnt ;
          //printf("%d: %d %d\n", i, bvs.Access(i), Utils::BitRead(B, i)) ;
        }
      }
      printf("Access mismatch count: %d\n", mismatchCnt) ;

      mismatchCnt = 0 ;
      sum = 0 ;
      for (i = 0 ; i < n ; ++i)
      {
        if (Utils::BitRead(B, i) == 1)
          ++sum ;
        if (bvs.Rank(1, i) != sum)
        {
          ++mismatchCnt ;
          //printf("compare %d: %d %d\n", i, bvs.Rank(1, i), sum) ;
        }
      }
      printf("Rank mismatch count: %d\n", mismatchCnt) ;

      mismatchCnt = 0 ;
      k = 0 ;
      for (i = 0 ; i < n ; ++i)
      {
        if (Utils::BitRead(B, i) == 1)
        {
          if (bvs.Select(k + 1) != i)
            ++mismatchCnt ;
          ++k ;
        }
      }
      printf("Select mismatch count: %d\n", mismatchCnt) ;
      printf("Space usage (byptes): %d\n\n", (int)bvs.GetSpace()) ;
    } 
    
    //-----
    {
      PrintLog("Run-length bitvector:") ;
      Bitvector_RunLength bvr ;
      bvr.Init(B, n) ;  
      mismatchCnt = 0 ;
      for (i = 0 ; i < n ; ++i)
      {
        if (bvr.Access(i) != Utils::BitRead(B, i))
        {
          ++mismatchCnt ;
        }
      }
      printf("Access mismatch count: %d\n", mismatchCnt) ;

      mismatchCnt = 0 ;
      sum = 0 ;
      for (i = 0 ; i < n ; ++i)
      {
        if (Utils::BitRead(B, i) == 1)
          ++sum ;
        if (bvr.Rank(1, i) != sum)
        {
          ++mismatchCnt ;
          //printf("compare %d: %d %d\n", i, bvr.Rank(1, i), sum) ;
        }
      }
      printf("Rank mismatch count: %d\n", mismatchCnt) ;

      mismatchCnt = 0 ;
      k = 0 ;
      for (i = 0 ; i < n ; ++i)
      {
        if (Utils::BitRead(B, i) == 1)
        {
          if (bvr.Select(k + 1) != i)
          {
            //printf("compare %d: %d %d\n", k, bvr.Select(k + 1), i) ;
            ++mismatchCnt ;
          }
          ++k ;
        }
      }
      printf("Select mismatch count: %d\n", mismatchCnt) ;
      printf("Space usage (byptes): %d\n\n", (int)bvr.GetSpace()) ;
    }
    free(B) ;
  }
  else if (!strcmp(argv[1], "sequence"))
  {
    char abList[] = "ACGT" ;
    Alphabet abCode ;
    abCode.InitFromList(abList, strlen(abList)) ;

    size_t n = 1000000 ;
    FixedSizeElemArray S ;
    S.Malloc(2, n) ;
    
    if (1)
    {
      //FILE *fp = fopen("testdata/bwt_1M.out", "r") ;
      FILE *fp = fopen("testdata/bwt_7M-8M.out", "r") ;
      //FILE *fp = fopen("testdata/bwt_2M.out", "r") ;
      //FILE *fp = fopen("testdata/tmp.out", "r") ;
      int t ;
      for (i = 0 ; i < n ; ++i)
      {
        fscanf(fp, "%d", &t) ;
        S.Write(i, t) ;
      }
      fclose(fp) ;
    }
    else
    {
      /*srand(1) ;
      for (i = 0 ; i < n ; ++i)
      {
        S.Write(i, rand()%4) ;
      }*/
      size_t rlen = 5 ;
      if (argc > 2)
        rlen = atoi(argv[2]) ;
      uint8_t prevc = -1 ;
      for (i = 0 ; i < n ; i += rlen)
      {
        uint8_t c = rand() % 4;
        while (c == prevc)
          c = rand() % 4 ;
        for (size_t j = 0 ; j < rlen ; ++j)
          S.Write(i + j, c) ;
        prevc = c ;
      }
    }

    printf("Raw size: %d\n", (int)S.GetSpace()) ;
    
    if (0)
    {
      printf("\nPlain+Bitvector_Plain\n") ;
      Sequence_Plain<Bitvector_Plain> t ;
      t.SetAlphabet(abCode) ;
      t.Init(S, n, abList ) ;

      mismatchCnt = 0 ;
      for (i = 0 ; i < n ; ++i)
      {
        if (t.Access(i) != abList[S.Read(i)])
        {
          ++mismatchCnt ;
        }
      }
      printf("Access mismatch count: %d\n", mismatchCnt) ;

      size_t j = 0 ;
      mismatchCnt = 0 ;
      for (j = 0 ; abList[j] ; ++j)
      {
        uint64_t sum = 0 ;
        for (i = 0 ; i < n ; ++i)
        {
          if (S.Read(i) == j)
            ++sum ;
          if (t.Rank(abList[j], i) != sum)
            ++mismatchCnt ;
        }
      }
      printf("Rank mismatch count: %d\n", mismatchCnt) ;

      mismatchCnt = 0 ;
      for (j = 0 ; abList[j] ; ++j)
      {
        int cnt = 0 ;
        for (i = 0 ; i < n ; ++i)
          if (S.Read(i) == j)
          {
            ++cnt ;
            if (t.Select(abList[j], cnt) != i)
            {
              //printf("%d: %d %d %d\n", (int)j, cnt, t.Select(abList[j], cnt), i) ;
              ++mismatchCnt ;
            }
          }
      }
      printf("Select mismatch count: %d\n", mismatchCnt) ;

      printf("Space usage (byptes): %d\n\n", (int)t.GetSpace()) ;
    }

    if (1)
    {
      printf("\nWavelet tree + plain bitvector:\n") ;
      Sequence_WaveletTree<> t ;
      t.SetAlphabet(abCode) ;
      t.Init(S, n, abList ) ;

      mismatchCnt = 0 ;
      for (i = 0 ; i < n ; ++i)
      {
        if (t.Access(i) != abList[S.Read(i)])
        {
          ++mismatchCnt ;
        }
      }
      printf("Access mismatch count: %d\n", mismatchCnt) ;

      size_t j = 0 ;
      mismatchCnt = 0 ;
      for (j = 0 ; abList[j] ; ++j)
      {
        uint64_t sum = 0 ;
        for (i = 0 ; i < n ; ++i)
        {
          if (S.Read(i) == j)
            ++sum ;
          if (t.Rank(abList[j], i) != sum)
            ++mismatchCnt ;
        }
      }
      printf("Rank mismatch count: %d\n", mismatchCnt) ;

      mismatchCnt = 0 ;
      for (j = 0 ; abList[j] ; ++j)
      {
        int cnt = 0 ;
        for (i = 0 ; i < n ; ++i)
          if (S.Read(i) == j)
          {
            ++cnt ;
            //printf("%d: %d %d %d\n", j, cnt, t.Select(abList[j], cnt), i) ;
            if (t.Select(abList[j], cnt) != i)
            {
              ++mismatchCnt ;
            }
          }
      }
      printf("Select mismatch count: %d\n", mismatchCnt) ;

      printf("Space usage (byptes): %d\n\n", (int)t.GetSpace()) ;
    }
   
    if (0)
    {
      printf("\nsave/load:\n") ;
      Sequence_Hybrid t ;
      t.SetAlphabet(abCode) ;
      t.Init(S, n, abList) ;
      
      FILE *fp = fopen("tmp.out", "w") ;
      t.Save(fp) ;
      fclose(fp) ;
      
      fp = fopen("tmp.out", "r") ;
      t.Load(fp) ;
      fclose(fp) ;

      mismatchCnt = 0 ;
      for (i = 0 ; i < n ; ++i)
      {
        if (t.Access(i) != abList[S.Read(i)])
        {
          ++mismatchCnt ;
        }
      }
      printf("Access mismatch count: %d\n", mismatchCnt) ;

      size_t j = 0 ;
      mismatchCnt = 0 ;
      for (j = 0 ; abList[j] ; ++j)
      {
        uint64_t sum = 0 ;
        for (i = 0 ; i < n ; ++i)
        {
          if (S.Read(i) == j)
            ++sum ;
          if (t.Rank(abList[j], i) != sum)
            ++mismatchCnt ;
        }
      }
      printf("Rank mismatch count: %d\n", mismatchCnt) ;

      /*mismatchCnt = 0 ;
      for (j = 0 ; abList[j] ; ++j)
      {
        int cnt = 0 ;
        for (i = 0 ; i < n ; ++i)
          if (S.Read(i) == j)
          {
            ++cnt ;
            //printf("%d: %d %d %d\n", j, cnt, t.Select(abList[j], cnt), i) ;
            if (t.Select(abList[j], cnt) != i)
            {
              ++mismatchCnt ;
            }
          }
      }
      printf("Select mismatch count: %d\n", mismatchCnt) ;*/

      printf("Space usage (byptes): %d\n\n", (int)t.GetSpace()) ;
    }
    if (1) 
    {
      printf("\nWavelet tree + run-length bitvector:\n") ;
      Sequence_WaveletTree<Bitvector_RunLength> t ;
      t.SetAlphabet(abCode) ;
      t.Init(S, n, abList ) ;

      mismatchCnt = 0 ;
      for (i = 0 ; i < n ; ++i)
      {
        if (t.Access(i) != abList[S.Read(i)])
        {
          ++mismatchCnt ;
        }
      }
      printf("Access mismatch count: %d\n", mismatchCnt) ;

      size_t j = 0 ;
      mismatchCnt = 0 ;
      for (j = 0 ; abList[j] ; ++j)
      {
        uint64_t sum = 0 ;
        for (i = 0 ; i < n ; ++i)
        {
          if (S.Read(i) == j)
            ++sum ;
          if (t.Rank(abList[j], i) != sum)
            ++mismatchCnt ;
        }
      }
      printf("Rank mismatch count: %d\n", mismatchCnt) ;

      /*mismatchCnt = 0 ;
      for (j = 0 ; abList[j] ; ++j)
      {
        int cnt = 0 ;
        for (i = 0 ; i < n ; ++i)
          if (S.Read(i) == j)
          {
            ++cnt ;
            //printf("%d: %d %d %d\n", j, cnt, t.Select(abList[j], cnt), i) ;
            if (t.Select(abList[j], cnt) != i)
            {
              ++mismatchCnt ;
            }
          }
      }
      printf("Select mismatch count: %d\n", mismatchCnt) ;*/

      printf("Space usage (byptes): %d\n\n", (int)t.GetSpace()) ;
    }
    
    {
      printf("\nRun length:\n") ;
      Sequence_RunLength t ;
      //t.SetAlphabet(abCode) ;
      t.Init(S, n, abList ) ;

      mismatchCnt = 0 ;
      for (i = 0 ; i < n ; ++i)
      {
        if (t.Access(i) != abList[S.Read(i)])
        {
          ++mismatchCnt ;
        }
      }
      printf("Access mismatch count: %d\n", mismatchCnt) ;

      size_t j = 0 ;
      mismatchCnt = 0 ;
      for (j = 0 ; abList[j] ; ++j)
      {
        uint64_t sum = 0 ;
        for (i = 0 ; i < n ; ++i)
        {
          if (S.Read(i) == j)
            ++sum ;
          //printf("%d: %d %d %d\n", j, sum, t.Rank(abList[j], i)) ;
          if (t.Rank(abList[j], i) != sum)
          {
            //printf("ERROR\n") ;
            ++mismatchCnt ;
          }
        }
      }
      printf("Rank mismatch count: %d\n", mismatchCnt) ;

      /*mismatchCnt = 0 ;
      for (j = 0 ; abList[j] ; ++j)
      {
        int cnt = 0 ;
        for (i = 0 ; i < n ; ++i)
          if (S.Read(i) == j)
          {
            ++cnt ;
            //printf("%d: %d %d %d\n", j, cnt, t.Select(abList[j], cnt), i) ;
            if (t.Select(abList[j], cnt) != i)
            {
              ++mismatchCnt ;
            }
          }
      }
      printf("Select mismatch count: %d\n", mismatchCnt) ;*/

      printf("Space usage (byptes): %d\n\n", (int)t.GetSpace()) ;
    }
    
    {
      printf("\nHybrid:\n") ;
      Sequence_Hybrid t ;
      //t.SetAlphabet(abCode) ;
      t.SetBlockSize(8) ;
      t.Init(S, n, abList ) ;

      mismatchCnt = 0 ;
      for (i = 0 ; i < n ; ++i)
      {
        if (t.Access(i) != abList[S.Read(i)])
        {
          ++mismatchCnt ;
        }
      }
      printf("Access mismatch count: %d\n", mismatchCnt) ;

      size_t j = 0 ;
      mismatchCnt = 0 ;
      for (j = 0 ; abList[j] ; ++j)
      {
        uint64_t sum = 0 ;
        for (i = 0 ; i < n ; ++i)
        {
          if (S.Read(i) == j)
            ++sum ;
          //printf("%d: %d %d %d\n", i, j, sum, t.Rank(abList[j], i)) ;
          if (t.Rank(abList[j], i) != sum)
          {
            //printf("ERROR\n") ;
            ++mismatchCnt ;
          }
        }
      }
      printf("Rank mismatch count: %d\n", mismatchCnt) ;

      /*mismatchCnt = 0 ;
      for (j = 0 ; abList[j] ; ++j)
      {
        int cnt = 0 ;
        for (i = 0 ; i < n ; ++i)
          if (S.Read(i) == j)
          {
            ++cnt ;
            //printf("%d: %d %d %d\n", j, cnt, t.Select(abList[j], cnt), i) ;
            if (t.Select(abList[j], cnt) != i)
            {
              ++mismatchCnt ;
            }
          }
      }
      printf("Select mismatch count: %d\n", mismatchCnt) ;*/

      printf("Space usage (byptes): %d\n\n", (int)t.GetSpace()) ;
    }
    {
      printf("\nRunBlock:\n") ;
      Sequence_RunBlock t ;
      //t.SetAlphabet(abCode) ;
      t.Init(S, n, abList ) ;
      
      /*FILE *fp = fopen("tmp.out", "w") ;
      t.Save(fp) ;
      fclose(fp) ;
      
      fp = fopen("tmp.out", "r") ;
      t.Load(fp) ;
      fclose(fp) ;*/


      mismatchCnt = 0 ;
      for (i = 0 ; i < n ; ++i)
      {
        if (t.Access(i) != abList[S.Read(i)])
        {
          ++mismatchCnt ;
        }
      }
      printf("Access mismatch count: %d\n", mismatchCnt) ;

      size_t j = 0 ;
      mismatchCnt = 0 ;
      for (j = 0 ; abList[j] ; ++j)
      {
        uint64_t sum = 0 ;
        for (i = 0 ; i < n ; ++i)
        {
          if (S.Read(i) == j)
            ++sum ;
          //printf("%d: %d %d %d\n", i, j, sum, t.Rank(abList[j], i)) ;
          if (t.Rank(abList[j], i) != sum)
          {
            //printf("ERROR\n") ;
            ++mismatchCnt ;
          }
        }
      }
      printf("Rank mismatch count: %d\n", mismatchCnt) ;

      /*mismatchCnt = 0 ;
        for (j = 0 ; abList[j] ; ++j)
        {
        int cnt = 0 ;
        for (i = 0 ; i < n ; ++i)
        if (S.Read(i) == j)
        {
        ++cnt ;
      //printf("%d: %d %d %d\n", j, cnt, t.Select(abList[j], cnt), i) ;
      if (t.Select(abList[j], cnt) != i)
      {
      ++mismatchCnt ;
      }
      }
      }
      printf("Select mismatch count: %d\n", mismatchCnt) ;*/

      printf("Space usage (byptes): %d\n\n", (int)t.GetSpace()) ;
    }
  }
  else if (!strcmp(argv[1], "hash"))
  {
    const int n = 20 ;
    uint64_t array[n] ;

    for (i = 0 ; i < n ; ++i)
      array[i] = i ;
    UniversalHashGenerator uh ;
    uint64_t a, b ;
    int j ;
    uh.Init(2 * n, /*seed=*/0) ;
    printf("Universal hash:\n") ;
    for (j = 0 ; j < 3 ; ++j)
    {
      uh.Generate(a, b) ;
      printf("Hash%d %llu %llu\n", j, (long long unsigned)a, (long long unsigned)b) ;
      for (i = 0 ; i < n ; ++i)
        printf("%d ", (int)uh.Map(a, b, array[i])) ;
      printf("\n") ;
    }
    printf("\n") ;

    PerfectHash perfhash ;
    perfhash.Init(array, n, /*m=*/0) ;
    printf("Perfect hash:\n") ;
    for (i = 0 ; i < n ; ++i)
    {
      printf("%d ", (int)perfhash.Map(array[i])) ;
    }
    printf("\n") ;
    printf("Space usage (bytes): %d\n", (int)perfhash.GetSpace()) ;
  }
  else if (!strcmp(argv[1], "huffman"))
  {
    const int n = 4 ;
    uint64_t freq[n] = {5, 10, 100, 1};
    HuffmanCode huffmanCode ;
    huffmanCode.InitFromFrequency(freq, n) ;

    printf("Huffman code:\n") ;
    for (i = 0 ; i < n ; ++i)
    {
      int l = 0 ;
      WORD code = huffmanCode.Encode(i, l) ;
      printf("%d %d: %llu %d => %d\n", (int)i, (int)freq[i], (long long unsigned)code, l,
          huffmanCode.Decode(code, l)) ;
    }
    printf("Space usage (bytes): %d\n", (int)huffmanCode.GetSpace()) ;
  }
  else if (!strcmp(argv[1], "partialsum"))
  {
    const int n = 100 ;
    int array[n] ;//= {0, 0, 0};
    for (i = 0 ; i < n ; ++i)
      array[i] = i ;
    array[10] = 0 ;

    PartialSum psum ;
    psum.Init(array, n) ;
    
    printf("Succinct partial sum:\n") ;
    int s = 0 ;
    mismatchCnt = 0 ; 
    for (i = 0 ; i <= n ; ++i)
    {
      if (s != (int)psum.Sum(i))
      {
        ++mismatchCnt ;
      }
      s += array[i] ;
    }
    printf("Sum query mismatch count: %d\n", mismatchCnt) ;
    
    mismatchCnt = 0 ;
    int j ;
    s = 0 ;
    for (i = 0 ; i < n ; ++i)
    {
      for (j = 0 ; j < array[i] ; ++j)
      {
        if (psum.Search(s + j) != i)
          ++mismatchCnt ;
      }
      s += array[i] ; 
    }
    printf("Search mismatch count: %d\n", mismatchCnt) ;

    mismatchCnt = 0 ;
    for (i = 0 ; i < n ; ++i)
    {
      if (psum.AccessValue(i) != array[i])
        ++mismatchCnt ;
    }
    printf("AccessValue mismatch count: %d\n", mismatchCnt) ; 

    printf("Space usage (bytes): %d\n", (int)psum.GetSpace()) ;
  }
  else if (!strcmp(argv[1], "sa"))
  {
    FixedSizeElemArray s ;
    size_t n = 10000 ;
    s.Malloc(2, n) ;
    srand(1) ;
    for (i = 0 ; i < (size_t)n ; ++i)
    {
      s.Write(i, rand() % 4) ;
      //s.Write(i, i % 4) ;
      //printf("%d ", s.Read(i)) ;
    }
    //printf("\n") ;
    /*std::vector<size_t> truth ;
    for (j = 0 ; j <= 3 ; ++j)
    {
      int i ;
      for (i = n - 1 - (n - 1) % 4 + j ; i >= 0 ; i -= 4)
      {
        if (i >= n)
          continue ;
        truth.push_back(i) ;
      }
    }*/
    
    // Check the cuts
    SuffixArrayGenerator saGenerator ;
    size_t cutCnt = saGenerator.Init(s, n, n / 4, /*diffcov_v=*/4096, 4) ;
    /*for (i = 0 ; i < cutCnt ; ++i)
    {
      std::vector< std::vector<size_t> > pos = saGenerator.GetChunksPositions(s, n, i, i) ;
      //printf("%d\n", pos[0].size()) ;
      int size = pos[0].size() ;
      int j ;
      for (j = 0 ; j < size ; ++j)
        printf("%d ", pos[0][j]) ;
      printf("\n") ;
    }*/
    size_t *sa = (size_t *)malloc(sizeof(size_t) * n);
    size_t calculated = 0 ;
    for (i = 0 ; i < cutCnt ; ++i)
    {
      std::vector< std::vector<size_t> > pos ;
      saGenerator.GetChunksPositions(s, n, i, i, 0, n - 1, pos) ;
      int size = pos[0].size() ;
      printf("%lu %d. %lu\n", i, size, calculated) ;
      saGenerator.SortSuffixByPos(s, n, pos[0].data(), size, sa + calculated) ;
      calculated += size ;
    }
    printf("Validate result: %d\n", saGenerator.ValidateSA(s, n, sa)) ;
    free(sa) ;
    
    /*mismatchCnt = 0 ;
    for (i = 0 ; i < (size_t)n ; ++i)
    {
      if (truth[i] != sa[i])
        ++mismatchCnt ;
      //printf("%d %d\n", (int)truth[i], (int)sa[i]) ;
    }
    printf("SA mismatch: %d\n", mismatchCnt) ;*/
    
  }
  else if (!strcmp(argv[1], "fm"))
  {
    FixedSizeElemArray s ;
    const size_t n = 10000 ;
    const size_t testLen = 50 ;
    s.Malloc(2, n) ;
    char strs[n + 1] ;
    srand(1) ;
    char abList[] = "ACGT" ;
    for (i = 0 ; i < (size_t)n ; ++i)
    {
      int r = rand() % 4 ;
      s.Write(i, r) ;
      strs[i] = abList[r] ; 
    }
    //s.Print(stdout) ;
    //printf("%s\n", strs) ;
    struct _FMBuilderParam param ;
    struct _FMIndexAuxData fmAuxData ;
    param.threadCnt = 4 ;
    param.saBlockSize = n / 4 ;
    FixedSizeElemArray BWT ;
    param.precomputeWidth = testLen > 10 ? 10 : testLen ;
    param.maxLcp = 17 ;
    
    size_t firstISA = 0 ;
    param.selectedISA[0] = 0 ; 
    param.selectedISA[1] = 0 ; 
    FMBuilder::Build(s, n, 4, BWT, firstISA, param) ;

    Sequence_RunBlock t ;
    t.Init(BWT, n, abList) ;

    //BWT.Print(stdout) ;
    //
    //printf("%d %d\n", precomputedRange[0].first, precomputedRange[0].second) ;
    /*for (i = 0 ; i < 100 ; ++i)
    {
      if (precomputedRange[i].second > 0)
        printf("%d %d\n", precomputedRange[i].first, precomputedRange[i].second) ;
    }*/
    size_t count = 0 ;
    for (i = 0 ; i < n ; i += WORDBITS)
      count += Utils::Popcount(param.semiLcpGreater[i / WORDBITS]) ;
    printf("Number of 1s in semiLcpGreater: %lu\n", count) ;
    
    count = 0 ;
    for (i = 0 ; i < n ; i += WORDBITS)
      count += Utils::Popcount(param.semiLcpEqual[i / WORDBITS]) ;
    printf("Number of 1s in semiLcpEqual: %lu\n", count) ;
    
    FMIndex< Sequence_WaveletTree<Bitvector_Plain> > fmIndex ;
    //FMIndex< Sequence_Plain<Bitvector_Plain> > fmIndex ;
    //FMIndex< Sequence_RunBlock > fmIndex ;
    fmIndex.Init(BWT, n, firstISA, 
        param,
        abList, strlen(abList)) ;
    printf("firstISA = %lu; lastISA = %lu\n", firstISA, fmIndex.GetLastISA()) ;

    size_t sp, ep, l ;
    char test[testLen + 1] ;
    test[testLen] = '\0' ;
    size_t k ;
    size_t mismatchCnt = 0 ;
    size_t compareCnt = 0 ;
    for (k = 0 ; k + testLen <= n; ++k)
    {
      memcpy(test, strs + k, testLen) ;
      //strcpy(test, "GATGGAGATG") ;
      //printf("test: %s\n", test) ;
      l = fmIndex.BackwardSearch(test, strlen(test), sp, ep) ;
      //printf("Backward search %d %d %d\n", l, sp, ep) ;
      if (sp < ep)
        continue ;
      ++compareCnt ;
      for (i = sp ; i <= ep ; ++i)
      {
        size_t sa = fmIndex.BackwardToSampledSA(i, l) ; 
        if (sa + l != k)
        {
          ++mismatchCnt ;
          printf("SA[%lu] = %lu+%lu. %lu\n", i, sa, l, k) ;
        }
      }
    }
    printf("Mismatch count: %lu out of %lu\n", mismatchCnt, compareCnt) ;

    /*FILE *fp = fopen("tmp.out", "w") ;
    fmIndex.Save(fp) ;
    fclose(fp) ;
    
    fp = fopen("tmp.out", "r") ;
    fmIndex.Load(fp) ;
    fclose(fp) ;

    printf("Save/Load:\n") ;
    l = fmIndex.BackwardSearch(test, strlen(test), sp, ep) ;
    printf("Backward search %d %d %d\n", l, sp, ep) ;
    for (i = sp ; i <= ep ; ++i)
    {
      size_t sa = fmIndex.BackwardToSampledSA(i, l) ; 
      printf("SA[%d] = %d+%d\n", i, sa, l) ;
    }*/

    
    //free(sampledSa) ;
    //free(precomputedRange) ;
    //free(semiLcpGreater) ;
    //free(semiLcpEqual) ;
  }
  else if (!strcmp(argv[1], "diffcover"))
  {
    DifferenceCover dc ;
    unsigned int v = 4096 ;
    dc.Init(v) ;
    size_t j ;
    mismatchCnt = 0 ;
    for (i = 0 ; i < v ; ++i)
      for (j = 0 ; j < v ; ++j)
      {
        int d = dc.Delta(i, j) ;
        if (!dc.IsInDC(i + d) || !dc.IsInDC(j + d))
        {
          ++mismatchCnt ;
        }
      }
    printf("%d\n", mismatchCnt) ;
  }
  else if (!strcmp(argv[1], "permutation"))
  {
    const int n = 1000;
    size_t *perm = new size_t[n];
    size_t *inv = new size_t[n] ;
    for (i = 0 ; i < n ; ++i)
      perm[i] = (i + 1)%n ;
    printf("Raw permutation size %lu\n", n * sizeof(perm[0])) ;
    for (i = 0 ; i < n ; ++i)
    {
      if (0)
      {
        size_t tmp ;
        size_t j = i + rand() % (n - i) ;
        tmp = perm[j] ;
        perm[j] = perm[i] ;
        perm[i] = tmp ;
      }
      //perm[i] = (i*10001+1)%n ;
      inv[ perm[i] ] = i ;
    }
    //for (i = 0 ; i < n ; ++i)
    //  printf("%d ", perm[i]) ;
    //printf("\n") ;

    {
      printf("\ninverse permutation\n") ;
      DS_InvPermutation invperm ;
      invperm.Init(perm, n) ;
      mismatchCnt = 0 ;
      for (i = 0 ; i < n ; ++i)
        if (invperm.Query(perm, i) != inv[i])
          ++mismatchCnt ;
      printf("Inverse mismatch count %d\n", mismatchCnt) ;
      printf("Space usage: %d\n", (int)invperm.GetSpace()) ;
    }

    {
      printf("\ncompressed permutation\n") ;
      Permutation cperm ;
      cperm.Init(perm, n) ;
      
      mismatchCnt = 0 ;
      for (i = 0 ; i < n ; ++i)
      {
        if (cperm.Next(i) != perm[i])
          ++mismatchCnt ;
      }
      printf("Next mismatch count %d\n", mismatchCnt) ;
      
      mismatchCnt = 0 ;
      for (i = 0 ; i < n ; ++i)
      {
        if (cperm.Prev(i) != inv[i])
          ++mismatchCnt ;
      }
      printf("Prev mismatch count %d\n", mismatchCnt) ;
      
      printf("Space usage: %d\n", (int)cperm.GetSpace()) ;
    }

    delete[] perm ;
    delete[] inv ;
  }
  else if (!strcmp(argv[1], "invindex"))
  {
    size_t n = 10000 ;
    FixedSizeElemArray a ;
    a.Malloc(3, n) ;
    size_t i ;
    srand(17) ;
    int stride = 5 ;
    for (i = 0 ; i < n ; ++i)
    {
      a.Write(i, i%stride) ;
    }

    InvertedIndex idx ;
    idx.Init(a, n, false) ;
    size_t mismatchCnt = 0 ;
    printf("Raw sequence space usage %lu\n", a.GetSpace()) ;
    int label = 1 ;
    for (i = label ; i < n ; i += stride)
    {
      if (idx.Search(label, i / stride) != i)
      {
        //printf("%lu %lu\n", i, idx.Search(label, i / stride)) ;
        ++mismatchCnt ;
      }
    }
    printf("Inverted index mismatch count %lu\n", mismatchCnt) ;
    printf("Inverted index based on permutation space usage %lu\n", idx.GetSpace()) ;
  }
#endif
  else if (!strcmp(argv[1], "rmmtree"))
  {
    int n = 1000000 ; 
    int i, j ;
    WORD *B = Utils::MallocByBits(n) ;
    printf("Raw representation space usage: %lu\n", Utils::BitsToWordBytes(n)) ;
    srand(1) ;
    for (i = 0 ; i < n ; ++i)
    {
      if (rand() & 1)
      //if (i < n / 2)
      //if (i % 2 == 0)
        Utils::BitSet(B, i) ;
    }
    
    DS_RangeMinMaxTree rmmTree ;
    rmmTree.SetBlockSize(32) ;
    rmmTree.Init(B, n) ;
    
    // Test forward and backward search
    if (0)
    {
      int d ;
      int stride = 11 ;
      for (d = -stride ; d <= stride ; d += 2 * stride)
      {
        mismatchCnt = 0 ;
        for (i = 0 ; i < n ; ++i)
        {
          int excess = 0 ;
          for (j = i ; j < n ; ++j)
          {
            excess += 2 * Utils::BitRead(B, j) - 1 ;
            if (excess == d)
              break ;
          }
          int truth = j ;
          j = rmmTree.FwdSearch(i, d, B, n) ;
          if (j != truth)
            ++mismatchCnt ;
          //if (j != truth)
          //  printf("%d %d %d\n", i, truth, j) ;
        }
        printf("FwdSearch %d mismatch count %u\n", d, mismatchCnt) ;
      }

      for (d = -stride ; d <= stride ; d += 2 * stride)
      {
        mismatchCnt = 0 ;
        for (i = 0 ; i < n ; ++i)
        {
          int excess = 0 ;
          for (j = i ; j >= 0 ; --j)
          {
            excess -= (2 * Utils::BitRead(B, j) - 1) ;
            if (excess == d)
              break ;
          }
          int truth = j ;
          if (truth == -1)
            truth = n ;
          j = rmmTree.BwdSearch(i, d, B, n) ;
          if (j != truth)
            ++mismatchCnt ;
          //if (j != truth)
          //  printf("%d %d %d\n", i, truth, j) ;
        }
        printf("BwdSearch %d mismatch count %u\n", d, mismatchCnt) ;
      }
    }
    
    // Test rmq and rMq
    {
      int len = 10000 ;
      mismatchCnt = 0 ;
      for (i = 0 ; i + len <= n ; ++i)
      {
        int excess = 0 ;
        int min = 2 ;
        int mintag = i, maxtag = i;
        int max = -2 ;
        int minCnt = 0 ;
        int lastMinTag = 0 ;
        for (j = i ; j < i + len ; ++j)
        {
          excess += (2 * Utils::BitRead(B, j) - 1) ;
          if (excess < min)
          {
            min = excess ;
            mintag = j ;
            minCnt = 1 ;
            lastMinTag = j ;
          }
          else if (excess == min)
          {
            ++minCnt ;
            lastMinTag = j ;
          }

          if (excess > max)
          {
            max = excess ;
            maxtag = j ; 
          }
        }
        
        //printf("%d %d\n", rmmTree.ExtremeExcess(B, n, i, i + len - 1, 0), min) ;
        if (rmmTree.ExtremeExcess(i, i + len - 1, 0, B, n) != min)
        {
          ++mismatchCnt ;
          //printf("min mismatch %d\n", i) ;
        }
        if (rmmTree.ExtremeExcess(i, i + len - 1, 1, B, n) != max)
        {
          ++mismatchCnt ;
          //printf("max mismatch %d\n", i) ;
        }

        if ((int)rmmTree.Rmq(i, i + len - 1, B, n) != mintag)
        {
          ++mismatchCnt ;
        }
        
        if ((int)rmmTree.RMq(i, i + len - 1, B, n) != maxtag)
        {
          ++mismatchCnt ;
        }

        if ((int)rmmTree.MinCount(i, i + len - 1, B, n) != minCnt)
        {
          ++mismatchCnt ;
          //printf("min count mismatch %d: %d %d\n", i, min, minCnt) ;
        }

        if ((int)rmmTree.MinSelect(i, i + len - 1, minCnt, B, n) != lastMinTag)
        {
          ++mismatchCnt ;
          //printf("min select mismatch %d: %d %d %d\n", i, min, minCnt, lastMinTag) ;
        }
      }
      printf("extreme excess mismatch count %u\n", mismatchCnt) ;
    }

    printf("rmmTree space usage (bytes): %lu\n", rmmTree.GetSpace(true)) ;
  }
  else if (!strcmp(argv[1], "tree"))
  {
    // Test example: tree with child count 2 (from root), 3, 4, 5, ....
    // Or a binary tree
    int i, j ;
    Tree_Plain tree ;
    tree.Init() ;
    
    int internalN = 10000 ; 
    srand(1) ;
    for (i = 0 ; i < internalN ; ++i)
    {
      int childCnt = rand() % 4 + 1 ;
      for (j = 0 ; j < childCnt ; ++j)
      {
        size_t tid = tree.AddNode(i) ;
        tree.SetLabel(tid, childCnt) ;
      }
    }

    size_t *map = new size_t[tree.GetSize()] ;
    
    if (0)
    {
      // This test is for binary tree.
      mismatchCnt = 0 ;
      for (i = 0 ; i < internalN ; ++i)
      {
        if (tree.ChildrenCount(i) != 2 || (int)tree.FirstChild(i) != (2 * i + 1) 
            || (int)tree.LastChild(i) != (2 * i + 2))
          ++mismatchCnt ;
      }
      printf("plain tree mismatch count %u\n", mismatchCnt) ;
      //printf("plain tree space usage (bytes): %lu\n", tree.GetSpace(true)) ;
    }
    printf("plain tree space usage (bytes): %lu\n", tree.GetSpace(true)) ;

    {
      Tree_LOUDS t ;
      mismatchCnt = 0 ;
    
      t.Init(tree.GetTreeData().data(), tree.GetSize(), map) ;
      for (i = 0 ; i < internalN ; ++i)
      {
        if (tree.IsLeaf(i))
          continue ;
        size_t v = t.NodeSelect(map[i]) ;
        if (t.ChildrenCount(v) != tree.ChildrenCount(i) 
            || t.NodeMap(t.FirstChild(v)) != map[tree.FirstChild(i)] 
            || t.NodeMap(t.LastChild(v)) != map[tree.LastChild(i)]
            || t.ChildRank(v) != tree.ChildRank(i)
            || (!tree.IsLastChild(i) && t.NodeMap(t.NextSibling(v)) != map[tree.NextSibling(i)])
            || (!tree.IsFirstChild(i) && t.NodeMap(t.PrevSibling(v)) != map[tree.PrevSibling(i)])
            || t.NodeMap(t.Parent(v)) != map[tree.Parent(i)]
            || t.NodeMap(t.LCA(v, t.NodeSelect(map[internalN]))) != map[tree.LCA(i, internalN)]
            )
          ++mismatchCnt ;
      }
      printf("\nLOUDS tree mismatch count %u\n", mismatchCnt) ;
      printf("LOUDS tree space usage (bytes): %lu\n", t.GetSpace()) ;
    }
    
    {
      Tree_BP t ;
      mismatchCnt = 0 ;

      t.Init(tree.GetTreeData().data(), tree.GetSize(), map) ;
      for (i = 0 ; i < internalN ; ++i)
      {
        if (tree.IsLeaf(i))
          continue ;
        size_t v = t.NodeSelect(map[i]) ;
        if (t.ChildrenCount(v) != tree.ChildrenCount(i) 
            || t.NodeMap(t.FirstChild(v)) != map[tree.FirstChild(i)] 
            || t.NodeMap(t.LastChild(v)) != map[tree.LastChild(i)]
            || t.ChildRank(v) != tree.ChildRank(i)
            || t.NodeMap(t.Parent(v)) != map[tree.Parent(i)]
            || t.NodeMap(t.ChildSelect(v, 1)) != map[tree.ChildSelect(i, 1)]
            || (!tree.IsLastChild(i) && t.NodeMap(t.NextSibling(v)) != map[tree.NextSibling(i)])
            || (!tree.IsFirstChild(i) && t.NodeMap(t.PrevSibling(v)) != map[tree.PrevSibling(i)])
            || (t.IsAncestor(v, t.NodeSelect(map[internalN])) != tree.IsAncestor(i, internalN))
            || t.Depth(v) != tree.Depth(i)
            || t.SubTreeSize(v) != tree.SubTreeSize(i)
            || t.LeafCountInSubTree(v) != tree.LeafCountInSubTree(i)
            || t.NodeMap(t.LCA(v, t.NodeSelect(map[internalN]))) != map[tree.LCA(i, internalN)]
           )
        {
          ++mismatchCnt ;
          //printf("%d %d. %d. %d %d. %d\n", t.LeafCountInSubTree(v), tree.LeafCountInSubTree(i), tree.GetSize(), v, i, tree.ChildrenCount(i)) ;
        }
      }
      
      for (i = internalN ; i < (int)tree.GetSize(); ++i)
      {
        size_t v = t.NodeSelect(map[i]) ;
        if (t.LeafCountInSubTree(v) != tree.LeafCountInSubTree(i))
          //|| (int)t.LeafRank(v) != i - internalN + 1)
          ++mismatchCnt ;
      }

      printf("\nBP tree mismatch count %u\n", mismatchCnt) ;
      printf("BP tree space usage (bytes): %lu\n", t.GetSpace()) ;
    }
    
    {
      Tree_DFUDS t ;
      mismatchCnt = 0 ;

      t.Init(tree.GetTreeData().data(), tree.GetSize(), map) ;
      for (i = 0 ; i < internalN ; ++i)
      {
        if (tree.IsLeaf(i))
          continue ;
        size_t v = t.NodeSelect(map[i]) ;
        if (t.ChildrenCount(v) != tree.ChildrenCount(i) 
            || t.NodeMap(t.FirstChild(v)) != map[tree.FirstChild(i)] 
            || t.NodeMap(t.LastChild(v)) != map[tree.LastChild(i)]
            || t.ChildRank(v) != tree.ChildRank(i)
            || t.NodeMap(t.Parent(v)) != map[tree.Parent(i)]
            || t.NodeMap(t.ChildSelect(v, 1)) != map[tree.ChildSelect(i, 1)]
            || (!tree.IsLastChild(i) && t.NodeMap(t.NextSibling(v)) != map[tree.NextSibling(i)])
            || (!tree.IsFirstChild(i) && t.NodeMap(t.PrevSibling(v)) != map[tree.PrevSibling(i)])
            || (t.ChildrenCount(v) > 1 && t.NodeMap(t.ChildSelect(v, 2)) != map[tree.ChildSelect(i, 2)])
            || (t.IsAncestor(v, t.NodeSelect(map[internalN])) != tree.IsAncestor(i, internalN))
            //|| t.Depth(v) != tree.Depth(i)
            || t.SubTreeSize(v) != tree.SubTreeSize(i)
            || t.LeafCountInSubTree(v) != tree.LeafCountInSubTree(i)
            || t.NodeMap(t.LCA(v, t.NodeSelect(map[internalN]))) != map[tree.LCA(i, internalN)]
           )
          ++mismatchCnt ;
        //printf("%d %d %d\n", v,
        //    t.NodeSelect(map[internalN]),
        //    t.LCA(v, t.NodeSelect(map[internalN]))) ;
        //printf("%d %d. %d. %d %d. %d\n", t.LeafCountInSubTree(v), tree.LeafCountInSubTree(i), tree.GetSize(), v, i, tree.ChildrenCount(i)) ;
      }

      for (i = internalN ; i < (int)tree.GetSize(); ++i)
      {
        size_t v = t.NodeSelect(map[i]) ;
        if (t.LeafCountInSubTree(v) != tree.LeafCountInSubTree(i))
          ++mismatchCnt ;
      }
      printf("\nDFUDS tree mismatch count %u\n", mismatchCnt) ;
      printf("DFUDS tree space usage (bytes): %lu\n", t.GetSpace()) ;
    }
    
    {
      Tree_Labeled<> t ;
      mismatchCnt = 0 ;

      t.Init(tree.GetTreeData().data(), tree.GetSize(), map) ;
      for (i = 0 ; i < internalN ; ++i)
      {
        if (tree.IsLeaf(i))
          continue ;
        size_t v = t.NodeSelect(map[i]) ;
        if (t.ChildrenCount(v) != tree.ChildrenCount(i) 
            || t.NodeMap(t.FirstChild(v)) != map[tree.FirstChild(i)] 
            || t.NodeMap(t.LastChild(v)) != map[tree.LastChild(i)]
            || t.ChildRank(v) != tree.ChildRank(i)
            || t.NodeMap(t.Parent(v)) != map[tree.Parent(i)]
            || t.NodeMap(t.ChildSelect(v, 1)) != map[tree.ChildSelect(i, 1)]
            || (!tree.IsLastChild(i) && t.NodeMap(t.NextSibling(v)) != map[tree.NextSibling(i)])
            || (!tree.IsFirstChild(i) && t.NodeMap(t.PrevSibling(v)) != map[tree.PrevSibling(i)])
            || (t.ChildrenCount(v) > 1 && t.NodeMap(t.ChildSelect(v, 2)) != map[tree.ChildSelect(i, 2)])
            || (t.IsAncestor(v, t.NodeSelect(map[internalN])) != tree.IsAncestor(i, internalN))
            //|| t.Depth(v) != tree.Depth(i)
            || t.SubTreeSize(v) != tree.SubTreeSize(i)
            || t.LeafCountInSubTree(v) != tree.LeafCountInSubTree(i)
            || t.NodeMap(t.LCA(v, t.NodeSelect(map[internalN]))) != map[tree.LCA(i, internalN)]
           )
          ++mismatchCnt ;
        //printf("%d %d %d\n", v,
        //    t.NodeSelect(map[internalN]),
        //    t.LCA(v, t.NodeSelect(map[internalN]))) ;
        //printf("%d %d. %d. %d %d. %d\n", t.LeafCountInSubTree(v), tree.LeafCountInSubTree(i), tree.GetSize(), v, i, tree.ChildrenCount(i)) ;
      }

      // Labels
      for (i = 0 ; i < internalN ; ++i)
      {
        size_t v = t.NodeSelect(map[i]) ;
        if (t.ChildrenLabeled(v, 3) != tree.ChildrenLabeled(i, 3)
            || (tree.ChildrenLabeled(i, 3 ) > 0 && t.NodeMap(t.LabeledChildSelect(v, 3, 2)) != map[tree.LabeledChildSelect(i, 3, 2)])
            || t.ChildLabel(v) != tree.ChildLabel(i)
            )
        {
         // printf("%d %d %d: %d %d %d\n", i, map[i], v, t.ChildrenLabeled(v, 1), tree.ChildrenLabeled(i, 1), 
        //      tree.ChildrenCount(i)) ;
          ++mismatchCnt ;
        }
      }
      
      printf("\nLabeled tree mismatch count %u\n", mismatchCnt) ;
      printf("Labeled tree space usage (bytes): %lu\n", t.GetSpace()) ;
    }
    
    delete[] map ;
  }
  else if (!strcmp(argv[1], "patternrs"))
  {
    int k = 0 ;
    unsigned int sum ;
    WORD *B ;
    size_t n = 1000000  ;
    n = 65536*32 ;

    B = Utils::MallocByBits(n) ;

    //for (i = 0 ; i*1 < n ; ++i )
    //  Utils::BitSet(B, i*1) ;
    /*for (i = 0 ; i < n ; )
      {
      int rlen = 5 ;
      if (argc > 2)
      rlen = atoi(argv[2]) ;
      for (int j = 0 ; j < rlen ; ++j)
      Utils::BitSet(B, i + j) ;
      i += 4 * rlen ;
      }*/
    /*for (i = 0 ; i < n ; ++i)
    {
      if (rand() & 1)
        Utils::BitSet(B, i) ;
    }*/
    DS_Parenthesis tmp ;
    tmp.GenerateRandomBalanceParenthesis(B, n) ;

    printf("Raw size: %d\n", (int)DIV_CEIL(n, 8)) ;

    WORD pat = 2 ;  // binary 10
    int patLen = 2 ;
    DS_PatternRankSelect patrs ;
    mismatchCnt = 0 ;
    sum = 0 ;
    patrs.Init(B, n, pat, patLen) ;
    for (i = 0 ; i < n ; ++i)
    {
      if (patrs.IsPattern(i, B, n))
        ++sum ;
      if (patrs.Rank(i, B, n) != sum)
        ++mismatchCnt ;
    }
    printf("Rank mismatch count: %d\n", mismatchCnt) ;

    mismatchCnt = 0 ;
    k = 0 ;
    for (i = 0 ; i < n ; ++i)
    {
      if (patrs.IsPattern(i, B, n))
      {
        size_t s = patrs.Select(k + 1, B, n) ;
        if (s != i)
        {
          ++mismatchCnt ;
          //printf("mismatch %d: %d %d\n", k + 1, s, i) ;
        }
        ++k ;
      }
    }
    printf("Select mismatch count: %d (%d)\n", mismatchCnt, k) ;
    printf("DS_PatternRankSelect space: %lu\n", patrs.GetSpace()) ;
  }
  else if (!strcmp(argv[1], "cardtree")) // cardinal tree
  {
    // Test example: tree with child count 2 (from root), 3, 4, 5, ....
    // Or a binary tree
    int i, j ;
    int c = 4 ; // cardinality

    Tree_Cardinal_Plain tree ;
    tree.Init(c) ;
    
    int internalN = 10000 ; 
    srand(1) ;
    for (i = 0 ; i < internalN ; ++i)
    {
      //int childCnt = rand() % c + 2 ;
      //int childCnt = c ;
      int step = rand() % c + 1 ;
      //step = 1 ;
      for (j = 0 ; j < c ; j += step)
        tree.AddNode(i, j) ;
    }
    
    size_t *map = new size_t[tree.GetSize()] ;
    if (0)
    {
      // This test is for binary tree.
      mismatchCnt = 0 ;
      for (i = 0 ; i < internalN ; ++i)
      {
        if (tree.ChildrenCount(i) != 2 || (int)tree.FirstChild(i) != (2 * i + 1) 
            || (int)tree.LastChild(i) != (2 * i + 2))
          ++mismatchCnt ;
      }
      printf("plain cardinal tree mismatch count %u\n", mismatchCnt) ;
    }
    printf("plain cardinal tree space usage (bytes): %lu\n", tree.GetSpace(true)) ;
    
    if (1)
    {
      Tree_Cardinal_LOUDS<> t ;
      mismatchCnt = 0 ;

      t.Init(tree.GetTreeData().data(), tree.GetSize(), c, map) ;
      for (i = 0 ; i < internalN ; ++i)
      {
        size_t v = t.NodeSelect(map[i]) ;
        if (t.ChildrenCount(v) != tree.ChildrenCount(i) 
            || t.NodeMap(t.FirstChild(v)) != map[tree.FirstChild(i)] 
            || t.NodeMap(t.LastChild(v)) != map[tree.LastChild(i)]
            || t.NodeMap(t.Parent(v)) != map[tree.Parent(i)]
            || (t.ChildrenCount(v) > 1 && t.NodeMap(t.ChildSelect(v, 2)) != map[tree.ChildSelect(i, 2)]) 
            || t.ChildRank(v) != tree.ChildRank(i)
            || t.NodeMap(t.LCA(v, t.NodeSelect(internalN))) != map[tree.LCA(i, internalN)]
            )
          ++mismatchCnt ;
        //printf("%d: %d %d\n", v, t.NodeMap(t.Parent(v)), tree.Parent(i)) ;
      }
      printf("\nLOUDS cardinal tree mismatch count %u\n", mismatchCnt) ;
      printf("LOUDS cardinal tree space usage (bytes): %lu\n", t.GetSpace()) ;
    }
    
    {
      Tree_Cardinal_Ordinal<> t ;
      mismatchCnt = 0 ;

      t.Init(tree.GetTreeData().data(), tree.GetSize(), c, map) ;
      for (i = 0 ; i < internalN ; ++i)
      {
        size_t v = t.NodeSelect(map[i]) ;
        if (t.ChildrenCount(v) != tree.ChildrenCount(i) 
            || t.NodeMap(t.FirstChild(v)) != map[tree.FirstChild(i)] 
            || t.NodeMap(t.LastChild(v)) != map[tree.LastChild(i)]
            || t.NodeMap(t.Parent(v)) != map[tree.Parent(i)]
            || (t.ChildrenCount(v) > 1 && t.NodeMap(t.ChildSelect(v, 2)) != map[tree.ChildSelect(i, 2)]) 
            || t.ChildRank(v) != tree.ChildRank(i)
            || t.NodeMap(t.LCA(v, t.NodeSelect(map[internalN]))) != map[tree.LCA(i, internalN)]
           )
          ++mismatchCnt ;
        //printf("%d %d %d: %d %d\n", v, i, map[i], t.NodeMap(t.LCA(v, t.NodeSelect(map[internalN]))), map[tree.LCA(i, internalN)]) ;
      }
      printf("\nDFUDS cardinal tree mismatch count %u\n", mismatchCnt) ;
      printf("DFUDS cardinal tree space usage (bytes): %lu\n", t.GetSpace()) ;
    }
    delete[] map ;
  }

  PrintLog("Done") ;
  return 0 ;
}

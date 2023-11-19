#ifndef _MOURISL_COMPACTDS_FM_BUILDER
#define _MOURISL_COMPACTDS_FM_BUILDER

// Build BWT and other auxiliary datas from text T using blockwise suffix array sorting

#include <utility>
#include <map>

#include <pthread.h> 

#include "Utils.hpp"
#include "SuffixArrayGenerator.hpp"

namespace compactds {
struct _FMBuilderParam
{
  size_t n ;

  size_t saBlockSize ;
  int saDcv ;  
  size_t threadCnt ;
  
  int sampleRate ;
  int sampleStrategy ; // on SA, on T or on the ends of BWT runs.
  size_t sampleSize ;
  size_t *sampledSA ;

  int precomputeWidth ;
  size_t precomputeSize ;
  std::pair<size_t, size_t> *precomputedRange ; 
  
  bool printLog ;

  size_t maxLcp ; // only consider LCP up to this point

  std::map<size_t, size_t> selectedISA ;  
  std::map<size_t, size_t> selectedSA ; // reverse selectedISA

  WORD *semiLcpGreater ; // The LCP is between current suffix and its previous one
  WORD *semiLcpEqual ;

  size_t adjustedSA0 ; // specialized sampled SA.

  FILE *dumpSaFp ; // dump SA to this file.

  _FMBuilderParam()
  {
    sampleStrategy = 0 ;
    saBlockSize = 1<<24 ;
    saDcv = 4096 ;
    sampleRate = 1<<5 ;
    threadCnt = 1 ; // the number of threads for sorting. 
    precomputeWidth = 10 ;
    adjustedSA0 = 0 ;

    printLog = true ;

    maxLcp = 0 ;
    dumpSaFp = NULL ; 
    
    // The memory for these arrays shall handled explicitly outside.
    sampledSA = NULL ;
    precomputedRange = NULL ;
    semiLcpGreater = NULL ;
    semiLcpEqual = NULL ;
  }

  // Use this free with caution,
  //   as some of the pointer will be 
  //   used in FMIndexAuxData.
  // Use this only when generating BWT string.
  void Free() 
  {
    if (sampledSA != NULL)
      free(sampledSA) ;
    if (precomputedRange != NULL)
      free(precomputedRange) ;
    if (semiLcpGreater != NULL)
      free(semiLcpGreater) ;
    if (semiLcpEqual != NULL)
      free(semiLcpEqual) ;
  }
} ;

struct _FMBuilderChunkThreadArg
{
  int tid ;
  int threadCnt ;
  
  FixedSizeElemArray *T ;
  size_t n ;
  
  SuffixArrayGenerator *saGenerator ;
  size_t from, to ;
  std::vector< std::vector<size_t> > pos ;
} ;

struct _FMBuilderSASortThreadArg 
{
  int tid ;
  int threadCnt ;

  FixedSizeElemArray *T ;
  size_t n ;

  SuffixArrayGenerator *saGenerator ;
  size_t *sa ;
  size_t saSize ;
  
  size_t accuChunkSize ;

  WORD *semiLcpGreater ; // The LCP is between current suffix and its previous one
  WORD *semiLcpEqual ;
  size_t maxLcp ; // only consider LCP up to this point
} ;

class FMBuilder
{
private:
  // Return the LCP up until the specified value bewteen T[i,...], T[j,...]
  static size_t ComputeSemiLcp(FixedSizeElemArray &T, size_t n, size_t i, size_t j, size_t maxLCP)
  {
    size_t k ;
    if (i >= n || j >= n || i < 0 || j < 0)
      return 0 ;
    else
    {
      for (k = 0 ; k < maxLCP && i + k < n && j + k < n ; ++k)
      {
        if (T.Read(i + k) != T.Read(j + k))
          break ;
      }
      return k ;
    }
  }

  static void *PosInChunk_Thread(void *arg)
  {
    struct _FMBuilderChunkThreadArg *pArg = (_FMBuilderChunkThreadArg *)arg ;
    size_t segLen = DIV_CEIL(pArg->n, pArg->threadCnt) ;
    size_t s = segLen * pArg->tid ;
    size_t e = s + segLen - 1 ;
    pArg->saGenerator->GetChunksPositions(*(pArg->T), pArg->n, 
        pArg->from, pArg->to, s, e, pArg->pos) ;
    pthread_exit(NULL) ;
  }
  
  // Compare the semiLCP between T[sai...], and T[saj,...], write the result to semiLcp[biti]
  static void SetSemiLcpBit(FixedSizeElemArray &T, size_t n, size_t sai, size_t saj, size_t biti, size_t maxLcp, WORD *semiLcpGreater, WORD *semiLcpEqual)
  {
    size_t l = 0 ;
    l = ComputeSemiLcp(T, n, sai, saj, maxLcp + 1) ;
    if (l > maxLcp)
      Utils::BitSet(semiLcpGreater, biti) ;
    else if (l == maxLcp)
      Utils::BitSet(semiLcpEqual, biti) ;
  }

  static void *SortSA_Thread(void *arg)
  {
    struct _FMBuilderSASortThreadArg *pArg = (struct _FMBuilderSASortThreadArg *)arg ;
    pArg->saGenerator->SortSuffixByPos(*(pArg->T),pArg->n, 
        pArg->sa, pArg->saSize, pArg->sa ) ;
    //printf("TEST %d\n",  saSortThreadArgs[0][0].sa[0]) ;
    
    if (pArg->maxLcp > 0)
    {
      size_t i ;
      // The first element's LCP is between the last element from previous
      // chunk, need to process outside.
      for (i = 1 ; i < pArg->saSize ; ++i)
        SetSemiLcpBit(*(pArg->T), pArg->n, pArg->sa[i], pArg->sa[i - 1], pArg->accuChunkSize + i, 
            pArg->maxLcp, pArg->semiLcpGreater, pArg->semiLcpEqual) ;
    }
    
    pthread_exit(NULL) ;
  }

public:
  // Allocate and init the memorys for auxiliary data arrays in FM index
  // chrbit: number of bits for each character
  static void MallocAuxiliaryData(size_t chrbit, size_t n, struct _FMBuilderParam &param)
  {
    size_t i ;
    param.n = n ;
    
    param.sampleSize = DIV_CEIL(n, param.sampleRate) ;
    param.sampledSA = (size_t *)malloc(sizeof(size_t) * DIV_CEIL(n, param.sampleRate)) ;

    size_t size = 1ull<<(chrbit * param.precomputeWidth) ;
    param.precomputeSize = size ;  
    param.precomputedRange = (std::pair<size_t, size_t> *)malloc(
        sizeof(std::pair<size_t, size_t>) * size) ;
    for (i = 0 ; i < size ; ++i)
    {
      param.precomputedRange[i].first = 0 ;
      param.precomputedRange[i].second = 0 ;
    }

    if (param.maxLcp > 0)
    {
      param.semiLcpGreater = Utils::MallocByBits(n) ;
      param.semiLcpEqual = Utils::MallocByBits(n) ;
    }
  }

  // Determine the parameters for block size and difference cover size
  //   based on memory requirement (bytes).
  // Assume mem is quite large.
  static void InferParametersGivenMemory(size_t n, int alphabetSize, size_t memory,
      struct _FMBuilderParam &param)
  {
    size_t logBlockSize ;
    size_t dcv ;
    size_t alphabetBits = Utils::Log2Ceil(alphabetSize) ; 

    size_t bestTime = POSITIVE_INF ;
    size_t bestBlockSize = 0 ;
    size_t bestDcv = 0 ;
    
    if (2 * n * alphabetBits / 8 > memory)
      return ;
    
    memory -= 2 * n * alphabetBits / 8 ;
    for (dcv = 512 ; dcv <= 8196 ; dcv *= 2)
    {
      size_t dcSize = DIV_CEIL(n, dcv) * DifferenceCover::EstimateCoverSize(dcv) ;
      for (logBlockSize = 24 ; logBlockSize <= 50 ; ++logBlockSize)
      {
        size_t blockSize = 1ull<<logBlockSize ;
        //if (blockSize >= n / param.threadCnt)
        //  break ;

        size_t space = (param.threadCnt * blockSize 
            + dcSize + DIV_CEIL(n, param.sampleRate)
            + (1ull<<(alphabetBits * param.precomputeWidth))*2
            ) * WORDBYTES ;  
        
        if (space <= memory)
        {
          size_t iterations = DIV_CEIL(n, (param.threadCnt * blockSize)) ;
          size_t time = dcSize * Utils::Log2Ceil(n) // sort difference cover 
            + iterations * n * dcv // making cuts
            + iterations * (blockSize * Utils::Log2Ceil(blockSize) + dcv * blockSize) ; // sort block
          //printf("%lu(%lu) %lu %lu. %lu %lu. %lu\n", dcv, dcSize, blockSize, iterations, 
          //    space, time, memory) ; 
          if (time < bestTime)
          {
            bestBlockSize = blockSize ;
            bestDcv = dcv ;
            bestTime = time ;
          }
        }
        else
          break ;
      }
    }

    if (bestDcv != 0)
    {
      param.saBlockSize = bestBlockSize ;
      param.saDcv = bestDcv ;

      if (param.printLog)
      {
        Utils::PrintLog("Estimated block size: %lu; dcv:%d",
            param.saBlockSize, param.saDcv) ;
      }
    }
  }

  // T: text
  // n: len(text)
  // firstISA: ISA[0]
  // Returned information is in BWT, firstISA, which are important in the F column. param holds all the other allocated array.
  static void Build(FixedSizeElemArray &T, size_t n, int alphabetSize, 
      FixedSizeElemArray &BWT, size_t &firstISA,
      struct _FMBuilderParam &param)
  {
    size_t i, j, k ;
    SuffixArrayGenerator saGenerator ;
    MallocAuxiliaryData(Utils::Log2Ceil(alphabetSize), n, param) ; 
    BWT.Malloc(Utils::Log2Ceil(alphabetSize), n) ;
    if (param.printLog)
      Utils::PrintLog("Generate difference cover and chunks.") ;
    size_t cutCnt = saGenerator.Init(T, n, param.saBlockSize, param.saDcv, alphabetSize) ;
    if (param.printLog)
      Utils::PrintLog("Found %llu chunks.", cutCnt) ;
    size_t bwtFilled = 0 ;
   
    pthread_t *threads = (pthread_t *)malloc(sizeof(*threads) * param.threadCnt) ;
    struct _FMBuilderChunkThreadArg *chunkThreadArgs ;
    struct _FMBuilderSASortThreadArg *saSortThreadArgs ; 
    pthread_attr_t attr ;

    pthread_attr_init( &attr ) ;
    pthread_attr_setdetachstate( &attr, PTHREAD_CREATE_JOINABLE ) ;
    size_t **sa ; // suffix array chunks
    size_t *saChunkSize ; // actual size
    size_t *saChunkCapacity ; // the memory capacity
    
    chunkThreadArgs = new struct _FMBuilderChunkThreadArg[param.threadCnt] ;
    for (i = 0 ; i < param.threadCnt ; ++i)
    {
      chunkThreadArgs[i].tid = i ;
      chunkThreadArgs[i].threadCnt = param.threadCnt ;
      chunkThreadArgs[i].saGenerator = &saGenerator ;
      chunkThreadArgs[i].T = &T ;
      chunkThreadArgs[i].n = n ;
    }

    sa = (size_t **)malloc(sizeof(sa[0]) * param.threadCnt) ;
    saChunkSize = (size_t *)malloc(sizeof(saChunkSize) * param.threadCnt) ;
    saChunkCapacity = (size_t *)malloc(sizeof(saChunkCapacity) * param.threadCnt) ;
    saSortThreadArgs = (struct _FMBuilderSASortThreadArg*)malloc(sizeof(struct _FMBuilderSASortThreadArg) * param.threadCnt) ; 
    for (i = 0 ; i < param.threadCnt ; ++i)
    {
      sa[i] = NULL ;
      saChunkSize[i] = 0 ;
      saChunkCapacity[i] = 0 ;

      saSortThreadArgs[i].tid = i ;
      saSortThreadArgs[i].threadCnt = param.threadCnt ;
      saSortThreadArgs[i].saGenerator = &saGenerator ;
      saSortThreadArgs[i].T = &T ;
      saSortThreadArgs[i].n = n ;

      saSortThreadArgs[i].maxLcp = param.maxLcp ;
      saSortThreadArgs[i].semiLcpGreater = param.semiLcpGreater ;
      saSortThreadArgs[i].semiLcpEqual = param.semiLcpEqual ;
    }

    size_t lastSA = 0 ; // record the last SA from previous batch or chunk
    size_t accuChunkSizeForSort = 0 ; // accumulated chunk size
    i = 0 ;
    
    if (param.dumpSaFp)
      fwrite(&n, sizeof(size_t), 1, param.dumpSaFp) ;

    // Start the core iterations
    for (i = 0 ; i < cutCnt ; i += param.threadCnt)
    {
      // Load positions for current batch
      if (param.printLog)
        Utils::PrintLog("Extract %d chunks. (%lu/%lu chunks finished)", param.threadCnt, i, cutCnt) ;
      for (j = 0 ; j < param.threadCnt ; ++j)
      {
        chunkThreadArgs[j].from = i ;
        chunkThreadArgs[j].to = (i + param.threadCnt - 1 < n ? i + param.threadCnt - 1 : n - 1) ;
        pthread_create(&threads[j], &attr, PosInChunk_Thread, (void *)(chunkThreadArgs + j)) ;
        //PosInChunk_Thread((void *)(chunkThreadArgs + j)) ;
      }

      if (param.printLog)
        Utils::PrintLog("Wait for the chunk extraction to finish.") ;
      for (j = 0 ; j < param.threadCnt ; ++j)
        pthread_join(threads[j], NULL) ;

      size_t chunkCnt = param.threadCnt ;
      if (i + chunkCnt >= cutCnt)
        chunkCnt = cutCnt - i ;

      // concatenate the pos in the chunks
      for (j = 0 ; j < chunkCnt ; ++j)
      {
        size_t totalSize = 0 ;
        for (k = 0 ; k < param.threadCnt ; ++k)
          totalSize += chunkThreadArgs[k].pos[j].size() ;
        saChunkSize[j] = totalSize ; 
        if (totalSize > saChunkCapacity[j])
        {
          free(sa[j]) ;
          saChunkCapacity[j] = totalSize ;
          sa[j] = (size_t *)malloc(sizeof(sa[j]) * totalSize) ;
        }

        totalSize = 0 ;
        for (k = 0 ; k < param.threadCnt ; ++k)
        {
          memcpy(sa[j] + totalSize, chunkThreadArgs[k].pos[j].data(),
              sizeof(sa[j][0]) * chunkThreadArgs[k].pos[j].size()) ;
          totalSize += chunkThreadArgs[k].pos[j].size() ;
          std::vector<size_t>().swap(chunkThreadArgs[k].pos[j]) ;
        }
      }

      // Submit the batch of chunks to sorting 
      if (param.printLog)
        Utils::PrintLog("Submit %d chunks.", chunkCnt) ;
      for (j = 0 ; j < chunkCnt ; ++j)
      {
        if (param.printLog)
          Utils::PrintLog("Chunk %d elements: %llu", j, saChunkSize[j]) ;
        saSortThreadArgs[j].sa = sa[j] ;
        saSortThreadArgs[j].saSize = saChunkSize[j] ;
        saSortThreadArgs[j].accuChunkSize = accuChunkSizeForSort ;
        accuChunkSizeForSort += saChunkSize[j] ;
        pthread_create(&threads[j], &attr, SortSA_Thread, (void *)(saSortThreadArgs + j)) ;
        //SortSA_Thread( (void *)(saSortThreadArgs + j)) ;
      }

      // Wait for current batch to finish
      if (param.printLog)
        Utils::PrintLog("Wait for the chunk sort to finish.") ;
      for (j = 0 ; j < chunkCnt ; ++j)
        pthread_join(threads[j], NULL) ;

      // Process the information from the chunks. 
      if (param.printLog)
        Utils::PrintLog("Postprocess %d chunks.", chunkCnt) ;
      for (j = 0 ; j < chunkCnt ; ++j) 
      {
        size_t l ;
        size_t size = saSortThreadArgs[j].saSize ;
        size_t *saChunk = saSortThreadArgs[j].sa ;

        // Fill FM string
        //printf("%d %d %d %d\n", size, j, saSortThreadArgs[prevPosTag][j].pos->at(1),
        //    saChunk[0]) ;
        for (l = 0 ; l < size ; ++l)
        {
          if (saChunk[l] == 0)
          {
            firstISA = bwtFilled ;
            BWT.Write(bwtFilled, T.Read(n - 1)) ;
          }
          else
            BWT.Write(bwtFilled, T.Read( saChunk[l] - 1 ) ) ;

          if (param.sampledSA != NULL && bwtFilled % param.sampleRate == 0)
            param.sampledSA[bwtFilled / param.sampleRate] = saChunk[l] ;

          if (param.precomputedRange != NULL)
          {
            int width = param.precomputeWidth ;
            WORD w = 0 ;// word
            if (saChunk[l] + width <= n)
            {
              w = T.PackRead(saChunk[l], width) ;
              if (param.precomputedRange[w].second == 0)
                param.precomputedRange[w].first = bwtFilled ;
              ++param.precomputedRange[w].second ;
            } 
            /*else // ignore the case near the end of the string
            {
              w = T.PackRead(saChunk[l], n - saChunk[l]) ;
              //size_t used = n - saChunk[l] ;
              //w = (T.PackRead(0, w - used)) << used | w
            }*/

          }

          if (param.selectedISA.size() != 0 )
          {
            if (param.selectedISA.find(saChunk[l]) != param.selectedISA.end()) 
              param.selectedISA[saChunk[l]] = bwtFilled ;
          }

          ++bwtFilled ;
        }

        if (param.maxLcp > 0)
        {
          size_t offseti = bwtFilled - size ; // equiavlent to accuChunkSize
          if (i > 0 || j > 0) // ignore the very first SA in the whole array
            SetSemiLcpBit(T, n, saChunk[0], lastSA, offseti, param.maxLcp, 
                param.semiLcpGreater, param.semiLcpEqual) ;                 
        }

        // the last element from previous chunk. 
        lastSA = saChunk[size - 1] ;
        
        if (param.dumpSaFp)
          fwrite(saChunk, sizeof(saChunk[0]), size, param.dumpSaFp) ;
      }
    } // end of the main while loop for populating BWTs
    
    // Fill in the selectedSA
    for (std::map<size_t, size_t>::iterator iter = param.selectedISA.begin() ;
        iter != param.selectedISA.end(); ++iter)
    {
      param.selectedSA[iter->second] = iter->first ;
    }
    std::map<size_t, size_t>().swap(param.selectedISA) ; // ISA will not be useful

    free(threads) ;
    pthread_attr_destroy(&attr) ;
    delete[] chunkThreadArgs ;
    for (j = 0 ; j < param.threadCnt ; ++j)
    {
      if (sa[j] != NULL)
      {
        free(sa[j]) ;
      }
    }
    free(sa) ;
    free(saChunkSize) ;
    free(saChunkCapacity) ;
    free(saSortThreadArgs) ;
  }
} ;
}

#endif

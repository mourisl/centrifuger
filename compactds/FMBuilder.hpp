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

  bool hasEndMarker ;
  size_t endMarkerCnt ; // the first alphabet is to mark the end of a sequence/genome.
                          // this record how many end markers (number of sequences) in the input
                          // 0 if endMarker is not used.
  size_t *endMarkerSA ;

  size_t adjustedSA0 ; // specialized sampled SA.

  FILE *dumpSaFp ; // dump SA to this file.

  bool hasCheckpointFile ; // Use checkpoint file, will resume the progress if the checkpoint exists 
  // The prefix for the checkpoint file
  //    .1: progress file. 1: difference cover finished; n+1: n-th chunk is finished
  //    .2: difference cover SAs
  //    .3: filled BWT. Since BWT is stored consecutively, we will update the remaining checkpoint files every 10% of chunks are processed.
  //        followed by firstISA, sampledSA, precomputedRange, selectedISA, endMarkerSA, lastSA, accuChunkSizeForSort
  std::string checkpointFilePrefix ; 

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
  
    hasEndMarker = false ;
    endMarkerCnt = 0 ;
    endMarkerSA = NULL ;

    hasCheckpointFile = false ;

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
    if (endMarkerSA != NULL)
      free(endMarkerSA) ;
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
} ;

struct _FMBuilderPostprocessThreadArg
{
  int tid ;
  int threadCnt ;

  FixedSizeElemArray *T ;
  FixedSizeElemArray *BWT ;
  size_t n ;

  size_t *saChunk ;
  size_t saSize ;
  size_t prevChunkLastSA ;
  size_t *pFirstISA ;

  size_t accuChunkSize ; // The start for this chunk

  int skippedBWT ; // The number of BWT entries that skipped becuase they overlap with the WORD from the previous chun/k

  WORD firstPrecomputeW ; // The first precompute w range might be the same as the last w in the previous chunk in parallel, so we store them here, and merge after the parallel execution 
  size_t firstPrecomputeWLen ; 

  struct _FMBuilderParam *builderParam ;
} ;

class FMBuilder
{
private:
  // Return the LCP up until the specified value bewteen T[i,...], T[j,...]
  static size_t ComputeSemiLcp(const FixedSizeElemArray &T, size_t n, size_t i, size_t j, size_t maxLCP)
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
  static void SetSemiLcpBit(const FixedSizeElemArray &T, size_t n, size_t sai, size_t saj, size_t biti, size_t maxLcp, WORD *semiLcpGreater, WORD *semiLcpEqual)
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
    pthread_exit(NULL) ;
  }

	// This is the thread handling filling BWT and other components
	// Some of the other componenets wil be filled after the parallel excution
	static void *Postprocess_Thread(void *arg)
  {
    struct _FMBuilderPostprocessThreadArg *pArg = ( struct _FMBuilderPostprocessThreadArg *)arg ;
    int tid = pArg->tid ;

    size_t i ;
    size_t size = pArg->saSize ;
    size_t *saChunk = pArg->saChunk ;
    struct _FMBuilderParam &param = *(pArg->builderParam) ;
    const FixedSizeElemArray &T = *(pArg->T) ;
    FixedSizeElemArray &BWT = *(pArg->BWT) ;
    size_t n = pArg->n ;

    // Fill FM string
    //printf("%d %d %d %d\n", size, j, saSortThreadArgs[prevPosTag][j].pos->at(1),
    //    saChunk[0]) ;
    size_t bwtFilled = pArg->accuChunkSize ;
    int skipLength = 0 ; // skip this amount of BWT as they may write to a word 
    if (tid > 0 && BWT.GetElemOffsetInWord(bwtFilled) > 0)
    {
      while (BWT.GetElemWordIndex(bwtFilled) == BWT.GetElemWordIndex(bwtFilled + skipLength))
      {
        ++skipLength ;
      }
    }
    pArg->skippedBWT = skipLength ;
    
    bool setFirstPrecomputeW = false ;
    for (i = 0 ; i < size ; ++i, ++bwtFilled)
    {
      if (i >= (size_t)skipLength)
      {
        if (saChunk[i] == 0)
        {
          *(pArg->pFirstISA) = bwtFilled ;
          BWT.Write(bwtFilled, T.Read(n - 1)) ;
        }
        else
          BWT.Write(bwtFilled, T.Read( saChunk[i] - 1 ) ) ;

        if (param.sampledSA != NULL && bwtFilled % param.sampleRate == 0)
          param.sampledSA[bwtFilled / param.sampleRate] = saChunk[i] ;
      }

      if (param.precomputedRange != NULL)
      {
        int width = param.precomputeWidth ;
        WORD w = 0 ;
        
        if (saChunk[i] + width <= n)
        {
          w = T.PackRead(saChunk[i], width) ;
          if (!setFirstPrecomputeW && tid > 0)
          {
            pArg->firstPrecomputeW = w ;
            pArg->firstPrecomputeWLen = 1 ;
            setFirstPrecomputeW = true ;
          }
          else
          {
            if (w == pArg->firstPrecomputeW && tid > 0)
            {
              ++pArg->firstPrecomputeWLen ;
            }
            else
            {
              if (param.precomputedRange[w].second == 0)
                param.precomputedRange[w].first = bwtFilled ;
              ++param.precomputedRange[w].second ;
            }
          }
        } 
        /*else // ignore the case near the end of the string
          {
          w = T.PackRead(saChunk[i], n - saChunk[i]) ;
        //size_t used = n - saChunk[i] ;
        //w = (T.PackRead(0, w - used)) << used | w
        }*/
      }

      // Since we are not inserting to the map structure, this should not affect the results
      if (param.selectedISA.size() != 0 )
      {
        if (param.selectedISA.find(saChunk[i]) != param.selectedISA.end()) 
          param.selectedISA[saChunk[i]] = bwtFilled ;
      }

      if (param.maxLcp > 0 && i > 0)
      {
        //TODO: need to address the racing for the semiLcpBits
        SetSemiLcpBit(T, n, saChunk[i], saChunk[i - 1], pArg->accuChunkSize + i, 
            param.maxLcp, param.semiLcpGreater, param.semiLcpEqual) ;
      }

      if (param.hasEndMarker && bwtFilled < param.endMarkerCnt)
      {
        // The endmarkers are always in the first of the SA, so we can directly 
        // use bwtFilled for endMarkerSA
        param.endMarkerSA[bwtFilled] = saChunk[i] ;
      }
    }

    if (param.maxLcp > 0 && pArg->accuChunkSize > 0) // ignore the very first SA in the whole array
      SetSemiLcpBit(T, n, saChunk[0], pArg->prevChunkLastSA, pArg->accuChunkSize, 
          param.maxLcp, param.semiLcpGreater, param.semiLcpEqual) ;                 
    pthread_exit(NULL) ;
  }

public:
  // Allocate and init the memorys for auxiliary data arrays in FM index
  // chrbit: number of bits for each character
  static void MallocAuxiliaryData(const FixedSizeElemArray &T, size_t chrbit, size_t n, struct _FMBuilderParam &param)
  {
    size_t i ;
    param.n = n ;
    
    param.sampleSize = DIV_CEIL(n, param.sampleRate) ;
    param.sampledSA = (size_t *)malloc(sizeof(size_t) * DIV_CEIL(n, param.sampleRate)) ;

    if (param.precomputeWidth > 0)
    {
      size_t size = 1ull<<(chrbit * param.precomputeWidth) ;
      param.precomputeSize = size ;  
      param.precomputedRange = (std::pair<size_t, size_t> *)malloc(
          sizeof(std::pair<size_t, size_t>) * size) ;
      for (i = 0 ; i < size ; ++i)
      {
        param.precomputedRange[i].first = 0 ;
        param.precomputedRange[i].second = 0 ;
      }
    }
    else
    {
      // Since we did not explicitly store the end of the text,
      //   we need an informative start position, otherwise the range
      //   may include the text end. (I think).
      Utils::PrintLog("precomputeWidth has to be greater than 0.\n") ;
      exit(1) ;
    }

    if (param.maxLcp > 0)
    {
      param.semiLcpGreater = Utils::MallocByBits(n) ;
      param.semiLcpEqual = Utils::MallocByBits(n) ;
    }

    if (param.hasEndMarker)
    {
      param.endMarkerCnt = 0 ;
      for (i = 0 ; i < n ; ++i)
        if (T[i] == 0)
          ++param.endMarkerCnt ;
      param.endMarkerSA = (size_t *)malloc(sizeof(param.endMarkerSA[0]) * param.endMarkerCnt) ;
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
		
    if (2 * n * alphabetBits / 8 > memory) // The input text and the output BWT
    {
      if (param.printLog)
        Utils::PrintLog("WARNING: memory is not enough for other block size and dcv values, will use the default.") ;
      return ;
    }
    
    memory -= 2 * n * alphabetBits / 8 ;
    for (dcv = 512 ; dcv <= 8196 ; dcv *= 2)
    {
      size_t dcSize = DIV_CEIL(n, dcv) * DifferenceCover::EstimateCoverSize(dcv) ;
      for (logBlockSize = 24 ; logBlockSize <= 50 ; ++logBlockSize)
      {
        size_t blockSize = 1ull<<logBlockSize ;
        //if (blockSize >= n / param.threadCnt)
        //  break ;
        size_t space = (2 * param.threadCnt * blockSize // SA position, SA result
            + dcSize // SA value for difference cover 
            + SuffixArrayGenerator::EstimateChunkCount(n, blockSize, dcv) * dcv // _cutLCP 
            + DIV_CEIL(n, param.sampleRate) // sampledSA
            + (1ull<<(alphabetBits * param.precomputeWidth))*2 // precompted width
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
    else if (param.printLog)
        Utils::PrintLog("WARNING: memory is not enough for other block size and dcv values, will use the default.") ;
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
   
    size_t checkpointStep = 0 ; // start from the beginning
    if (param.hasCheckpointFile)
    {
      FILE *fp = fopen((param.checkpointFilePrefix + ".1").c_str(), "r") ;
      if (fp != NULL)
      {
        fscanf(fp, "%lu", &checkpointStep) ;
        fclose(fp) ;
      }
    }

    size_t cutCnt = 0 ;
    if (checkpointStep < 1)
    {
      if (param.printLog)
        Utils::PrintLog("Generate difference cover and chunks.") ;
      cutCnt = saGenerator.Init(T, n, param.saBlockSize, param.saDcv, alphabetSize) ;
      if (param.printLog)
        Utils::PrintLog("Found %lu chunks.", cutCnt) ;
      
      if (param.hasCheckpointFile)
      {
        FILE *fp = fopen((param.checkpointFilePrefix + ".2").c_str(), "w") ;
        saGenerator.Save(fp) ;
        fclose(fp);

        fp = fopen((param.checkpointFilePrefix + ".1").c_str(), "w") ;
        fprintf(fp, "1") ;
        fclose(fp);
      }
    }
    else // checkpointStep >= 1
    {
      FILE *fp = fopen((param.checkpointFilePrefix + ".2").c_str(), "r") ;
      saGenerator.Load(fp) ;
      fclose(fp);

      cutCnt = saGenerator.GetChunkCount() ;
      if (param.printLog)
        Utils::PrintLog("Load %lu chunks from the checkpoint file.", cutCnt) ;
    }

    size_t alphabetBits = Utils::Log2Ceil(alphabetSize) ;
    MallocAuxiliaryData(T, alphabetBits, n, param) ; 
    BWT.Malloc(alphabetBits, n) ;
       
    pthread_t *threads = (pthread_t *)malloc(sizeof(*threads) * param.threadCnt) ;
    struct _FMBuilderChunkThreadArg *chunkThreadArgs ;
    struct _FMBuilderSASortThreadArg *saSortThreadArgs ; 
    struct _FMBuilderPostprocessThreadArg *postprocessThreadArgs ;
    pthread_attr_t attr ;

    pthread_attr_init( &attr ) ;
    pthread_attr_setdetachstate( &attr, PTHREAD_CREATE_JOINABLE ) ;
    size_t **sa ; // suffix array chunks
    size_t *saChunkSize ; // actual size
    size_t *saChunkCapacity ; // the memory capacity

    size_t pthreadStackSize, estimatedStackSize ;
    pthread_attr_getstacksize(&attr, &pthreadStackSize) ;
    estimatedStackSize = param.saDcv * 500 ;
    if (estimatedStackSize > pthreadStackSize) 
    {
      if (!pthread_attr_setstacksize(&attr, estimatedStackSize))
      {
        Utils::PrintLog("Increased the pthread stack size from %lu to %lu.", pthreadStackSize, estimatedStackSize) ;
      }
      else
      {
        Utils::PrintLog("WARNING: Could not set the pthread stack size with estimated stack size %lu. Please try a different difference cover and block size parameters if the program crashes later.", estimatedStackSize) ;
      }
    }

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
    }

    postprocessThreadArgs = (struct _FMBuilderPostprocessThreadArg*)malloc(sizeof(*postprocessThreadArgs) * param.threadCnt) ;
    for (i = 0 ; i < param.threadCnt ; ++i)
    {
      postprocessThreadArgs[i].tid = i ;
      postprocessThreadArgs[i].threadCnt = param.threadCnt ;
      postprocessThreadArgs[i].T = &T ;
      postprocessThreadArgs[i].BWT = &BWT ;
      postprocessThreadArgs[i].n = n ;
      postprocessThreadArgs[i].pFirstISA = &firstISA ;
      postprocessThreadArgs[i].builderParam = &param ;
    }

    size_t lastSA = 0 ; // record the last SA from previous batch or chunk
    size_t accuChunkSizeForSort = 0 ; // accumulated chunk size
    i = 0 ;
    
    if (param.dumpSaFp)
      fwrite(&n, sizeof(size_t), 1, param.dumpSaFp) ;

    // Load from the check point
    size_t chunkStart = 0 ;
    if (checkpointStep > 1)
    {
      chunkStart = (checkpointStep - 1) ;
      
      // .3: filled BWT. firstISA. 
      //  followed by sampledSA, precomputedRange, selectedISA, endMarkerSA, lastSA, accuChunkSizeForSort
      FILE *fp = fopen((param.checkpointFilePrefix + ".3").c_str(), "r") ;
      BWT.Load(fp) ;
      LOAD_VAR(fp, firstISA) ;
      
      LOAD_ARR(fp, param.sampledSA, param.sampleSize) ;
      LOAD_ARR(fp, param.precomputedRange, param.precomputeSize) ;
      size_t size ;
      LOAD_VAR(fp, size) ;
      for (i = 0 ; i < size ; ++i)
      {
        size_t a, b ;
        LOAD_VAR(fp, a) ;
        LOAD_VAR(fp, b) ;
        param.selectedISA[a] = b ;
      }
      if (param.hasEndMarker)
        LOAD_ARR(fp, param.endMarkerSA, param.endMarkerCnt) ;
      LOAD_VAR(fp, lastSA) ;
      LOAD_VAR(fp, accuChunkSizeForSort) ;  

      fclose(fp) ;
      if (param.printLog)
        Utils::PrintLog("Load processed data of the first %lu chunks from the checkpoint file.", chunkStart) ;
    }

    // Start the core iterations
    for (i = chunkStart ; i < cutCnt ; i += param.threadCnt)
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
      {
        pthread_join(threads[j], NULL) ;
        
        if (param.dumpSaFp)
          fwrite(saSortThreadArgs[j].sa, sizeof(saSortThreadArgs[j].sa[0]), 
              saSortThreadArgs[j].saSize, param.dumpSaFp) ;
      }

      // Process the information from the chunks. 
      if (param.printLog)
        Utils::PrintLog("Postprocess %d chunks.", chunkCnt) ;
      
      for (j = 0 ; j < chunkCnt ; ++j)
      {
        postprocessThreadArgs[j].saChunk = sa[j] ;
        postprocessThreadArgs[j].saSize = saChunkSize[j] ;
        postprocessThreadArgs[j].accuChunkSize = saSortThreadArgs[j].accuChunkSize ;
        
        // Variables that needed to simplify the overlap between previous chunk and current chunk.
        postprocessThreadArgs[j].skippedBWT = 0 ;
        postprocessThreadArgs[j].firstPrecomputeW = 0 ;
        postprocessThreadArgs[j].firstPrecomputeWLen = 0 ;
        
        // the last element from previous chunk. 
        postprocessThreadArgs[j].prevChunkLastSA = lastSA ;
        lastSA = sa[j][ saChunkSize[j] - 1 ] ;
      
        pthread_create(&threads[j], &attr, Postprocess_Thread, (void *)(postprocessThreadArgs + j)) ;
      }

      for (j = 0 ; j < chunkCnt ; ++j)
        pthread_join(threads[j], NULL) ;
      
      // Fill in the overlapped portion between SA chunks in BWT
      // Start from the second chunk.
      for (j = 1 ; j < chunkCnt ; ++j) 
      {
        int l ; 
        //size_t size = postprocessThreadArgs[j].saSize ;
        size_t *saChunk = postprocessThreadArgs[j].saChunk ;
        size_t accuChunkSize = postprocessThreadArgs[j].accuChunkSize ;

        for (l = 0 ; l < postprocessThreadArgs[j].skippedBWT ; ++l)
        {
          // Fill in the BWT
          size_t bwtFilled = accuChunkSize + l ;
          if (saChunk[l] == 0)
          {
            firstISA = bwtFilled ;
            BWT.Write(bwtFilled, T.Read(n - 1)) ;
          }
          else
            BWT.Write(bwtFilled, T.Read( saChunk[l] - 1 ) ) ;

          if (param.sampledSA != NULL && bwtFilled % param.sampleRate == 0)
            param.sampledSA[bwtFilled / param.sampleRate] = saChunk[l] ;
        }

        // Fill the precomputew
        if (param.precomputedRange != NULL)   
        {
          WORD w = postprocessThreadArgs[j].firstPrecomputeW ;
          size_t wlen = postprocessThreadArgs[j].firstPrecomputeWLen ;
          if (param.precomputedRange[w].second > 0)
          {
            param.precomputedRange[w].second += wlen ;
          }
          else // This w is the first one to show up.
          {
            // This also handles that the same precompute w spans more than one chunk.
            param.precomputedRange[w].first = accuChunkSize ;
            param.precomputedRange[w].second = wlen ;
          }
        }

        // TODO: Fill the lcp structure
      }

      // i+param.threadCnt is the number of finished blocks at this point
      if (param.hasCheckpointFile
          && i > 0 && 
          (i + param.threadCnt)/ (cutCnt / 10 + 1) > i / (cutCnt / 10 + 1))
      {
        // Reset the main checkpoint tracker to 1 just in case program crashes during the output
        FILE *fp = fopen((param.checkpointFilePrefix + ".1").c_str(), "w") ;
        fprintf(fp, "1") ;
        fclose(fp);

        // .3: filled BWT. firstISA 
        //  followed by sampledSA, precomputedRange, selectedISA, endMarkerSA, firstISA, lastSA, accuChunkSizeForSort
        fp = fopen((param.checkpointFilePrefix + ".3").c_str(), "w") ;
        BWT.Save(fp) ;
        SAVE_VAR(fp, firstISA) ;
        
        SAVE_ARR(fp, param.sampledSA, param.sampleSize) ;
        SAVE_ARR(fp, param.precomputedRange, param.precomputeSize) ;
        size_t size = param.selectedISA.size() ;
        SAVE_VAR(fp, size) ;
        for (std::map<size_t, size_t>::iterator iter = param.selectedISA.begin() ;
            iter != param.selectedISA.end() ; ++iter)
        {
          size_t a, b ;
          a = iter->first ;
          b = iter->second ;
          SAVE_VAR(fp, a) ;
          SAVE_VAR(fp, b) ;
        }
        if (param.hasEndMarker > 0)
          SAVE_ARR(fp, param.endMarkerSA, param.endMarkerCnt) ;
        SAVE_VAR(fp, lastSA) ;
        SAVE_VAR(fp, accuChunkSizeForSort) ;  
        fclose(fp) ;
        
        fp = fopen((param.checkpointFilePrefix + ".1").c_str(), "w") ;
        fprintf(fp, "%lu", 1 + (i + param.threadCnt > cutCnt ? cutCnt : i + param.threadCnt) ) ; // 1 + number_of_finished_chunks
        fclose(fp);
      } 
    } // end of the main while loop for populating BWTs
    
    // Fill in the selectedSA from selectedISA.
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
    free(postprocessThreadArgs) ;
  }
} ;
}

#endif

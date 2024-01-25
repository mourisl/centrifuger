#include <stdio.h>
#include <time.h>
#include <getopt.h>

#include "argvdefs.h"
#include "ReadFiles.hpp"
#include "compactds/Sequence_Hybrid.hpp"
#include "compactds/FMBuilder.hpp"
#include "compactds/FMIndex.hpp"
#include "compactds/Alphabet.hpp"

#include "Taxonomy.hpp"
#include "Classifier.hpp"
#include "ResultWriter.hpp"

//#define CENTRIFUGER_VERSION "Centrifuger v1.0.0"

char usage[] = "./centrifuger [OPTIONS]:\n"
  "Required:\n"
  "\t-x FILE: index prefix\n"
  "\t-1 FILE -2 FILE: paired-end read\n"
  "\t-u FILE: single-end read\n"
  //"\t--sample-sheet FILE: \n"
  "Optional:\n"
  //"\t-o STRING: output prefix [centrifuger]\n"
  "\t-t INT: number of threads [1]\n"
  "\t-k INT: report upto <int> distinct, primary assignments for each read pair [1]\n"
  "\t-v: print the version information and quit\n"
  "\t--min-hitlen INT: minimum length of partial hits [auto]\n"
  "\t--hitk-factor INT: resolve at most <int>*k entries for each hit [40; use 0 for no restriction]\n"
  ;

static const char *short_options = "x:1:2:u:o:t:k:v" ;
static struct option long_options[] = {
  { "min-hitlen", required_argument, 0, ARGV_MIN_HITLEN},
  { "hitk-factor", required_argument, 0, ARGV_MAX_RESULT_PER_HIT_FACTOR},
  { (char *)0, 0, 0, 0} 
} ;

struct _inputThreadArg
{
  ReadFiles *reads, *mateReads ;
  struct _Read *readBatch, *readBatch2 ;
  int maxBatchSize ;
  int *pBatchSize ;
} ;

struct _threadArg 
{
	struct _Read *readBatch, *readBatch2 ;
  
	int threadCnt ;
	int batchSize ;

	Classifier *classifier ;
  struct _classifierResult *results ;

	int tid ;
} ;

int GetReadBatch(ReadFiles &reads, struct _Read *readBatch, 
    ReadFiles &mateReads, struct _Read *readBatch2, int maxBatchSize)
{
  int fileInd1, fileInd2 ;
  
  int batchSize = reads.GetBatch( readBatch, maxBatchSize, fileInd1, true, true ) ;
  if ( readBatch2 != NULL )
  {
    int tmp = mateReads.GetBatch( readBatch2, maxBatchSize, fileInd2, true, true ) ;
    if ( tmp != batchSize )
    {
      Utils::PrintLog("ERROR: The two mate-pair read files have different number of reads." ) ;
      exit(EXIT_FAILURE) ;
    }
  }

  return batchSize ;
}

void *LoadReads_Thread(void *pArg)
{
  struct _inputThreadArg &arg = *((struct _inputThreadArg *)pArg);
  *(arg.pBatchSize) = GetReadBatch(*(arg.reads), arg.readBatch, 
      *(arg.mateReads), arg.readBatch2, arg.maxBatchSize) ;

  pthread_exit(NULL) ;
}

void *ClassifyReads_Thread(void *pArg)
{
  int i ;
  struct _threadArg &arg = *((struct _threadArg *)pArg);
  for (i = 0 ; i < arg.batchSize ; ++i)
  {
    if (i % arg.threadCnt != arg.tid)
      continue ;
    arg.classifier->Query(arg.readBatch[i].seq, arg.readBatch2 ? arg.readBatch2[i].seq : NULL, arg.results[i]) ;
  }
  pthread_exit(NULL) ;
}

int main(int argc, char *argv[])
{
  int i ;

	if ( argc <= 1 )
	{
		fprintf( stderr, "%s", usage ) ;
		return 0 ;
  }
	
  int c, option_index ;
	option_index = 0 ;
  
  char outputPrefix[1024] = "centrifuger" ;
  char *idxPrefix = NULL ;
  int threadCnt = 1 ;
  Classifier classifier ;
  struct _classifierParam classifierParam ;
  ReadFiles reads ;
  ReadFiles mateReads ;
  bool hasMate = false ;
  ResultWriter resWriter ;

  while (1)
  {
		c = getopt_long( argc, argv, short_options, long_options, &option_index ) ;
		
		if (c == -1)
			break ;
  
    if (c == 'x') // reference genome file
    {
      idxPrefix = strdup(optarg) ;
    }
    else if (c == 'u')
    {
      reads.AddReadFile(optarg, false) ;
    }
    else if (c == '1')
    {
      reads.AddReadFile(optarg, true) ;
      hasMate = true ;
    }
    else if (c == '2')
    {
			mateReads.AddReadFile( optarg, true ) ;
			hasMate = true ;
    }
    else if (c == 'o')
    {
      strcpy(outputPrefix, optarg) ;
    }
    else if (c == 't')
    {
      threadCnt = atoi(optarg) ;
    }
    else if (c == 'k')
    {
      classifierParam.maxResult = atoi(optarg) ;
    }
    else if (c == 'v')
    {
      printf("Centrifuger v" CENTRIFUGER_VERSION "\n") ;
      exit(0) ; 
    }
    else if (c == ARGV_MIN_HITLEN)
    {
      classifierParam.minHitLen = atoi(optarg) ;
    }
    else if (c == ARGV_MAX_RESULT_PER_HIT_FACTOR)
    {
      classifierParam.maxResultPerHitFactor = atoi(optarg) ;
    }
		else
		{
      Utils::PrintLog("Unknown parameter found.\n%s", usage ) ;
			return EXIT_FAILURE ;
		}
  }

  Utils::PrintLog("Centrifuger starts." ) ;
  if (idxPrefix == NULL)
  {
    Utils::PrintLog("Need to use -x to specify index prefix.") ;
    return EXIT_FAILURE ;
  }

  classifier.Init(idxPrefix, classifierParam) ;
  resWriter.OutputHeader() ;

  const int maxBatchSize = 1024 * threadCnt ;
  int batchSize ;
    
  
  bool useLoadOutputThreads = false ;

  int classificationThreadCnt = threadCnt ;
  if (threadCnt > 12)
  {
    classificationThreadCnt = threadCnt - 2 ;
    useLoadOutputThreads = true ;
  }

  pthread_t *threads = (pthread_t *)malloc( sizeof( pthread_t ) * classificationThreadCnt ) ;
  struct _threadArg *args = (struct _threadArg *)malloc( sizeof( struct _threadArg ) * classificationThreadCnt ) ;
  pthread_attr_t attr ;
  pthread_attr_init( &attr ) ;
  pthread_attr_setdetachstate( &attr, PTHREAD_CREATE_JOINABLE ) ;
  
  //useLoadOutputThreads = false ;
  if (!useLoadOutputThreads)
  {
    struct _Read *readBatch = NULL, *readBatch2 = NULL ;
    readBatch = ( struct _Read *)calloc( sizeof( struct _Read ), maxBatchSize ) ;
    if ( hasMate )
      readBatch2 = ( struct _Read *)calloc( sizeof( struct _Read ), maxBatchSize ) ;
    struct _classifierResult *classifierBatchResults = new struct _classifierResult[maxBatchSize] ;
    
    for ( i = 0 ; i < classificationThreadCnt ; ++i )
    {
      args[i].threadCnt = classificationThreadCnt ;
      args[i].tid = i ;
      args[i].readBatch = readBatch ;
      args[i].readBatch2 = readBatch2 ;
      args[i].results = classifierBatchResults ;
      args[i].classifier = &classifier ;
    }
    
    while ( 1 )
    {
      batchSize = GetReadBatch(reads, readBatch, mateReads, readBatch2, maxBatchSize) ;

      if ( batchSize == 0 )
        break ; 

      for ( i = 0 ; i < classificationThreadCnt ; ++i )
      {
        args[i].batchSize = batchSize ;
        pthread_create( &threads[i], &attr, ClassifyReads_Thread, (void *)&args[i] ) ;	
      }

      for ( i = 0 ; i < classificationThreadCnt ; ++i )
        pthread_join( threads[i], NULL ) ;

      for (i = 0 ; i < batchSize ; ++i)
        resWriter.Output(readBatch[i].id, 
            NULL, classifierBatchResults[i]) ;
    }
    
    reads.FreeBatch(readBatch, maxBatchSize) ;
    if (hasMate)
      mateReads.FreeBatch(readBatch2, maxBatchSize) ;
    free(readBatch) ;
    if (hasMate)
      free(readBatch2) ;
    delete[] classifierBatchResults ;
  }
  else
  {
    pthread_t inputThread ;
    struct _inputThreadArg inputThreadArg ;

    struct _Read *readBatch[3] ;
    struct _Read *readBatch2[3] ;
    struct _classifierResult *classifierBatchResults[3] ;
    
    for (i = 0 ; i < 3 ; ++i)
    {
      readBatch[i] = ( struct _Read *)calloc( sizeof( struct _Read ), maxBatchSize ) ;
      if ( hasMate )
        readBatch2[i] = ( struct _Read *)calloc( sizeof( struct _Read ), maxBatchSize ) ;
      else
        readBatch2[i] = NULL ;
      classifierBatchResults[i] = new struct _classifierResult[maxBatchSize] ;
    }
    int batchSize[3] ;
    
    bool started = false ;
    // Load in the first batch
    batchSize[0] = GetReadBatch(reads, readBatch[0], mateReads, readBatch2[0], maxBatchSize) ;
    
    int tag = 0 ; // which batch to use
    inputThreadArg.reads = &reads ;
    inputThreadArg.mateReads = &mateReads ;
    inputThreadArg.maxBatchSize = maxBatchSize ;
    for ( i = 0 ; i < classificationThreadCnt ; ++i )
    {
      args[i].threadCnt = classificationThreadCnt ;
      args[i].tid = i ;
      args[i].classifier = &classifier ;
    }

    while (1)
    {
      int nextTag = (tag + 1) % 3 ;
      int prevTag = tag >= 1 ? tag - 1 : 2 ;
        
      if (started)
        pthread_join(inputThread, NULL) ;
      if (batchSize[tag] > 0)
      {
        // Load in the next batch
        inputThreadArg.readBatch = readBatch[nextTag] ;
        inputThreadArg.readBatch2 = readBatch2[nextTag] ;
        inputThreadArg.pBatchSize = &batchSize[nextTag] ;
        pthread_create(&inputThread, &attr, LoadReads_Thread, (void *)&inputThreadArg) ;
      }

      // Process the current batch
      if (started)
      {
        for (i = 0 ; i < classificationThreadCnt ; ++i)
          pthread_join(threads[i], NULL) ;
      }
      
      if (batchSize[tag] > 0)
      {
        for ( i = 0 ; i < classificationThreadCnt ; ++i )
        {
          args[i].readBatch = readBatch[tag] ;
          args[i].readBatch2 = readBatch2[tag] ;
          args[i].results = classifierBatchResults[tag] ;
          args[i].batchSize = batchSize[tag] ;

          pthread_create( &threads[i], &attr, ClassifyReads_Thread, (void *)&args[i] ) ;	
        }
      }

      // Output the previous batch
      if (started)
      {
        for (i = 0 ; i < batchSize[prevTag] ; ++i)
          resWriter.Output(readBatch[prevTag][i].id, 
              NULL, classifierBatchResults[prevTag][i]) ;
      }
      
      if (batchSize[tag] == 0)
        break ;
      
      started = true ;
      tag = nextTag ;
    }

    for (i = 0 ; i < 3 ; ++i)
    {
      reads.FreeBatch(readBatch[i], maxBatchSize) ;
      free(readBatch[i]) ;
      if (hasMate)
      {
        mateReads.FreeBatch(readBatch2[i], maxBatchSize) ;
        free(readBatch2[i]) ;
      }
      delete[] classifierBatchResults[i] ;
    }
  } // end of if-else for use input output thread
  
  pthread_attr_destroy( &attr ) ;
  free( threads ) ;
  free( args ) ;
  free(idxPrefix) ;

  Utils::PrintLog("Centrifuger finishes." ) ;
  return 0 ;
}

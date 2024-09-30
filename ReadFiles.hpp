// The class for reading reads from a file
#ifndef _MOURISL_READS
#define _MOURISL_READS

#include <stdio.h>
#include <string.h>
#include <zlib.h>

#include <glob.h>

#include <vector>
#include <string>

#include "defs.h"
#include "kseq.h"

KSEQ_INIT( gzFile, gzread ) ;

struct _Read
{
  char *id ;
  char *seq ;
  char *qual ;
} ;

class ReadFiles
{
  private:
    std::vector<std::string> fileNames ;
    std::vector<bool> hasMate ;

    gzFile gzFp ;
    kseq_t *inSeq ;
    int fileType ; // 0-FASTA, 1-FASTQ
    int fileCnt ;
    int currentFpInd ;
    
    bool opened ;

    void GetFileBaseName(const char *in, char *out ) 
    {
      int i, j, k ;
      int len = (int)strlen( in ) ;
      for ( i = len ; i >= 0 && in[i] != '.' && in[i] != '/' ; --i )
        ;
      for ( j = len ; j >= 0 && in[j] != '/' ; --j )
        ;
      if ( i >= 0 && in[i] == '.' )
      {
        for (k = j + 1 ; k < i ; ++k)
          out[k - (j + 1)] = in[k] ;
        out[k - (j + 1)] = '\0' ;
      }
      else
      {
        strcpy( out, in + j + 1 ) ;
      }
    }

    void OpenFile(int fileInd)
    {
      if (opened)
      {
        kseq_destroy(inSeq) ;
        gzclose( gzFp ) ;
      }

      opened = true ;
      gzFp = gzopen( fileNames[fileInd].c_str(), "r" ) ;
      inSeq = kseq_init( gzFp ) ;

      kseq_read( inSeq ) ;
      if ( inSeq->qual.l == 0 )
      {
        fileType = 0 ;
        //qual[0] = '\0' ;
      }
      else 
      {
        fileType = 1 ;
      }
      /*else
      {
        fprintf( stderr, "\"%s\"'s format is wrong.\n", file ) ;
        exit( 1 ) ;
      }*/

      //printf( "%s %s\n", inSeq[fileCnt]->name.s, inSeq[fileCnt]->comment.s ) ;
      gzrewind( gzFp ) ;
      kseq_rewind( inSeq ) ;
      //kseq_read( inSeq[ fileCnt ] ) ;
      //printf( "%s %s\n", inSeq[fileCnt]->name.s, inSeq[fileCnt]->comment.s ) ;
      //gzrewind( gzFp[ fileCnt ]) ;
      //kseq_rewind( inSeq[ fileCnt] ) ;
    }

    void RemoveReadIdSuffix(char *id)
    {
      int len = strlen( id ) ;
      if ( ( id[len - 1] == '1' || id[len - 1] == '2' )
          && id[len - 2] == '/' )
      {
        id[len - 2] = '\0' ;
      }
    }
  public:
    char *id ;
    char *seq ;
    char *qual ;

    ReadFiles(): fileCnt(0), currentFpInd(0), opened(false)
    {
      id = seq = qual = NULL ;
    }

    ~ReadFiles()
    {
      Free() ;
    }

    void Free()
    {
      if ( id != NULL )
        free( id ) ;
      if ( seq != NULL )
        free( seq ) ;
      if ( qual != NULL )
        free( qual ) ;

      if (opened)
      {
        kseq_destroy( inSeq) ;
        gzclose( gzFp ) ;
      
        opened = false ;
      }
    }

    void AddReadFile( char *file, bool fileHasMate )
    {
      //std::string s(file) ;
      unsigned int i ;
      bool needGlob = false ;
      for (i = 0 ; file[i] ; ++i)
        if (file[i] == '*')
          needGlob = true ;
      
      int addFileCnt = 1 ;
      if (!needGlob)
      {
        fileNames.push_back(file) ;
        hasMate.push_back(fileHasMate) ;
      }
      else
      {
        glob_t globResult ;
        memset(&globResult, 0, sizeof(globResult)) ;
        
        int globRet = glob(file, GLOB_TILDE, NULL, &globResult) ;
        if (globRet != 0)
        {
          globfree(&globResult) ;
          fprintf(stderr, "glob() failed with return value %d.\n", globRet) ;
        }

        for (i = 0 ; i < globResult.gl_pathc ; ++i)
        {
          fileNames.push_back(globResult.gl_pathv[i]) ;
          hasMate.push_back(fileHasMate) ;
        }

        addFileCnt = globResult.gl_pathc ;
        globfree(&globResult) ;
      }

      if (fileCnt == 0)
        OpenFile(0) ;
      fileCnt += addFileCnt ;
    }

    bool HasQuality()
    {
      return ( fileType != 0  ) ;
    }

    void Rewind() 
    {
      if ( id != NULL )
        free( id ) ;
      if ( seq != NULL )
        free( seq ) ;
      if ( qual != NULL )
        free( qual ) ;
      id = seq = qual = NULL ;
      currentFpInd = 0 ;

      OpenFile(0) ;
    }

    int Next() 
    {
      //int len ;
      //char buffer[2048] ;
      while ( currentFpInd < fileCnt && ( kseq_read( inSeq ) < 0 ) )
      {
        ++currentFpInd ;
        if (currentFpInd < fileCnt)
          OpenFile(currentFpInd) ;
      }
      if ( currentFpInd >= fileCnt )
        return 0 ;
      /*printf( "%s %s\n", id, inSeq[currentFpInd ]->comment.s ) ;
        if ( inSeq[currentFpInd]->comment.l )
        {
        id[ inSeq[currentFpInd ]->name.l ] = ' ' ;
        strcpy( &id[ inSeq[currentFpInd]->name.l + 1], inSeq[ currentFpInd]->comment.s ) ;
        }*/
      if ( id != NULL )	
        free( id ) ;
      if ( seq != NULL )
        free( seq ) ;
      if ( qual != NULL )
        free( qual ) ;

      id = strdup( inSeq->name.s ) ;
      RemoveReadIdSuffix(id) ;
      seq = strdup( inSeq->seq.s ) ;
      if ( inSeq->qual.l )
        qual = strdup( inSeq->qual.s ) ;
      else
        qual = NULL ;

      return 1 ;
    }

    int NextWithBuffer( char **id, char **seq, char **qual, bool removeReturn = true, bool stopWhenFileEnds = false ) 
    {
      //int len ;
      //char buffer[2048] ;
      while ( currentFpInd < fileCnt && ( kseq_read( inSeq ) < 0 ) )
      {
        ++currentFpInd ;
        if (currentFpInd < fileCnt)
          OpenFile(currentFpInd) ;
        if ( stopWhenFileEnds )
          return -1 ;
      }
      if ( currentFpInd >= fileCnt )
        return 0 ;

      if ( *id != NULL )
        free( *id ) ;
      if ( *seq != NULL )
        free( *seq ) ;
      if ( *qual != NULL )
        free( *qual ) ;

      *id = strdup( inSeq->name.s ) ;
      RemoveReadIdSuffix(*id) ;
      *seq = strdup( inSeq->seq.s ) ;
      /*if ( removeReturn )
        {
        int i ;
        for ( i = strlen( *seq ) - 1 ; i >= 0 ; --i )
        if ( (*seq)[i] != '\n' )
        break ;
        (*seq)[i + 1] = '\0' ;
        }*/
      if ( inSeq->qual.l )
        *qual = strdup( inSeq->qual.s ) ;
      else
        *qual = NULL ;

      return 1 ;
    }

    // Get a batch of reads, it terminates until the buffer is full or 
    // the file ends.
    int GetBatch( struct _Read *readBatch, int maxBatchSize, int &fileInd, bool trimReturn, bool stopWhenFileEnds )
    {
      int batchSize = 0 ;
      while ( batchSize < maxBatchSize ) 
      {
        int tmp = NextWithBuffer( &readBatch[ batchSize].id, &readBatch[batchSize].seq,
            &readBatch[batchSize].qual, trimReturn, stopWhenFileEnds ) ;
        if ( tmp == -1 && batchSize > 0 )
        {
          fileInd = currentFpInd - 1 ;
          return batchSize ; // Finished read current file. The next file is open
        }
        else if ( tmp == -1 && batchSize == 0 )
          continue ; // The current read file is empty
        else if ( tmp == 0 && batchSize == 0 )
        {
          fileInd = currentFpInd ;
          return 0 ; // Finished reading	
        }

        ++batchSize ;
      }

      fileInd = currentFpInd ;
      return batchSize ;
    }

    void FreeBatch(struct _Read *readBatch, int batchSize)
    {
      int i ;
      for (i = 0 ; i < batchSize ; ++i)
      {
        free(readBatch[i].id) ;
        free(readBatch[i].seq) ;
        if (readBatch[i].qual)
          free(readBatch[i].qual) ;
      }
    }

    int GetCurrentFileInd()
    {
      return currentFpInd ;
    }

    std::string GetFileName(int fileInd)
    {
      return fileNames[fileInd] ;
    }

    int GetFileCount()
    {
      return fileCnt ;
    }
} ;

// The class handling read in a batch of reads
/*class ReadBatch
{
} ;*/

#endif

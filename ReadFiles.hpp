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
  char *comment ;
} ;

class ReadFiles
{
  private:
    std::vector<std::string> fileNames ;
    std::vector<bool> hasMate ;
    std::vector<bool> interleaved ; // it is also interleaved 

    gzFile gzFp ;
    kseq_t *inSeq ;
    int fileCnt ;
    int currentFpInd ;
    bool needComment ;
    
    bool opened ;

    std::string specialReadId ;
    bool addSpecialReadForFileEnd ; 
    int fileEndSpecialReadFlag ; // flag:0 hasn't output the special read yet, 1 already output the speical read, so should move to the next file.

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
      if (fileNames[fileInd] != "-")
        gzFp = gzopen( fileNames[fileInd].c_str(), "r" ) ;
      else
        gzFp = gzdopen(fileno(stdin), "r") ;
      inSeq = kseq_init( gzFp ) ;
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
    char *comment ;
    char *seq ;
    char *qual ;

    ReadFiles(): fileCnt(0), currentFpInd(0), opened(false)
    {
      needComment = false ;
      id = comment = seq = qual = NULL ;
      addSpecialReadForFileEnd = false ;
    }

    ~ReadFiles()
    {
      Free() ;
    }

    void Free()
    {
      if ( id != NULL )
        free( id ) ;
      if ( comment != NULL )
        free( comment ) ;
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

    void SetNeedComment(bool in)
    {
      needComment = in ;
    }

    // interleaved file is not frequently set
    void AddReadFile(const char *file, bool fileHasMate, int fileInterleaved = false)
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
        interleaved.push_back(fileInterleaved) ;
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
          interleaved.push_back(fileInterleaved) ;
        }

        addFileCnt = globResult.gl_pathc ;
        globfree(&globResult) ;
      }

      if (fileCnt == 0)
        OpenFile(0) ;
      fileCnt += addFileCnt ;
    }

    void Rewind() 
    {
      if ( id != NULL )
        free( id ) ;
      if ( comment != NULL )
        free( comment ) ;
      if ( seq != NULL )
        free( seq ) ;
      if ( qual != NULL )
        free( qual ) ;
      id = seq = qual = NULL ;
      currentFpInd = 0 ;

      OpenFile(0) ;
    }

    void SetSpecialReadToMarkFileEnd(const char *readId)
    {
      addSpecialReadForFileEnd = true ;
      specialReadId = readId ;
      fileEndSpecialReadFlag = 0 ;
    }
    
    bool HasMate()
    {
      return hasMate[0] ;
    }

    bool IsInterleaved()
    {
      return interleaved[0] ;
    }

    int Next() 
    {
      //int len ;
      //char buffer[2048] ;
      while ( currentFpInd < fileCnt && ( kseq_read( inSeq ) < 0 ) )
      {
        if (addSpecialReadForFileEnd)
        {
          if (fileEndSpecialReadFlag == 0)
          {
            if ( id != NULL )	
              free( id ) ; 
            id = strdup(specialReadId.c_str()) ;
            if (seq != NULL)
              seq[0] = '\0' ;
            fileEndSpecialReadFlag = 1 ;
            return 1 ;
          }
          else
            fileEndSpecialReadFlag = 0 ;
        }
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
      if ( comment != NULL )
        free( comment ) ;
      if ( seq != NULL )
        free( seq ) ;
      if ( qual != NULL )
        free( qual ) ;

      id = strdup( inSeq->name.s ) ;
      RemoveReadIdSuffix(id) ;
      seq = strdup( inSeq->seq.s ) ;
      if ( needComment && inSeq->comment.l )
        comment = strdup(inSeq->comment.s) ;
      else
        comment = NULL ;
      if ( inSeq->qual.l )
        qual = strdup( inSeq->qual.s ) ;
      else
        qual = NULL ;

      return 1 ;
    }

    int NextWithBuffer( char **id, char **seq, char **qual, char **comment, bool removeReturn = true, bool stopWhenFileEnds = false ) 
    {
      //int len ;
      //char buffer[2048] ;
      while ( currentFpInd < fileCnt && ( kseq_read( inSeq ) < 0 ) )
      {
        if (addSpecialReadForFileEnd)
        {
          if (fileEndSpecialReadFlag == 0)
          {
            if ( *id != NULL )	
              free( *id ) ; 
            *id = strdup(specialReadId.c_str()) ;
            if (*seq != NULL)
              (*seq)[0] = '\0' ;
            fileEndSpecialReadFlag = 1 ;
            return 1 ;
          }
          else
            fileEndSpecialReadFlag = 0 ;
        }

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
      if ( *comment != NULL )
        free( *comment ) ;
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
      if ( needComment && inSeq->comment.l )
        *comment = strdup(inSeq->comment.s) ;
      else
        *comment = NULL ;
      if ( inSeq->qual.l )
        *qual = strdup( inSeq->qual.s ) ;
      else
        *qual = NULL ;

      return 1 ;
    }

    // Get a batch of reads, it terminates until the buffer is full or 
    // the file ends.
    // readBatch2 can be for interleaved file. 
    int GetBatch( struct _Read *readBatch, int maxBatchSize, int &fileInd, bool trimReturn, bool stopWhenFileEnds, struct _Read *readBatch2 = NULL)
    {
      int batchSize = 0 ;
      while ( batchSize < maxBatchSize ) 
      {
        int tmp = NextWithBuffer( &readBatch[ batchSize].id, &readBatch[batchSize].seq,
            &readBatch[batchSize].qual, &readBatch[batchSize].comment,
            trimReturn, stopWhenFileEnds ) ;
        
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
        
        if (readBatch2 != NULL)
          tmp = NextWithBuffer( &readBatch2[ batchSize].id, &readBatch2[batchSize].seq,
              &readBatch2[batchSize].qual, &readBatch2[batchSize].comment,
              trimReturn, stopWhenFileEnds ) ;

        ++batchSize ;
      }

      fileInd = currentFpInd ;
      return batchSize ;
    }

    void CopyBatch(struct _Read *to, struct _Read *from, int batchSize)
    {
      int i ;
      for (i = 0 ; i < batchSize ; ++i)
      {
        free(to[i].id) ;
        free(to[i].seq) ;
        if (needComment && to[i].comment)
          free(to[i].comment) ;
        if (to[i].qual)
          free(to[i].qual) ;
        
        to[i].id = strdup(from[i].id) ;
        to[i].seq = strdup(from[i].seq) ;
        if (needComment && from[i].comment)
          to[i].comment = strdup(from[i].comment) ;
        else
          to[i].comment = NULL ;
        if (from[i].qual)
          to[i].qual = strdup(from[i].qual) ;
        else
          to[i].qual = NULL ;
      }
    }

    void FreeBatch(struct _Read *readBatch, int batchSize)
    {
      int i ;
      for (i = 0 ; i < batchSize ; ++i)
      {
        free(readBatch[i].id) ;
        free(readBatch[i].seq) ;
        if (readBatch[i].comment)
          free(readBatch[i].comment) ;
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

#endif

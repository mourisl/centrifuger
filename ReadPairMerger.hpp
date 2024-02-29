#ifndef _MOURISL_READ_PAIR_MERGER
#define _MOURISL_READ_PAIR_MERGER

#include <stdio.h>
#include <string.h>

class ReadPairMerger
{
private:
  char _compChar[256] ;
  bool _checkReadThrough ;

	int IsMateOverlap( char *fr, int flen, char *sr, int slen, int minOverlap, int &offset, int &bestMatchCnt,
		bool checkTandem = true )
  {
    int i, j, k ;		
    bestMatchCnt = -1 ;
    int offsetCnt = 0 ;
    int overlapSize = -1 ;
    for ( j = 0 ; j < flen - minOverlap ; ++j ) // The overlap start position in first read
    {
      // Whether the overlap works.
      int matchCnt = 0 ;
      bool flag = true ;

      double similarityThreshold = 0.95 ;
      if ( flen - j >= 100 )
        similarityThreshold = 0.85 ;
      else if ( flen - j >= 50 )
        similarityThreshold = 0.85 + ( flen - j - 50 ) / 50.0 * 0.1 ;

      for ( k = 0 ; j + k < flen && k < slen ; ++k )
      {
        if ( fr[j + k] == sr[k] )
          ++matchCnt ;
        if ( matchCnt + ( flen - ( j + k ) - 1 ) < int( ( flen - j ) * similarityThreshold ) )
        {
          flag = false ;
          break ;
        }
      }

      if ( flag ) 
      {
        offset = j ;
        ++offsetCnt ;
        overlapSize = k ;
        bestMatchCnt = matchCnt ;
      }
    }

    if ( offsetCnt != 1 )
      return -1 ;

    // If the overlap size is near minOverlap, the overlap could still be ambiguous. 
    // i- the repeat size
    if ( checkTandem && overlapSize <= minOverlap * 2 )
    {
      for ( i = 1 ; i <= overlapSize / 2 ; ++i )
      {
        bool tandem = true ;
        for ( j = i ; j + i - 1 < overlapSize ; j += i )
        {
          for ( k = j ; k <= j + i - 1 ; ++k )
          {
            if ( sr[k -j] != sr[k] )
              break ;
          }
          if ( k <= j + i - 1 )
          {
            tandem = false ;
            break ;
          }
        }

        if ( tandem )
          return -1 ;
      }
    }

    return overlapSize ;
  }

  void ReverseBuffer(char *buffer, int len)
  {
    int i, j ;
    for (i = 0, j = len - 1 ; i < j ; ++i, --j )
    {
      char tmp = buffer[i] ;
      buffer[i] = buffer[j] ;
      buffer[j] = tmp ;
    }
  }

  void ComplementBuffer(char *buffer, int len)
  {
    int i ;
    for (i = 0 ; i < len ; ++i)
      buffer[i] = _compChar[ (int)buffer[i] ] ;
  }

public:
  ReadPairMerger()
  {
    int i ;
    for (i = 0 ; i < 256 ; ++i)
      _compChar[i] = 'N' ;
    _compChar['A'] = 'T' ;
    _compChar['C'] = 'G' ;
    _compChar['G'] = 'C' ;
    _compChar['T'] = 'A' ;

    _checkReadThrough = true ;
  }

  ~ReadPairMerger()
  {
  }

  void SetCheckReadThrough(bool check)
  {
    _checkReadThrough = check ;
  }

  // r1,q1: read and quality for mate 1.
  // r2, q2: for read 2
  // rm, qm: merged read and quality score
  // @ret: 0 no merge. rm and qm will be NULL
  //       1 regular merge
  //       2 read through
  //       Also overlapSize, offset, bestMatchCnt describes the overlap statistics
  int Merge(char *r1, char *q1, char *r2, char *q2, char **rm, char **qm, int &overlapSize, int &offset, int&bestMatchCnt) 
  {
    int i ;
    int len1 = strlen(r1) ;
    int len2 = strlen(r2) ;
    
    *rm = NULL ;
    *qm = NULL ;
    if (r2 == NULL)
      return 0 ;
    
    char *rcr2 = strdup(r2) ;
    ReverseBuffer(rcr2, len2) ;
    ComplementBuffer(rcr2, len2) ;
    
    char *rcq2 = NULL ;
    if (q2 != NULL)
    {
      rcq2 = strdup(q2) ;
      ReverseBuffer(rcq2, len2) ;
    }

    int minOverlap = ( len1 + len2 ) / 10 ; // overlap required for read trhough
    int minOverlap2 = ( len1 + len2 ) / 10 ; // overlap required for directly overlap
    if ( minOverlap > 31 )
      minOverlap = 31 ;
    if ( minOverlap2 > 31 )
      minOverlap2 = 31 ;
    offset = -1 ;
    bestMatchCnt = -1 ;

    // Read through
    overlapSize = IsMateOverlap(rcr2, len2, r1, len1, minOverlap, offset, bestMatchCnt, false ) ;
    if (overlapSize >= 0)
    {
      *rm = (char *)malloc(sizeof(char) * (overlapSize + 1)) ;
      memcpy(*rm, r1, overlapSize) ;
      (*rm)[overlapSize] = '\0' ;
      
      if (q1 != NULL)
      {
        *qm = (char *)malloc(sizeof(char) * (overlapSize + 1)) ;
        memcpy(*qm, q1, overlapSize) ;
        (*qm)[overlapSize] ='\0' ;

        for (i = 0 ; i < overlapSize ; ++i)
        {
          if (rcq2[i + offset] > q1[i] || (*rm)[i] == 'N')
          {
            (*rm)[i] = rcr2[i + offset] ;
            (*qm)[i] = rcq2[i + offset] ;
          }
        }
      }
      
      free(rcr2) ;
      if (rcq2)
        free(rcq2) ;
      return 2 ;
    }
    
    // Simple overlap
    overlapSize = IsMateOverlap(r1, len1, rcr2, len2, minOverlap2, offset, bestMatchCnt, true) ;
    if (overlapSize >= 0)
    {
      *rm = (char *)malloc(sizeof(char) * (len1 + len2 - overlapSize + 1)) ;
      if (rcq2)
        *qm = (char *)malloc(sizeof(char) * (len1 + len2 - overlapSize + 1)) ;
      
      for (i = 0 ; i < len2 ; ++i)
      {
        (*rm)[offset + i] = rcr2[i] ;
        if (rcq2 != NULL)
          (*qm)[offset + i] = rcq2[i] ;
      }
      
      int len = offset + i ; // Potentially for r2 is a substirng of r1
      for (i = 0 ; i < len1 && i < len ; ++i)
      {
        if (i < offset || (q1 != NULL && q1[i] >= (*qm)[i] - 14) || (*rm)[i] == 'N')        
        {
          (*rm)[i] = r1[i] ;
          if (q1 != NULL)
            (*qm)[i] = q1[i] ;
        }
      }
      (*rm)[len] = '\0' ;
      if (rcq2)
        (*qm)[len] = '\0' ;

      free(rcr2) ;
      if (rcq2)
        free(rcq2) ;

      return 1 ;
    }
    
    free(rcr2) ;
    if (rcq2)
      free(rcq2) ;
    return 0 ;
  }

  // A wrapper that don't return the overlap statistics
  int Merge(char *r1, char *q1, char *r2, char *q2, char **rm, char **qm) 
  {
    int overlapSize, matchCnt, offset ;
    return Merge(r1, q1, r2, q2, rm, qm, overlapSize, matchCnt, offset) ;
  }
} ;

#endif

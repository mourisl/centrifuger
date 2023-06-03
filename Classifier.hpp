#ifndef _MOURISL_CLASSIFIER_HEADER
#define _MOURISL_CLASSIFIER_HEADER

#include <string.h>

#include "Taxonomy.hpp"
#include "compactds/FMIndex.hpp"
#include "compactds/Sequence_Hybrid.hpp"

struct _classifierParam 
{
  int maxResult ; // the number of entries in the results    
  int minHitLen ;
  _classifierParam()
  {
    maxResult = 1 ;
    minHitLen = 15 ;
  }
} ;

struct _classifierResult
{
  int score ;
  int secondaryScore ;
  std::vector<size_t> taxHits ;
  std::vector<size_t> seqHits ;
} ;

class Classifier
{
private:
  FMIndex<Sequence_Hybrid> _fm ;
  Taxonomy _taxonomy ;
  std::map<size_t, size_t> _seqLength ;
  _classifierParam _param ;

public:
  Classifier() {}
  ~Classifier() {Free() ;}

  void Free()
  {
    _fm.Free() ;
    _taxonomy.Free() ;
    _seqLength.clear() ;
  }

  void Init(char *idxPrefix, struct _classifierParam param)
  {
    FILE *fp ;
    char *nameBuffer = (char *)malloc(sizeof(char) * (strlen(idxPrefix) + 17))  ;  
    
    sprintf(nameBuffer, "%s.1.cfr", idxPrefix) ;
    fp = fopen(nameBuffer, "r") ;
    _fm.Load(fp) ;
    fclose(fp) ;

    // .2.cfr file is for taxonomy structure
    sprintf(nameBuffer, "%s.2.cfr", idxPrefix) ;
    fp = fopen(nameBuffer, "r") ;
    _taxonomy.Load(fp) ;
    fclose(fp) ;

    // .3.cfr file is for sequence length
    sprintf(nameBuffer, "%s.3.cfr", idxPrefix) ;
    fp = fopen(nameBuffer, "r") ;
    size_t tmp[2] ;
    while (fread(tmp, sizeof(tmp[0]), 2, fp))
    {
      _seqLength[tmp[0]] = tmp[1] ;
    }
    fclose(fp) ;
    
    _param = param ;

    free(nameBuffer) ;
  }

  // Main function to return the classification results 
  void Query(char *r1, char *r2, struct _classifierResult results)
  {
    size_t i ;
    int readLen = strlen(r1) ;
    size_t sp = 0, ep = 0, l = 0 ;
    l = _fm.BackwardSearch(r1, readLen, sp, ep) ;
    printf("Backward search %d %d %d\n", (int)l, (int)sp, (int)ep) ;
    for (i = sp ; i <= ep ; ++i)
    {
      size_t searchl ;
      size_t sa = _fm.BackwardToSampledSA(i, searchl) ;
      printf("SA[%d] = %d+%d\n", (int)i, (int)sa, (int)searchl) ;
    }
  }
} ;

#endif 

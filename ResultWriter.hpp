#ifndef _MOURISL_RESULTWRITER_HEADER 
#define _MOURISL_RESULTWRITER_HEADER 

#include "Classifier.hpp" 

class ResultWriter
{
private:
  FILE *fpClassification ;
public:
  //ResultWriter(const Taxonomy taxonomy): _taxonomy(taxonomy)  
  ResultWriter() 
  {
    fpClassification = stdout ;
  }

  ~ResultWriter() 
  {
    if (fpClassification != stdout)
      fclose(fpClassification) ;
  }
  
  void SetClassificationOutput(const char *filename)
  {
    fpClassification = fopen(filename, "w") ;
  }

  void OutputHeader()
  {
    printf("readID\tseqID\ttaxID\tscore\t2ndBestScore\thitLength\tqueryLength\tnumMatches\n") ;
  }

  void Output(char *readid, char *barcode, const struct _classifierResult &r)
  {
    int i ;
    int matchCnt = r.taxIds.size() ;
    if (matchCnt > 0)
    {
      for (i = 0 ; i < matchCnt ; ++i)
      {
        fprintf(fpClassification,
            "%s\t%s\t%lu\t%lu\t%lu\t%d\t%d\t%d\n",
            readid, r.seqStrNames[i].c_str(), r.taxIds[i],
            r.score, r.secondaryScore, r.hitLength, r.queryLength, matchCnt) ;
      }
    }
    else
    {
      fprintf(fpClassification,
          "%s\tunclassified\t0\t0\t0\t0\t%d\t1\n", readid, r.queryLength) ;
    }
  }

  void Finalize()
  {
  }
} ;

#endif

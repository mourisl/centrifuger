#ifndef _MOURISL_RESULTWRITER_HEADER 
#define _MOURISL_RESULTWRITER_HEADER 

#include "Classifier.hpp" 

#include "ReadFormatter.hpp"
#include "BarcodeCorrector.hpp"
#include "BarcodeTranslator.hpp"
#include "ReadFiles.hpp"

class ResultWriter
{
private:
  FILE *_fpClassification ;
  bool _hasBarcode ;
  bool _hasUmi ;
  bool _outputUnclassified ;
  bool _outputClassified ;
  gzFile _gzFpUnclassified[4] ; 
  gzFile _gzFpClassified[4] ;

  size_t _classifiedCnt ;
  size_t _totalCnt ;

  // Variable to handle output from sample-sheet
  std::vector< std::string > _multiOutputFileList ;
  bool _hasSpecialReadIdForFileEnd ;
  int _currentMultiOutputFile ;
  std::map< std::string, int > _multiOutputFileMap ;

  void PrintExtraCol(const char *s) 
  {
    if (s == NULL)
      fprintf(_fpClassification, "\t") ;
    else
      fprintf(_fpClassification, "\t%s", s) ;
  }
public:
  //ResultWriter(const Taxonomy taxonomy): _taxonomy(taxonomy)  
  ResultWriter() 
  {
    _fpClassification = stdout ;
    _outputUnclassified = false ;
    _outputClassified = false ;
    _hasBarcode = false ;
    _hasUmi = false ;
    _hasSpecialReadIdForFileEnd = false ;

    int i ;
    for (i = 0 ; i < 4 ; ++i)
    {
      _gzFpUnclassified[i] = NULL ;
      _gzFpClassified[i] = NULL ;
    }

    _classifiedCnt = _totalCnt = 0 ;
  }

  ~ResultWriter() 
  {
    if (_fpClassification != stdout && _fpClassification != NULL)
      fclose(_fpClassification) ;
    int i ;
    for (i = 0 ; i < 4 ; ++i)
    {
      if (_gzFpUnclassified[i])
        gzclose(_gzFpUnclassified[i]) ;
      if (_gzFpClassified[i])
        gzclose(_gzFpClassified[i]) ;
    }
  }

  void SetMultiOutputFileList(std::vector< std::string > &filenames)
  {
    _multiOutputFileList = filenames ;
    SetClassificationOutput(filenames[0].c_str(), "w") ;

    _currentMultiOutputFile = 0 ;
    _multiOutputFileMap[ _multiOutputFileList[0] ] = 1 ;
    _hasSpecialReadIdForFileEnd = true ;
  }

  //@return: open file mode
  char NextMultiOutputFile()
  {
    if (_fpClassification != NULL)
    {
      fclose( _fpClassification ) ;
      _fpClassification = NULL ;
    }
    ++_currentMultiOutputFile ;
    
    if (_currentMultiOutputFile >= (int)_multiOutputFileList.size()) 
      return 'e' ; // end

    char mode[2] = "w" ;
    if (_multiOutputFileMap.find(_multiOutputFileList[ _currentMultiOutputFile ]) != _multiOutputFileMap.end() )
      mode[0] = 'a' ;
    SetClassificationOutput(_multiOutputFileList[ _currentMultiOutputFile ].c_str(), mode ) ;
    if (mode[0] == 'w') 
    {
      _multiOutputFileMap[ _multiOutputFileList[ _currentMultiOutputFile ] ] = 1 ;
    }
    return mode[0] ;
  }
  
  void SetClassificationOutput(const char *filename, const char *mode)
  {
    _fpClassification = fopen(filename, mode) ;
  }

  // category: 0: unclassified reads, 1: classified reads
  void SetOutputReads(const char *prefix, bool hasMate, bool hasBarcode, bool hasUmi, int category)
  {
    int len = strlen(prefix) ;      
    char extension[10] = "" ;
    char *name = (char *)malloc(sizeof(char) * (len + 1 + 10)) ;
    
    gzFile *gzFps = _gzFpUnclassified ;
    if (category == 0)
      _outputUnclassified = true ;
    else 
    {
      gzFps = _gzFpClassified ;
      _outputClassified = true ;
    }
    
    // Add "fa" or "fq" to the name
    // Now always fq.
    extension[0] = '.' ;
    extension[1] = 'f' ;
    extension[2] = 'q' ;

    // Add "gz" to the name
    extension[3] = '.' ;
    extension[4] = 'g' ;
    extension[5] = 'z' ;
    extension[6] = '\0' ;

    if (hasMate)
    {
      sprintf(name, "%s_1%s", prefix, extension) ;
      gzFps[0] = gzopen(name, "w1") ;
     
      sprintf(name, "%s_2%s", prefix, extension) ;
      gzFps[1] = gzopen(name, "w1") ;
    }
    else
    {
      sprintf(name, "%s%s", prefix, extension) ;
      gzFps[0] = gzopen(name, "w1") ;
    }

    extension[2] = 'a' ; // always 'fa' for barcode and umi
    if (hasBarcode)
    {
      sprintf(name, "%s_bc%s", prefix, extension) ;
      gzFps[2] = gzopen(name, "w1") ;
    }
    if (hasUmi)
    {
      sprintf(name, "%s_um%s", prefix, extension) ;
      gzFps[3] = gzopen(name, "w1") ;
    }

    free(name) ;
  }

  void SetHasBarcode(bool s)
  {
    _hasBarcode = s ;
  }

  void SetHasUmi(bool s) 
  {
    _hasUmi = s ;
  }

  void OutputHeader()
  {
    fprintf(_fpClassification, "readID\tseqID\ttaxID\tscore\t2ndBestScore\thitLength\tqueryLength\tnumMatches") ;
  
    if (_hasBarcode)
      fprintf(_fpClassification, "\tbarcode") ;
    if (_hasUmi)
      fprintf(_fpClassification, "\tUMI") ;
    fprintf(_fpClassification, "\n") ;
  }

  void Output(const char *readid, 
      const char *seq1, const char *qual1, const char *seq2, const char *qual2,
      const char *barcode, const char *umi, const struct _classifierResult &r)
  {
    if (_hasSpecialReadIdForFileEnd && !strcmp(readid, SAMPLE_SHEET_SEPARATOR_READ_ID))
    {
      if (NextMultiOutputFile() == 'w') 
        OutputHeader() ;
      return ;
    }

    int i ;
    int matchCnt = r.taxIds.size() ;
    ++_totalCnt ;
    if (matchCnt > 0)
    {
      ++_classifiedCnt ;
      for (i = 0 ; i < matchCnt ; ++i)
      {
        fprintf(_fpClassification,
            "%s\t%s\t%lu\t%lu\t%lu\t%d\t%d\t%d",
            readid, r.seqStrNames[i].c_str(), r.taxIds[i],
            r.score, r.secondaryScore, r.hitLength, r.queryLength, matchCnt) ;
        if (_hasBarcode)
          PrintExtraCol(barcode) ;
        if (_hasUmi)
          PrintExtraCol(umi) ;
        fprintf(_fpClassification, "\n") ;
      }
    }
    else
    {
      fprintf(_fpClassification,
          "%s\tunclassified\t0\t0\t0\t0\t%d\t1", readid, r.queryLength) ;
      if (_hasBarcode)
        PrintExtraCol(barcode) ;
      if (_hasUmi)
        PrintExtraCol(umi) ;
      fprintf(_fpClassification, "\n") ;
    }

    for (i = 0 ; i <= 1 ; ++i)
    {
      gzFile *gzFps = NULL ;
      if ( i == 0 && matchCnt == 0 && _outputUnclassified)
        gzFps = _gzFpUnclassified ;
      else if (i == 1 && matchCnt > 0 && _outputClassified)
        gzFps = _gzFpClassified ;
      else
        continue ;
      
      if (qual1 == NULL)
        gzprintf(gzFps[0], ">%s\n%s\n", readid, seq1) ;
      else
        gzprintf(gzFps[0], "@%s\n%s\n+\n%s\n", readid, seq1, qual1) ;
      
      if (seq2 != NULL) 
      {
        if (qual2 == NULL)
          gzprintf(gzFps[1], ">%s\n%s\n", readid, seq2) ;
        else
          gzprintf(gzFps[1], "@%s\n%s\n+\n%s\n", readid, seq2, qual2) ;
      }

      if (_hasBarcode)
      {
        gzprintf(gzFps[2], ">%s\n%s\n", readid, barcode) ;
      }

      if (_hasUmi)
      {
        gzprintf(gzFps[3], ">%s\n%s\n", readid, umi) ;
      }
    }
  }

  void Finalize()
  {
    Utils::PrintLog("Processed %lu read fragments, and %lu (%.2lf\%) can be classified.",
        _totalCnt, _classifiedCnt, (double)_classifiedCnt / (double)_totalCnt * 100.0) ;
  }
} ;

#endif

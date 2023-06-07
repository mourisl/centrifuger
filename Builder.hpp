#ifndef _MOURISL_BUILDER_HEADER
#define _MOURISL_BUILDER_HEADER

#include "ReadFiles.hpp"
#include "compactds/Sequence_Hybrid.hpp"
#include "compactds/FMBuilder.hpp"
#include "compactds/FMIndex.hpp"
#include "compactds/Alphabet.hpp"
#include "compactds/SequenceCompactor.hpp"
#include "Taxonomy.hpp"

// Holds various method regarding building index
using namespace compactds ; 

class Builder
{
private:
  FMIndex<Sequence_Hybrid> _fmIndex ;
  Taxonomy _taxonomy ;
  std::map<size_t, size_t> _seqLength ; // we use map here is for the case that a seq show up in the conversion table but not in the actual genome file.

  void TransformSampledSAToSeqId(struct _FMBuilderParam &fmBuilderParam, std::vector<size_t> genomeSeqIds,
      std::vector<size_t> genomeLens)
  {
    size_t i ;
    PartialSum lenPsum ;
    lenPsum.Init(genomeLens.data(), genomeLens.size()) ;
    for (i = 0 ; i < fmBuilderParam.sampleSize ; ++i)
    {
      fmBuilderParam.sampledSA[i] = genomeSeqIds[lenPsum.Search(fmBuilderParam.sampledSA[i])] ;   
    }
    fmBuilderParam.adjustedSA0 = genomeSeqIds[0] ;
  }

public: 
  Builder() {}
  ~Builder() 
  {
    _fmIndex.Free() ;
    _taxonomy.Free() ;
  }

  void Build(ReadFiles &refGenomeFile, char *taxonomyFile, char *nameTable, char *conversionTable, struct _FMBuilderParam &fmBuilderParam, const char *alphabetList)
  {
    const int alphabetSize = strlen(alphabetList) ;
  
    _taxonomy.Init(taxonomyFile, nameTable, conversionTable)  ; 

    FixedSizeElemArray genomes ;
    SequenceCompactor seqCompactor ;
    seqCompactor.Init(alphabetList, genomes, 1000000) ;
    
    std::vector<size_t> genomeSeqIds ;
    std::vector<size_t> genomeLens ; 
    while (refGenomeFile.Next())
    {
      size_t seqid = _taxonomy.SeqNameToId(refGenomeFile.id) ;
      if (seqid >= _taxonomy.GetSeqCount())
      {
        fprintf(stderr, "WARNING: taxonomy id doesn't exist for %s!\n", refGenomeFile.id) ;
        seqid = _taxonomy.AddExtraSeqName(refGenomeFile.id) ;
      }

      size_t len = seqCompactor.Compact(refGenomeFile.seq, genomes) ;
      
      _seqLength[seqid] = len ;
      genomeSeqIds.push_back(seqid) ;
      genomeLens.push_back(len) ;
    }

    FixedSizeElemArray BWT ;
    size_t firstISA ;
    FMBuilder::Build(genomes, genomes.GetSize(), alphabetSize, BWT, firstISA, fmBuilderParam) ;
    
    TransformSampledSAToSeqId(fmBuilderParam, genomeSeqIds, genomeLens) ;
    _fmIndex.Init(BWT, genomes.GetSize(), 
        firstISA, fmBuilderParam, alphabetList, alphabetSize) ;
  }

  void Save(const char *outputPrefix)
  {
    char outputFileName[1024] ; 
    // Convert the sampled point to seqID.
    FILE *fpOutput ;
    // .1.cfr file is for the index
    sprintf(outputFileName, "%s.1.cfr", outputPrefix) ;
    fpOutput = fopen(outputFileName, "w") ;
    _fmIndex.Save(fpOutput) ;
    fclose(fpOutput) ;

    // .2.cfr file is for taxonomy structure
    sprintf(outputFileName, "%s.2.cfr", outputPrefix) ;
    fpOutput = fopen(outputFileName, "w") ;
    _taxonomy.Save(fpOutput) ;
    fclose(fpOutput) ;

    // .3.cfr file is for sequence length
    sprintf(outputFileName, "%s.3.cfr", outputPrefix) ;
    fpOutput = fopen(outputFileName, "w") ;
    for (std::map<size_t, size_t>::iterator iter = _seqLength.begin() ; 
        iter != _seqLength.end() ; ++iter)
    {
      size_t tmp[2] = {iter->first, iter->second} ;
      fwrite(tmp, sizeof(tmp[0]), 2, fpOutput) ;
    }
    fclose(fpOutput) ;

  }
} ;

#endif

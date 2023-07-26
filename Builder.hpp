#ifndef _MOURISL_BUILDER_HEADER
#define _MOURISL_BUILDER_HEADER

#include "ReadFiles.hpp"
#include "compactds/Sequence_Hybrid.hpp"
#include "compactds/Sequence_RunBlock.hpp"
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
  FMIndex<Sequence_RunBlock> _fmIndex ;
  Taxonomy _taxonomy ;
  std::map<size_t, size_t> _seqLength ; // we use map here is for the case that a seq show up in the conversion table but not in the actual genome file.

  // SampledSA need to be processed before FMIndex.Init() because the sampledSA is represented by FixedElemLengthArray, which requires the largest element size
  void TransformSampledSAToSeqId(struct _FMBuilderParam &fmBuilderParam, std::vector<size_t> genomeSeqIds,
      std::vector<size_t> genomeLens, size_t n)
  {
    size_t i ;
    PartialSum lenPsum ;
    lenPsum.Init(genomeLens.data(), genomeLens.size()) ;
    for (i = 0 ; i < fmBuilderParam.sampleSize ; ++i)
    {
      // The precomputeWidth + 1 here to handle the fuzzy boundary
      if (fmBuilderParam.sampledSA[i] + fmBuilderParam.precomputeWidth + 1 < n)
        fmBuilderParam.sampledSA[i] = genomeSeqIds[lenPsum.Search(
            fmBuilderParam.sampledSA[i] + fmBuilderParam.precomputeWidth + 1)] ;  
      else
        fmBuilderParam.sampledSA[i] = genomeSeqIds[lenPsum.Search(
            fmBuilderParam.sampledSA[i])] ;
    }
    fmBuilderParam.adjustedSA0 = genomeSeqIds[0] ;

    for (std::map<size_t, size_t>::iterator iter = fmBuilderParam.selectedSA.begin() ;
        iter != fmBuilderParam.selectedSA.end() ; ++iter)
    {
      iter->second = genomeSeqIds[lenPsum.Search(iter->second + fmBuilderParam.precomputeWidth + 1)] ; // the selected SA stores the fuzzy start position for the next genome, so we need to plus the adjusted boundary.
    }
  }

public: 
  Builder() {}
  ~Builder() 
  {
    _fmIndex.Free() ;
    _taxonomy.Free() ;
  }

  void Build(ReadFiles &refGenomeFile, char *taxonomyFile, char *nameTable, char *conversionTable, uint64_t subsetTax, struct _FMBuilderParam &fmBuilderParam, const char *alphabetList)
  {
    size_t i ;
    const int alphabetSize = strlen(alphabetList) ;
  
    _taxonomy.Init(taxonomyFile, nameTable, conversionTable)  ; 

    FixedSizeElemArray genomes ;
    SequenceCompactor seqCompactor ;
    seqCompactor.Init(alphabetList, genomes, 1000000) ;
    
    std::map<size_t, int> selectedTaxIds ;
    if (subsetTax != 0)
      _taxonomy.GetChildrenTax(_taxonomy.CompactTaxId(subsetTax), selectedTaxIds) ; 
    std::vector<size_t> genomeSeqIds ;
    std::vector<size_t> genomeLens ; 
    while (refGenomeFile.Next())
    {
      size_t seqid = _taxonomy.SeqNameToId(refGenomeFile.id) ;
      
      if (subsetTax != 0)
      {
        size_t taxid = _taxonomy.SeqIdToTaxId(seqid) ;
        if (selectedTaxIds.find(taxid) == selectedTaxIds.end())
          continue ;
      }
      
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
    
    // Put in the genome boundary information for selected SA in our FM index
    size_t genomeCnt = genomeLens.size() ;
    if (genomeCnt == 0)
    {
      fprintf(stderr, "ERROR: found 0 genomes in the input or after filtering.\n") ;
      exit(EXIT_FAILURE) ;
    }
    size_t psum = 0 ; // genome length partial sum
    for (i = 0 ; i < genomeCnt - 1 ; ++i)
    {
      psum += genomeLens[i] ;
      if (psum < (size_t)fmBuilderParam.precomputeWidth + 1ull) // default: 12bp fuzzy boundary
        continue ;
      fmBuilderParam.selectedISA[psum - fmBuilderParam.precomputeWidth - 1] ;
    }

    Utils::PrintLog("Found %lu sequences with total length %lu bp", 
        genomeCnt, genomes.GetSize()) ;

    FMBuilder::Build(genomes, genomes.GetSize(), alphabetSize, BWT, firstISA, fmBuilderParam) ;
    TransformSampledSAToSeqId(fmBuilderParam, genomeSeqIds, genomeLens, genomes.GetSize()) ;
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

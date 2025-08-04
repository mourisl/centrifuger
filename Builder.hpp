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

  void SetRBBWTBlockSize(size_t b)
  {
    _fmIndex.SetSequenceExtraParameter((void *)b) ;
  }

  void Build(ReadFiles &refGenomeFile, char *taxonomyFile, char *nameTable, char *conversionTable, bool conversionTableAtFileLevel, bool concatSameTaxIdSeqs, bool ignoreUncategorizedSeqs, uint64_t subsetTax, size_t memoryConstraint, struct _FMBuilderParam &fmBuilderParam, const char *alphabetList)
  {
    size_t i ;
    const int alphabetSize = strlen(alphabetList) ;
  
    _taxonomy.Init(taxonomyFile, nameTable, conversionTable, conversionTableAtFileLevel)  ; 

    FixedSizeElemArray genomes ;
    std::map<size_t, FixedSizeElemArray *> taxIdGenomes ; // the genomes from each tax ID. For the concatSameTaxIdSeqs option.
    SequenceCompactor seqCompactor ;
    seqCompactor.Init(alphabetList, genomes, 1000000) ;
    
    std::map<size_t, int> selectedTaxIds ;
    if (subsetTax != 0)
      _taxonomy.GetChildrenTax(_taxonomy.CompactTaxId(subsetTax), selectedTaxIds) ; 
    std::vector<size_t> genomeSeqIds ;
    std::vector<size_t> genomeLens ; 
    while (refGenomeFile.Next())
    {
      size_t seqid = 0 ;
      char fileNameBuffer[1024] ;
      if (conversionTableAtFileLevel)
      {
        Utils::GetFileBaseName(refGenomeFile.GetFileName( refGenomeFile.GetCurrentFileInd() ).c_str(), 
            "fna|fa|fasta|faa", fileNameBuffer) ;
        seqid = _taxonomy.SeqNameToId(fileNameBuffer) ;
      }
      else
        seqid = _taxonomy.SeqNameToId(refGenomeFile.id) ;

      if (subsetTax != 0)
      {
        size_t taxid = _taxonomy.SeqIdToTaxId(seqid) ;
        if (selectedTaxIds.find(taxid) == selectedTaxIds.end())
          continue ;
      }

      if (seqid >= _taxonomy.GetSeqCount())
      {
        fprintf(stderr, "WARNING: taxonomy id doesn't exist for %s!\n", 
            conversionTableAtFileLevel ? fileNameBuffer : refGenomeFile.id) ;
        if (!ignoreUncategorizedSeqs)
          seqid = _taxonomy.AddExtraSeqName(conversionTableAtFileLevel ? fileNameBuffer : refGenomeFile.id) ;
        else
          continue ;
      }

      if (!concatSameTaxIdSeqs)
      {
        size_t len = seqCompactor.Compact(refGenomeFile.seq, genomes) ;
        if (len < fmBuilderParam.precomputeWidth + 1ull) // A genome too short
        {
          fprintf(stderr, "WARNING: %s is filtered due to its short length (could be from masker)!\n", refGenomeFile.id) ;
          size_t size = genomes.GetSize() ;
          genomes.SetSize(size - len) ;
          continue ;
        }

        if (_seqLength.find(seqid) == _seqLength.end()) // Assume there is no duplicated seqid
        {
          _seqLength[seqid] = len ;
          genomeSeqIds.push_back(seqid) ;
          genomeLens.push_back(len) ;
        }
        else
        {
          _seqLength[seqid] += len ;
          genomeLens[ genomeLens.size() - 1 ] += len ;
        }
      }
      else // concatenate seuqences with the same tax ID
      {
        size_t taxid = _taxonomy.SeqIdToTaxId(seqid) ;
        if (taxIdGenomes.find(taxid) == taxIdGenomes.end())
        {
          FixedSizeElemArray *a = new FixedSizeElemArray ;
          seqCompactor.Init(alphabetList, *a, 10000) ;
          taxIdGenomes[taxid] = a ;
        }

        FixedSizeElemArray *a = taxIdGenomes[taxid] ;
        size_t len = seqCompactor.Compact(refGenomeFile.seq, *a) ;
        if (len < fmBuilderParam.precomputeWidth + 1ull) // A genome too short
        {
          fprintf(stderr, "WARNING: %s is filtered due to its short length (could be from masker)!\n", refGenomeFile.id) ;
          size_t size = a->GetSize() ;
          a->SetSize(size - len) ;
          continue ;
        }
      }
    }

    if (concatSameTaxIdSeqs)
    {
      genomes.SetSize(0) ;
      
      // In this case, seqId essentially is taxId
      _taxonomy.SetTaxIdAsSeqId() ;

      for ( std::map<size_t, FixedSizeElemArray *>::iterator iter = taxIdGenomes.begin() ;
          iter != taxIdGenomes.end() ; ++iter)
      {
        if (iter->second->GetSize() == 0)
        {
          delete iter->second ;
          continue ;
        }
        genomes.PushBack(*(iter->second), iter->second->GetSize()) ;
        genomeSeqIds.push_back(iter->first) ;
        genomeLens.push_back(iter->second->GetSize()) ;
        _seqLength[iter->first] = iter->second->GetSize() ;
        delete iter->second ;
      }
      Utils::PrintLog("Finish concatenating genomes") ; 
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
    
    size_t totalGenomeSize = genomes.GetSize() ;
    Utils::PrintLog("Found %lu sequences with total length %lu bp.", 
        genomeCnt, totalGenomeSize) ;
    
    if (memoryConstraint != 0)
      FMBuilder::InferParametersGivenMemory(genomes.GetSize(), alphabetSize, memoryConstraint,fmBuilderParam) ;
    /*{
      size_t memoryCost = totalGenomeSize / 4 / WORDBYTES + DIV_CEIL(totalGenomeSize, fmBuilderParam.sampleRate) * Utils::Log2Ceil(genomeCnt) ; // We need to substract other portion out, like the space for the rank structure and the reduced sampled SA. ;
      if (memoryConstraint > adjustMemoryCost) 
      {
        FMBuilder::InferParametersGivenMemory(totalGenomeSize, alphabetSize, 
            memoryConstraint, fmBuilderParam) ; 
      }
      else
      {
        Utils::PrintLog("No enough memory, please increase the memory allocation.") ;
        exit(0) ;
      }
    }*/

    FMBuilder::Build(genomes, totalGenomeSize, alphabetSize, BWT, firstISA, fmBuilderParam) ;
    genomes.Free() ;
    Utils::PrintLog("Start to transform sampled SA to sequence ID.") ;
    TransformSampledSAToSeqId(fmBuilderParam, genomeSeqIds, genomeLens, totalGenomeSize) ;
    Utils::PrintLog("Start to compress BWT with RBBWT.") ;
    _fmIndex.Init(BWT, totalGenomeSize, 
        firstISA, fmBuilderParam, alphabetList, alphabetSize) ;
    Utils::PrintLog("centrifuger-build finishes.") ;
  }

  void OutputBuilderMeta(FILE *fp, const FMIndex<Sequence_RunBlock> &fm) 
  {
    fprintf(fp, "version\t" CENTRIFUGER_VERSION "\n") ;
    fprintf(fp, "SA_sample_rate\t%d\n", fm._auxData.sampleRate) ;

    time_t mytime = time(NULL) ;
    struct tm *localT = localtime( &mytime ) ;
    char stime[500] ;
    strftime( stime, sizeof( stime ), "%c", localT ) ;
    fprintf(fp, "build_date\t%s", stime) ;
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

    // .4.cfr file is tsv file for some version information
    sprintf(outputFileName, "%s.4.cfr", outputPrefix) ;
    fpOutput = fopen(outputFileName, "w") ;
    OutputBuilderMeta(fpOutput, _fmIndex) ;
    fclose(fpOutput) ;
  }
} ;

#endif

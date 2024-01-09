#ifndef _MOURISL_CLASSIFIER_HEADER
#define _MOURISL_CLASSIFIER_HEADER

#include <string.h>

#include "Taxonomy.hpp"
#include "compactds/FMIndex.hpp"
#include "compactds/Sequence_Hybrid.hpp"
#include "compactds/Sequence_RunBlock.hpp"
#include "compactds/SimpleVector.hpp"

using namespace compactds ;

//#define LI_DEBUG

struct _classifierParam 
{
  int maxResult ; // the number of entries in the results    
  int minHitLen ;
  int maxResultPerHitFactor ; // Get the SA/tax id for at most maxREsultPerHitsFactor * maxResult entries for each hit 
  _classifierParam()
  {
    maxResult = 1 ;
    minHitLen = 0 ;
    maxResultPerHitFactor = 40 ;
  }
} ;

// classification result for one read
struct _classifierResult
{
  size_t score ;
  size_t secondaryScore ;
  int hitLength ;
  int queryLength ;
  std::vector<std::string> seqStrNames ; // sequence names 
  std::vector<uint64_t> taxIds ; // taxonomy ids (original, not compacted)

  void Clear()
  {
    score = secondaryScore = 0 ;
    hitLength = queryLength = 0 ;
    seqStrNames.clear() ;
    taxIds.clear() ;
  }
} ;

// The hit result for each sequence ID
struct _seqHitRecord
{
  size_t seqId ;
  size_t score ;
  int hitLength ;
} ;

// Each individual hit on BWT string
struct _BWTHit
{
  size_t sp, ep ; //[sp, ep] range on BWT 
  int l ; // hit length
  int strand ; // -1: minus strand, 0: unkonwn, 1: plus strand
  int offset ; // 0-based offset to the end of the read (because the search is in backward fashion)
  
  _BWTHit(size_t isp, size_t iep, int il, int ioffset, int istrand)
  {
    sp = isp ;
    ep = iep ;
    l = il ;
    offset = ioffset ;
    strand = istrand ;
  }
} ;

class Classifier
{
private:
  FMIndex<Sequence_RunBlock> _fm ;
  Taxonomy _taxonomy ;
  std::map<size_t, size_t> _seqLength ;
  _classifierParam _param ;
  int _scoreHitLenAdjust ;
  char _compChar[256] ;
  
  void ReverseComplement(char *r, int len)
  {
    int i, j ;
    for (i = 0, j = len - 1 ; i < j ; ++i, --j)
    {
      char tmp ;
      tmp = r[i] ;
      r[i] = r[j] ;
      r[j] = tmp ;
    }
    for (i = 0 ; i < len ; ++i)
      r[i] = _compChar[(int)r[i]] ; 
  }

  void InferMinHitLen()
  {
    int mhl = 23 ; // Though centrifuge uses 22, but internally it filter length <= 22, so in our implementation, it should corresponds to 23.
    int alphabetSize = _fm.GetAlphabetSize() ; 
    uint64_t kmerspace = Utils::PowerInt(alphabetSize, mhl)/ 2 ;
    uint64_t n = _fm.GetSize() ;
    for ( ; mhl <= 32 ; ++mhl)
    {
      if (kmerspace >= 100 * n)
        break ;
      kmerspace *= alphabetSize ;
    }
    _param.minHitLen = mhl ;
  }

  //l: hit length
  size_t CalculateHitScore(int l)
  {
    if (l < _param.minHitLen)
      return 0 ;
    return (l - _scoreHitLenAdjust) * (l - _scoreHitLenAdjust) ;
  }

  // one hit
  size_t CalculateHitScore(const struct _BWTHit &hit)
  {
    return CalculateHitScore(hit.l) ;
  }

  // hit list
  size_t CalculateHitsScore(const SimpleVector<struct _BWTHit> &hits)
  {
    int i ;
    int hitCnt = hits.Size() ;
    size_t score = 0 ;
    for (i = 0 ; i < hitCnt ; ++i)
    {
      score += CalculateHitScore(hits[i]) ; 
    }
    return score ;
  }

  //@return: the number of hits 
  size_t GetHitsFromRead(char *r, size_t len, SimpleVector<struct _BWTHit> &hits) 
  {
    size_t sp = 0, ep = 0 ;
    int l = 0 ;
    int remaining = len ;
    
    while (remaining >= _param.minHitLen)
    {
      l = _fm.BackwardSearch(r, remaining, sp, ep) ;
      if (l >= _param.minHitLen && sp <= ep)
      {
        struct _BWTHit nh(sp, ep, l, len - remaining, 0) ;
        hits.PushBack(nh) ;
      }

      // +1 is to skip the base
      remaining -= (l + 1) ;
    }
    return hits.Size() ;
  }

  // The hit search method has strand bias, so we shall use the other strand
  //   information to mitigate the bias. This is important if some strain's 
  //   sequence is reverse-complemented.
  //  
  // E.g.: 100bp read: 90bp-N-9bp from two strains, one forward, one rc 
  //   Forward search probably would be ~20bp random hits + ~80 real hit
  //   Reverse-complement search: will be 90bp real hit
  //   As a result, we will lose the forward candidate
  void AdjustHitBoundaryFromStrandHits(char *r, char *rc, int len, 
      SimpleVector<struct _BWTHit> *strandHits)
  {
    int i, j ;
    if (!strandHits[0].Size() || !strandHits[1].Size())
      return ;
    int hitSize[2] = {strandHits[0].Size(), strandHits[1].Size()} ;
  
    size_t sp, ep ;
    int l ;
    j = hitSize[0] - 1 ;
    for (i = 0 ; i < hitSize[1] ; ++i)
    {
      int left, right ; // range on the read, original read
      right = len - strandHits[1][i].offset - 1 ; 
      left = right - strandHits[1][i].l + 1 ;
      for ( ; j >= 0 ; --j)
      {
        int rcLeft, rcRight ;
        rcLeft = strandHits[0][j].offset ;
        rcRight = rcLeft + strandHits[0][j].l - 1 ;
        
        if (rcLeft >= right) // no overlap yet 
          continue ;
        if (left >= rcRight) // already passed
          break ;
        if (left == rcLeft && right == rcRight) // both hits are good
          break ;
        if (left < rcLeft && rcRight < right) // forward hit contains reverse hit
          break ;
        if (rcLeft < left && right < rcRight) // reverse hit contains forward hit
          break ;
        if (rcRight > right)
        {
          l = _fm.BackwardSearch(r, rcRight + 1, sp, ep) ;
          if (rcRight - l + 1 == left && sp <= ep)
          {
            struct _BWTHit nh(sp, ep, l, len - rcRight + 1, 1) ;
            strandHits[1][i] = nh ;
          }
        }

        if (left < rcLeft)
        {
          l = _fm.BackwardSearch(rc, len - left, sp, ep) ;
          if (left + l - 1 == rcRight && sp <= ep)
          {
            struct _BWTHit nh(sp, ep, l, left, -1) ;
            strandHits[0][j] = nh ;
          }
        }
      }
    }
  }

  // It seems the performance for not synchronize mate pair direction works better
  size_t SearchForwardAndReverseWithWeakMateDirection(char *r1, char *r2, SimpleVector<struct _BWTHit> &hits)
  {
    int i, k, ridx ;
    
    hits.Clear() ;
    SimpleVector<struct _BWTHit> strandHits[2] ; // 0: minus strand, 1: postive strand
    
    for (ridx = 0 ; ridx <= 1 ; ++ridx) //0-r1, 1-r2
    {
      if (ridx == 1 && r2 == NULL)
        break ;

      char *r = r1 ;
      if (ridx == 1)
        r = r2 ;
      char *rc = NULL ;

      int rlen = strlen(r) ;
      rc = strdup(r) ;
      ReverseComplement(rc, rlen) ;

      strandHits[0].Clear() ; 
      strandHits[1].Clear() ;
      //Notice that GetHitsFromRead will not clear the hits
      GetHitsFromRead(r, rlen, strandHits[1]) ;
      GetHitsFromRead(rc, rlen, strandHits[0]) ;
      AdjustHitBoundaryFromStrandHits(r, rc, rlen, strandHits) ;
      
      size_t strandScore[2] ;
      //int strandLongestHit[2] = {0, 0} ;
      for (k = 0 ; k < 2 ; ++k)
      {
        int size = strandHits[k].Size() ;
        for (i = 0 ; i < size ; ++i)
        {
          strandHits[k][i].strand = (2 * k - 1) * (ridx == 0 ? 1 : -1) ; // the strand is with respect to the template, not read
          //if (strandHits[k][i].l > strandLongestHit[k])
          //  strandLongestHit[k] = strandHits[k][i].l ;
        }
        strandScore[k] = CalculateHitsScore(strandHits[k]) ;
      }
      
      if (strandScore[1] >= strandScore[0])
        hits.PushBack(strandHits[1]) ;
      if (strandScore[0] >= strandScore[1]) // if equal, both strands will be added
        hits.PushBack(strandHits[0]) ;
      /*else
      {
        if (strandLongestHit[1] >= strandLongestHit[0])
          hits.PushBack(strandHits[1]) ;
        if (strandLongestHit[0] >= strandLongestHit[1]) 
          hits.PushBack(strandHits[0]) ;
      }*/

      free(rc) ;
    }
    return hits.Size() ;
  }

  //@return: the size of the hits after selecting the strand 
  size_t SearchForwardAndReverse(char *r1, char *r2, SimpleVector<struct _BWTHit> &hits)
  {
    int i, k ;
    char *rcR1 = NULL ;
    char *rcR2 = NULL ;
    int r1len = strlen(r1) ;
    rcR1 = strdup(r1) ;
    ReverseComplement(rcR1, r1len) ;
    
    SimpleVector<struct _BWTHit> strandHits[2] ; // 0: minus strand, 1: postive strand
    
    GetHitsFromRead(r1, r1len, strandHits[1]) ;
    GetHitsFromRead(rcR1, r1len, strandHits[0]) ;
    AdjustHitBoundaryFromStrandHits(r1, rcR1, r1len, strandHits) ;
    if (r2)
    {
      rcR2 = strdup(r2) ;
      int r2len = strlen(r2) ;
      ReverseComplement(rcR2, r2len) ;
      SimpleVector<struct _BWTHit> r2StrandHits[2] ; // 0: minus strand, 1: postive strand
      
      GetHitsFromRead(r2, r2len, r2StrandHits[1]) ;
      GetHitsFromRead(rcR2, r2len, r2StrandHits[0]) ;
      AdjustHitBoundaryFromStrandHits(r2, rcR2, r2len, r2StrandHits) ;
      for (i = 0 ; i <= 1 ; ++i)
        strandHits[i].PushBack(r2StrandHits[1 - i]) ;
    }
    
    size_t strandScore[2] ;
    for (k = 0 ; k < 2 ; ++k)
    {
      int size = strandHits[k].Size() ;
      for (i = 0 ; i < size ; ++i)
        strandHits[k][i].strand = 2 * k - 1 ; // the strand is with respect to the template, not read
      strandScore[k] = CalculateHitsScore(strandHits[k]) ;
    }
   
    if (strandScore[1] > strandScore[0])
      hits = strandHits[1] ;
    else if (strandScore[0] > strandScore[1])
      hits = strandHits[0] ;
    else
    {
      hits = strandHits[1] ;
      hits.PushBack(strandHits[0]) ;
    }
    
    free(rcR1) ;
    if (rcR2)
      free(rcR2) ;

    return hits.Size() ;
  }

  size_t GetClassificationFromHits(const SimpleVector<struct _BWTHit> &hits, struct _classifierResult &result)
  {
    int i, k ;
    size_t j ;
    int hitCnt = hits.Size() ;
    std::map<size_t, struct _seqHitRecord > seqIdStrandHitRecord[2] ;
    
    struct _seqHitRecord prevUniqHitRecord ; // record information from previous unique hit 
    prevUniqHitRecord.seqId = 0 ;
    prevUniqHitRecord.hitLength = 0 ;
    prevUniqHitRecord.score = 0 ;

    bool mixStrand = false ;
    for (i = 1 ; i < hitCnt ; ++i)
    {
      if (hits[i].strand != hits[i - 1].strand)
      {
        mixStrand = true ;
        break ;
      }
    }

    // The hit seqId need to consider the strand separately.
    //   Because sometimes a read can hit both the plus and minus strand and will artifically double the hit length.
    for (i = 0 ; i < hitCnt ; ++i)
    {
      size_t score = CalculateHitScore(hits[i]) ;
      std::map<size_t, int> localSeqIdHit ;
      k = (hits[i].strand + 1) / 2 ;
#ifdef LI_DEBUG
      printf("hit: %d sp-ep: %lu %lu %lu offset_l: %d %d\n", i, hits[i].sp, hits[i].ep, hits[i].ep - hits[i].sp + 1, hits[i].offset, hits[i].l) ;
#endif
      const size_t maxEntries = _param.maxResult * _param.maxResultPerHitFactor ;
      if (hits[i].ep - hits[i].sp + 1 <= maxEntries 
          || _param.maxResultPerHitFactor <= 0)
      {
        for (j = hits[i].sp ; j <= hits[i].ep ; ++j)
        {
          size_t backsearchL = 0 ;
          size_t seqId = _fm.BackwardToSampledSA(j, backsearchL) ;
#ifdef LI_DEBUG
          printf("%lu\n", _taxonomy.GetOrigTaxId( _taxonomy.SeqIdToTaxId(seqId) )) ;
#endif
          localSeqIdHit[seqId] = 1 ;
        }
      }
      else
      {
        // Since the first entry and last entry are likely to be more different
        //   taxonomy-wisely, we shall search "bidirectionally" to make sure 
        //   both end is covered
        size_t rangeSize = hits[i].ep - hits[i].sp + 1 ;
        size_t step = DIV_CEIL(rangeSize, maxEntries) ;
        size_t resolvedCnt = 0 ;
        for (j = hits[i].sp ; j <= hits[i].ep ; j += step)
        {
          size_t backsearchL = 0 ;
          size_t seqId = _fm.BackwardToSampledSA(j, backsearchL) ;
#ifdef LI_DEBUG
          printf("%lu\n", _taxonomy.GetOrigTaxId( _taxonomy.SeqIdToTaxId(seqId) )) ;
#endif
          localSeqIdHit[seqId] = 1 ;
          ++resolvedCnt ;
        }

        for (j = hits[i].ep ; j >= hits[i].sp && j <= hits[i].ep ; j -= step)
        {
          size_t backsearchL = 0 ;
          size_t seqId = _fm.BackwardToSampledSA(j, backsearchL) ;
#ifdef LI_DEBUG
          printf("%lu\n", _taxonomy.GetOrigTaxId( _taxonomy.SeqIdToTaxId(seqId) )) ;
#endif
          localSeqIdHit[seqId] = 1 ;
          ++resolvedCnt ;
          if (resolvedCnt >= maxEntries)
            break ;
        }
      }

      // Update the scores for each seqid
      for (std::map<size_t, int>::iterator iter = localSeqIdHit.begin() ;
          iter != localSeqIdHit.end() ; ++iter)
      {
        size_t seqId = iter->first ;
        if (!mixStrand && i > 0 && hits[i].ep == hits[i].sp && 
            hits[i - 1].ep == hits[i - 1].sp && 
            hits[i - 1].strand == hits[i].strand &&
            hits[i - 1].offset + hits[i - 1].l + 1 == hits[i].offset && // the other strand adjustication may cause overlaps of the hit regions. Make sure the two hits only separate by 1 base.
            seqId == prevUniqHitRecord.seqId) // Merge adjacent unique hits
        {
          seqIdStrandHitRecord[k][seqId].score -= prevUniqHitRecord.score ;

          prevUniqHitRecord.hitLength += hits[i].l ;
          prevUniqHitRecord.score = CalculateHitScore(prevUniqHitRecord.hitLength) ;
          seqIdStrandHitRecord[k][seqId].score += prevUniqHitRecord.score ;
          seqIdStrandHitRecord[k][seqId].hitLength += hits[i].l ;
        }
        else // Regularly update the score
        {
          if (seqIdStrandHitRecord[k].find(seqId) == seqIdStrandHitRecord[k].end())
          {
            seqIdStrandHitRecord[k][seqId].seqId = seqId ;
            seqIdStrandHitRecord[k][seqId].score = score ;
            seqIdStrandHitRecord[k][seqId].hitLength = hits[i].l ;
          }
          else
          {
            seqIdStrandHitRecord[k][seqId].score += score ;
            seqIdStrandHitRecord[k][seqId].hitLength += hits[i].l ;
          }
        
          if (hits[i].ep == hits[i].sp)
          {
            prevUniqHitRecord.seqId = seqId ;
            prevUniqHitRecord.score = score ;
            prevUniqHitRecord.hitLength = hits[i].l ;
          }
        }
      }
    }

    // Select the best score
    size_t bestScore = 0 ;
    size_t secondBestScore = 0 ;
    size_t bestScoreHitLength = 0 ;
    for (k = 0 ; k <= 1 ; ++k)
    {
      for (std::map<size_t, struct _seqHitRecord>::iterator iter = seqIdStrandHitRecord[k].begin() ; 
          iter != seqIdStrandHitRecord[k].end() ; ++iter)
      {
#ifdef LI_DEBUG
        printf("score: %lu %lu %d\n", _taxonomy.GetOrigTaxId( _taxonomy.SeqIdToTaxId(iter->first)), iter->second.score, iter->second.hitLength) ;
#endif
        if (iter->second.score > bestScore)
        {
          secondBestScore = bestScore ;
          bestScore = iter->second.score ;
          bestScoreHitLength = iter->second.hitLength ;
        }
        else if (iter->second.score > secondBestScore)
          secondBestScore = iter->second.score ;
      }
    }

    // Collect match corresponding to the best score.
    result.score = bestScore ;
    result.secondaryScore = secondBestScore ;
    result.hitLength = bestScoreHitLength ;

    SimpleVector<size_t> bestSeqIds ;
    std::map<size_t, int> bestSeqIdUsed ;
    for (k = 0 ; k <= 1 ; ++k)
    {
      for (std::map<size_t, struct _seqHitRecord>::iterator iter = seqIdStrandHitRecord[k].begin() ; 
          iter != seqIdStrandHitRecord[k].end() ; ++iter)
      {
        if (iter->second.score == bestScore && 
            bestSeqIdUsed.find(iter->first) == bestSeqIdUsed.end())
        {
          bestSeqIds.PushBack(iter->first) ;
          bestSeqIdUsed[iter->first] = 1 ;
        }
      }
    }

    if (bestSeqIds.Size() > 1)
      result.secondaryScore = bestScore ;

    if (bestSeqIds.Size() <= _param.maxResult)
    {
      int size = bestSeqIds.Size() ;
      for (i = 0 ; i < size ; ++i)
      {
        result.seqStrNames.push_back( _taxonomy.SeqIdToName(bestSeqIds[i]) ) ;
        result.taxIds.push_back( _taxonomy.GetOrigTaxId(_taxonomy.SeqIdToTaxId( bestSeqIds[i] )) ) ;
      }
    }
    else
    {
      int size = bestSeqIds.Size() ;
      SimpleVector<size_t> bestSeqTaxIds ;
      bestSeqTaxIds.Reserve(size) ;
      for (i = 0 ; i < size ; ++i)
        bestSeqTaxIds.PushBack( _taxonomy.SeqIdToTaxId(bestSeqIds[i]) ) ;

      SimpleVector<size_t> taxIds ;
      _taxonomy.ReduceTaxIds(bestSeqTaxIds, taxIds, _param.maxResult) ;

      size = taxIds.Size() ;
      for (i = 0 ; i < size ; ++i)
      {
        std::string rankName(_taxonomy.GetTaxRankString( _taxonomy.GetTaxIdRank(taxIds[i])) ) ;
        result.seqStrNames.push_back( rankName ) ;
        result.taxIds.push_back( _taxonomy.GetOrigTaxId(taxIds[i]) ) ;
      }
    }
    return result.taxIds.size() ;
  }
  
public:
  Classifier() 
  {
    _scoreHitLenAdjust = 15 ;
    int i ;
    for (i = 0 ; i < 256 ; ++i)
      _compChar[i] = 'N' ;
    _compChar['A'] = 'T' ;
    _compChar['C'] = 'G' ;
    _compChar['G'] = 'C' ;
    _compChar['T'] = 'A' ;
  }

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
 
    // .1.cfr file for FM index
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
    
    Utils::PrintLog("Finishes loading index.") ;
    
    _param = param ;
    if (_param.minHitLen <= 0)
    {
      InferMinHitLen() ;
      Utils::PrintLog("Inferred --min-hitlen: %d", _param.minHitLen) ;
    }

    free(nameBuffer) ;
  }

  // Main function to return the classification results 
  void Query(char *r1, char *r2, struct _classifierResult &result)
  {
    result.Clear() ;

    SimpleVector<struct _BWTHit> hits ;
    
    SearchForwardAndReverse(r1, r2, hits) ;
    GetClassificationFromHits(hits, result) ;
    result.queryLength = strlen(r1) ;
    if (r2)
      result.queryLength += strlen(r1) ;
  }

  const Taxonomy &GetTaxonomy()
  {
    return _taxonomy ;
  }
} ;

#endif 

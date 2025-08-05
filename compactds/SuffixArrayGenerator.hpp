#ifndef _MOURISL_COMPACTDS_SUFFIXARRAY_GENERATOR
#define _MOURISL_COMPACTDS_SUFFIXARRAY_GENERATOR

#include <vector>

#include "FixedSizeElemArray.hpp"
#include "DifferenceCover.hpp"

// The class handle the generation of suffix array by chunks
// The chunk creation is based the sampled difference cover (Algorithm 11.9 from the textbook is commented out)
namespace compactds {
class SuffixArrayGenerator
{
private:
  size_t _n ;
  size_t _space ;
  size_t _alphabetSize ;
  
  // The variables relate to generate the boundaries/_cuts 
  size_t _b ;
  size_t* _cuts ; // The index on T  
  size_t _cutCnt ;
  size_t **_cutLCP ; // self LCP for each cut

  // The variables relate to the difference cover and its ISAs
  DifferenceCover _dc ;
  size_t *_dcISA ; // The difference cover's index should be compacted when query this ISA 
  size_t _dcSize ; 
  
  // Relate to cut ============================================  
#if 0 // The commented out codes is for Algorithm 11.9, which might be too slow for very repetitive sequence (i.e: ACGTACGTACGT....), so we have another implementation now

  // The functions relate to generate cut
  int SuffixCompareCutString(const FixedSizeElemArray &T, size_t n, size_t i, const FixedSizeElemArray &s, size_t k) 
  {
    if (k == 0)
      return 0 ;
    else
      return T.SubrangeCompare(i, i + k - 1, s, 0, k - 1) ;
  }
 
  // Count the size for each alphabet following current cut string.
  // s: the cut string
  // k: length of the cut string
  void CountCutExtension(const FixedSizeElemArray &T, size_t n, const FixedSizeElemArray &s, size_t k, size_t *alphabetCounts)
  {
    size_t i ;
    memset(alphabetCounts, 0, sizeof(alphabetCounts[0]) * _alphabetSize) ;
    for (i = 0 ; i < n - k ; i += downsample)
    {
      if (!SuffixCompareCutString(T, n, i, s, k))
        alphabetCounts[ T.Read(i + k) ] += downsample ;
    }
  }
  
  // s: the cut string
  // k: length of the cut string
  void ExpandInCut(const FixedSizeElemArray &T, size_t n, FixedSizeElemArray &s, size_t k, size_t *chunkLens)
  {
    size_t* alphabetCounts = (size_t *)malloc(sizeof(size_t) * _alphabetSize) ;
    CountCutExtension(T, n, s, k, alphabetCounts) ;
    size_t c ;
    for (c = 0 ; c < _alphabetSize ; ++c)
    {
      s.Write(k, c) ; 
      if (alphabetCounts[c] <= b)
      {
        if (chunkLens[_cutCnt - 1] + alphabetCounts[c] > b)
        {
          _cuts[_cutCnt].InitFromOtherPrefix(s, k + 1) ;
          ++_cutCnt ;
          chunkLens[_cutCnt - 1] = 0 ;
        }
        chunkLens[_cutCnt - 1] += alphabetCounts[c] ;
      }
      else
        ExpandInCut(T, n, s, k + 1, chunkLens) ;
    }
    free(alphabetCounts) ;
  }
  
  size_t GenerateCuts(const FixedSizeElemArray &T, size_t n)
  {
    size_t m = DIV_CEIL(2 * n, b) ;
    _cuts = new FixedSizeElemArray[m + 1] ;
    size_t* chunkLens = (size_t *)malloc(sizeof(size_t) * (m + 1)) ;
    FixedSizeElemArray s ;
    s.Malloc(T.GetElemLength(), _n < _b ? _n : b) ;
   
    _cuts[0].Malloc(T.GetElemLength(), 0) ;
    chunkLens[0] = 0 ;
    _cutCnt = 1 ;
    ExpandInCut(T, n, s, 0, chunkLens) ;
    _cuts[_cutCnt].Malloc(T.GetElemLength(), 0) ;
    printf("%d\n", _cuts[1].GetSize()) ;
    free(chunkLens) ;
    return _cutCnt ;
  }
#endif

  // TODO: handle the case where we don't generate difference cover
  size_t GenerateCuts(size_t *_dcSA)
  {
    size_t i ;
    size_t blockCnt = DIV_CEIL(_n, _b) ;
    size_t stride = DIV_CEIL(_dcSize, blockCnt) ;
    blockCnt = DIV_CEIL(_dcSize, stride) ;
    _cutCnt = blockCnt ;
    _cuts = (size_t *)malloc(sizeof(size_t) * (_cutCnt + 1)) ;
    _cuts[0] = _n ;
    for (i = stride ; i < _dcSize ; i += stride)
    {
      //printf("_cuts %d = %d\n", i / stride, _dcSA[i]) ;
      _cuts[i / stride] = _dcSA[i] ;
    }
    _cuts[blockCnt] = _n ;
    return _cutCnt ;
  }

  // For each cut s, compute LCP(s, s[i:]) for i <= maxSize
  void ComputeCutLCP(const FixedSizeElemArray &T, size_t n, size_t maxSize) 
  {
    size_t i, j, l ;
    _cutLCP = (size_t **)malloc(sizeof(*_cutLCP) * _cutCnt) ;
    for (i = 0 ; i < _cutCnt ; ++i)
    {
      size_t jopenend = n - _cuts[i] ;
      if (jopenend > maxSize )
        jopenend = maxSize ;
      _cutLCP[i] = (size_t *)malloc(sizeof(_cutLCP[i][0]) * jopenend) ;
      if (jopenend == 0)
        continue ;
      _cutLCP[i][0] = jopenend ;
      for (j = 1 ; j < jopenend ; ++j) 
      {
        for (l = 0 ; l < jopenend - j ; ++l)
        {
          if (T.Read(_cuts[i] + l) != T.Read(_cuts[i] + j + l))
            break ;
        }
        _cutLCP[i][j] = l ;
      }
    }
  }
  
  // Compare the T[i,...] with a cut ci, and adjust other auxiliary data relating
  // rightmosti: the start position corresponding to the rightmost j 
  // rightj: the rightmost end position
  //   The auxiliary data is for reusing some information from the 
  //   preiouv T[i-1,...] comparison
  // return: sign represent T[i:n]-cut . 
  int CompareCutUsingCutLCP(const FixedSizeElemArray &T, size_t n, size_t i, size_t ci, 
      size_t &rightmosti, size_t &rightmostj)
  {
    size_t j ; // position on T
    size_t cpos ; // position on cut
    size_t overlap = 0 ;
    size_t cutLen = n - _cuts[ci] ;
    if (i == _cuts[ci])
      return 0 ; // return 0 for equal
    //if (i == 15 && ci == 6)
    //  printf("> %d %d %d\n", i, rightmosti, rightmostj) ;

    if (cutLen > (size_t)_dc.GetV())
      cutLen = _dc.GetV() ;
    if (rightmostj > 0 && i <= rightmostj)  // i<=rightmostj? 
    {
      overlap = _cutLCP[ci][i - rightmosti] ;  
      if (rightmostj <= i + overlap - 1) // we may need to update the range 
      {
        /*for (j = rightmostj + 1, cpos = rightmostj - i + 1 ; 
            j < n && cpos < cutLen; ++j, ++cpos)
        {
          if (T.Read(j) != T.Read(_cuts[ci] + cpos))
            break ;
        }
        --j ; --cpos ; 
        rightmostj = j ;
        rightmosti = i ;*/
        size_t localMatchCnt = T.PrefixMatchLen(rightmostj + 1, n - 1, 
            T, _cuts[ci] + rightmostj - i + 1, _cuts[ci] + cutLen - 1) ;
        rightmosti = i ;
        j = (rightmostj + 1) + localMatchCnt - 1 ;
        rightmostj = j ;
        cpos = j - i ;
      }
      else
      {
        return (int)T.Read(i + overlap) - T.Read(_cuts[ci] + overlap) ;
      }
    }
    else
    {
      /*for (j = i, cpos = 0 ; j < n && cpos < cutLen ; ++j, ++cpos)
        if (T.Read(j) != T.Read(_cuts[ci] + cpos))
          break ;
      //printf("%d %d %d. %d\n", i, ci, _cuts[ci], cpos) ;
      if (cpos == 0)
        return T.Read(j) - T.Read(_cuts[ci] + cpos) ;

      --j; --cpos ;
      rightmostj = j ;
      rightmosti = i ;*/

      size_t localMatchCnt = T.PrefixMatchLen(i, n - 1, 
          T, _cuts[ci], _cuts[ci] + cutLen - 1) ;
      if (localMatchCnt == 0)
        return T.Read(i) - T.Read(_cuts[ci]) ;
      j = i + localMatchCnt - 1 ;
      rightmostj = j ;
      cpos = localMatchCnt - 1 ;
      rightmosti = i ;
    }
    
    if (j == n - 1 && cpos == n - _cuts[ci] - 1)
      return 0 ;
    else if (j == n - 1)
      return -1 ;
    else if (cpos == n - _cuts[ci] - 1)
      return 1 ;
    else
    {
      int cmp = T.Read(j + 1) - T.Read(_cuts[ci] + cpos + 1) ;
      if (cmp == 0)
        return CompareSuffixWithDC(i, _cuts[ci], n) ;
      else
        return cmp ;
    }
  }

  // Relate to suffix sorting ============================================  
  void Swap(size_t &a, size_t &b)
  {
    size_t tmp ;
    tmp = a ; a = b ; b = tmp ;
  }

  // Compare T[a:] and T[b:] directly with DC, which assumes their first
  //  v prefix are matched
  // @return: sign(T[a:]-T[b:])
  int CompareSuffixWithDC(size_t a, size_t b, size_t n)
  {    
    if (a == b)
      return 0 ;
    int delta = _dc.Delta(a, b) ;
    int compare = 0 ;
    if (a + delta >= n)
      compare = 1 ;
    else if (b + delta >= n)
      compare = -1 ;
    else
    {
      size_t aisa = _dcISA[ _dc.CompactIndex(a + delta) ] ;
      size_t bisa = _dcISA[ _dc.CompactIndex(b + delta) ] ;
      if (aisa < bisa)
        compare = -1 ;
      else 
        compare = 1 ;
    }
    return compare ;
  }

  // Use difference cover to quick sort the SA. We don't need to pass T now  
  void QSortWithDC(size_t *sa, size_t m, size_t s, size_t e, size_t n)
  {
    if (s >= e)
      return ;
    // Partition
    Swap(sa[s], sa[(s + e)/2]) ;
    size_t pivot = sa[s] ; // pivot is the median element
    size_t pi, pj ; // partiation indexes
    pi = s + 1; // pi points to the current process element
    pj = e + 1 ; // pj points the first element of the second chunk
    while (pi < pj)
    {
      int comparePivot = CompareSuffixWithDC(sa[pi], pivot, n) ;

      if (comparePivot < 0)
        ++pi ;
      else
      {
        Swap(sa[pi], sa[pj - 1]) ;
        --pj ;
      }
    }
    Swap(sa[s], sa[pi - 1]) ;
    if (pi > 2)
      QSortWithDC(sa, m, s, pi - 2, n) ;
    QSortWithDC(sa, m, pi, e, n) ;
  }

  // Compare two SA's in the same h-group
  // sai: SA[i], saj: SA[j]
  // return: sign(T[sai..] - T[saj..]), i.e. sign(rank[sa[i+h]] - rank[sa[j+h]])
  int CompareMultikeyQSortForLSandDC(size_t *rank, size_t n, size_t sai, size_t saj, size_t h)
  {
    if (sai == saj)
      return 0 ;
    if (sai + h >= _n && saj + h >= _dcSize)
      return (sai + h > saj + h) ? -1 : 1 ; 
    else if (sai + h >= n && saj + h < n)
      return -1 ;
    else if (sai + h < n && saj + h >= n)
      return 1 ;
    else
    {
      size_t rih = rank[ _dc.CompactIndex(sai + h) ] ;
      size_t rjh = rank[ _dc.CompactIndex(saj + h) ] ;
      if (rih == rjh)
        return 0 ;
      else
        return rih > rjh ? 1 : -1 ; 
    }
  }

  // The internal same-rank q sort for the Larsson-Sadakane SA sorting's h groupd
  // Sort SA[s..e] inclusive, using their +h as the value to sort
  // This sort NO need to be stable! The entries with the same value will have the same rank.
  void MultikeyQSortForLSandDC(size_t *sa, size_t *rank, size_t n, size_t s, size_t e, size_t h, size_t *newRank)
  {
    if (s > e)
      return ;
    if (s == e)
    {
      newRank[ _dc.CompactIndex(sa[s]) ] = s ;
      return ;
    }

    size_t i ;
    /*if (e - s + 1 <= 5) // Use select sort for short intervals, but seems not useful in my implementation
    {
      size_t j ;
      size_t r = e ; // track the rank
      for (i = e ; i >= s && i <= e ; --i) // for right to left so the update of newRank is easier
      {
        for (j = s ; j < i ; ++j)
          if (CompareMultikeyQSortForLSandDC(rank, n, sa[i], sa[j], h) < 0)
            Swap(sa[i], sa[j]) ;

        if (i < e 
            && CompareMultikeyQSortForLSandDC(rank, n, sa[i], sa[i + 1], h) == 0)
          newRank[ _dc.CompactIndex(sa[i]) ] = r ;
        else
        {
          newRank[ _dc.CompactIndex(sa[i]) ] = i ;
          r = i ;
        }
      }
      return ;
    }*/

    // Partition
    size_t pivot = sa[(s+e)/2] ;
    size_t pi, pj, pk ; // partiation indexes
    pi = s ; // pi points to the first element of the middle chunk
    pj = s ; // pj is the current element
    pk = e ; // pk points the last element of the middle chunk
    while (pj <= pk)
    {
      int comparePivot = CompareMultikeyQSortForLSandDC(rank, n, sa[pj], pivot, h) ;

      if (comparePivot == -1) // less than pivot
      {
        Swap(sa[pi], sa[pj]) ;
        ++pi ; ++pj ;
      }
      else if (comparePivot == 1)
      {
        Swap(sa[pj], sa[pk]) ;
        if (pk == 0)
          break ;
        --pk ;
      }
      else
        ++pj ;
    }
    //printf("===%d %d. %d %d %d\n", s, e, pi, pj, pk) ;

    if (pi >= 1)
      MultikeyQSortForLSandDC(sa, rank, n, s, pi - 1, h, newRank) ;
    for (i = pi ; i <= pk ; ++i)
      newRank[ _dc.CompactIndex(sa[i]) ] = pk ;
    MultikeyQSortForLSandDC(sa, rank, n, pj, e, h, newRank) ;
  } 

  // Sort T, and only consider the positions in sa[s..e], and total size of sa is m 
  // s, e: the range for sa
  // d: the preifx already matched in T[s..e], kind of as depth.
  // dcStrategy: how to use the difference cover. 0-no _dc, 1-use _dc, 2-return when reach _dcv
  void MultikeyQSort(const FixedSizeElemArray &T, size_t n, size_t *sa, size_t m, size_t s, size_t e, size_t d, int dcStrategy, size_t *alphabetCounts)
  {
    if (s >= e)
      return ;
    if (dcStrategy != 0 && d >= (size_t)_dc.GetV())
    {
      if (dcStrategy == 2)
        return ;
      else if (dcStrategy == 1)
      {
        // We now use compare the suffix using difference cover
        QSortWithDC(sa, m, s, e, n) ;
        return ;
      }
    }

    size_t i, j ;
    size_t tmp ;
    
    size_t shortRangeDTarget = (size_t)_dc.GetV() ; // Set a large value, so it is ineffective when the range is not short.
    if (dcStrategy == 1 && s + 3 >= e) // very few element (4 for now), then we check whether the difference cover is effective, if so, we can use insert sort
    {
      size_t maxDelta = 0 ;
      for (i = s ; i < e ; ++i)
      {
        for (j = i + 1 ; j <= e ; ++j)
        {
          tmp = _dc.Delta(sa[i], sa[j]) ;
          if (tmp > maxDelta) 
            maxDelta = tmp ;
        }
      }
      
      shortRangeDTarget = maxDelta ;
    }

    // Find pivot
    size_t pivot = 0 ;
    
    // quick check whether every suffix is the same using blocks
    //  Now it will also reach to the point where the differnece happens
    const int alphabetBits = T.GetElemLength() ;
    const int block = WORDBITS / alphabetBits ;
    while (1)
    {
      if ((dcStrategy != 0 && d >= (size_t)_dc.GetV())
          || (dcStrategy == 1 && s + 3 >= e && d >= shortRangeDTarget))
        break ;
      bool passEnd = false ; // any suffix pass the end of the T
      WORD foundw = 0 ;
      if (sa[s] + d + block - 1 < n)
        foundw = T.PackRead(sa[s] + d, block) ;
      else
        passEnd = true ;
      
      int skipSize = block ;
      if (!passEnd)
      {
        for (i = s + 1 ; i <= e ; ++i)
        {
          if (sa[i] + d + block - 1 < n)
          {
            WORD w = T.PackRead(sa[i] + d, block) ;
            if (w != foundw)
            {
              int prefixLen = Utils::CountTrailingZeros(w ^ foundw) / alphabetBits;
              if (prefixLen < skipSize)
              {
                skipSize = prefixLen ;
                if (skipSize == 0) // No shared prefix at all for this range, so we don't need to check future elements
                  break ;
              }
            }
          }
          else
          {
            //if (passEnd == false) // why I allow the first one to pass the end in the original implementation?
            passEnd = true ;
            break ;
          }
        }
      }
      
      if (!passEnd && skipSize >= block)
        d += block ;
      else
      {
        if (!passEnd)
          d += skipSize ;
        break ;
      }
    }

    // We can use insert sort if the shared prefix (d) is long enough for the short range case
    // If the range already satisfy the d target, the while loop above
    // actually won't happen. So we don't need this insert sort before
    // the loop, with the same speed.
    if (dcStrategy == 1 && s + 3 >= e
        && d >= shortRangeDTarget)
    {
      for (i = s ; i < e ; ++i)
      {
        size_t minTag = i ;
        for (j = i + 1 ; j <= e ; ++j)
        {
          if (CompareSuffixWithDC(sa[j], sa[minTag], n) < 0)
            minTag = j ;
        }

        Swap(sa[i], sa[minTag]) ;
      }
      return ;    
    }

    // Real search
    while (1)
    {
      if (dcStrategy != 0 && d >= (size_t)_dc.GetV())
        break ;
      
      memset(alphabetCounts, 0, sizeof(alphabetCounts[0]) * (_alphabetSize + 1)) ;
      for (i = s ; i <= e ; ++i)
      {
        if (sa[i] + d < n)
          ++alphabetCounts[T.Read(sa[i] + d) + 1] ;
        else
          ++alphabetCounts[0] ;
      }
      if (alphabetCounts[0] == e - s + 1)
        return ;
      tmp = 0 ;
      for (i = 0 ; i <= _alphabetSize ; ++i)
      {
        if (alphabetCounts[i] > 0)
          ++tmp ;
      }
      if (tmp == 1) // the next character is the same for all the suffixes in the range
      {
        ++d ;
        continue ;
      }

      tmp = 0 ;
      for (i = 0 ; i <= _alphabetSize ; ++i)
      {
        tmp += alphabetCounts[i] ;
        if (tmp >= (e - s + 1) / 2)
          break ;
      }
      pivot = i ; // pivot character
      if (pivot > 0)
        --pivot ;
      break ;
    }
    
    // Partition
    size_t pi, pj, pk ; // partiation indexes
    pi = s ; // pi points to the first element of the middle chunk
    pj = s ; // pj is the current element
    pk = e ; // pk points the last element of the middle chunk
    while (pj <= pk)
    {
      int comparePivot = 0 ;
      if (sa[pj] + d >= n)
        comparePivot = -1 ;
      else
      {
        size_t c = T.Read(sa[pj] + d) ;
        if (c < pivot)
          comparePivot = -1 ;
        else if (c == pivot)
          comparePivot = 0 ;
        else
          comparePivot = 1 ;
      }

      if (comparePivot == -1)
      {
        Swap(sa[pi], sa[pj]) ;
        ++pi ; ++pj ;
      }
      else if (comparePivot == 1)
      {
        Swap(sa[pj], sa[pk]) ;
        if (pk == 0)
          break ;
        --pk ;
      }
      else
        ++pj ;
    }
    
    // Recursive sorting
    if (pi >= 1)
      MultikeyQSort(T, n, sa, m, s, pi - 1, d, dcStrategy, alphabetCounts) ;
    MultikeyQSort(T, n, sa, m, pi, pk, d + 1, dcStrategy, alphabetCounts) ;
    MultikeyQSort(T, n, sa, m, pj, e, d, dcStrategy, alphabetCounts) ;
  }

  // Sort the suffixes in the difference cover using Manber-Myers algorithm
  // @return: the suffix array of the difference cover
  size_t *SortSuffixInDCWithMM(const FixedSizeElemArray &T, size_t n)
  {
    size_t i, k ;
    size_t *sa ;
    size_t *rank ; // rank of a suffix consider the prefix of size k, allowing ties
    size_t *nextBuffer ; // buffer for next iteration (double expanded) information
    size_t *count ; // cumulative rank count 
    size_t maxRank ; // distinct ranks
    size_t dci ; // compacted difference cover index
    int v = _dc.GetV() ;
    
    sa = (size_t *)malloc(sizeof(size_t) * _dcSize) ;
    rank = (size_t *)malloc(sizeof(size_t) * _dcSize) ;
    nextBuffer = (size_t *)malloc(sizeof(size_t) * _dcSize) ;
    count = (size_t *)malloc(sizeof(size_t) * _dcSize) ;
    
    // Sort by their first v characters
    // It has to be at least v characters, and v is a multiple of 2. 
    // This is because the following Manber-Myers algorithm (or Larsson-Sadakane) 
    //  needs to start at k=v, so the +k is also in the difference cover.
    //Utils::PrintLog("SA sort start") ;
    _dc.GetDiffCoverList(n, sa) ;
    size_t *alphabetCounts = (size_t *)malloc(sizeof(size_t) * (_alphabetSize + 1)) ;
    MultikeyQSort(T, n, sa, _dcSize, 0, _dcSize - 1, 0, /*dcStrategy=*/2, alphabetCounts) ;
    free(alphabetCounts) ; 
    //Utils::PrintLog("SA sort MultikeyQSort finishes") ;
    
    maxRank = 0 ;
    count[0] = 0 ;
    for (i = 0 ; i < _dcSize ; ++i)
    {
      if (i > 0 && T.SubrangeCompare( sa[i - 1], sa[i - 1] + v - 1, T, sa[i], sa[i] + v - 1))
      {
        ++maxRank ;
        count[maxRank] = 0 ;
      }
      dci = _dc.CompactIndex( sa[i] ) ; 
      rank[dci] = maxRank ;
      ++count[maxRank] ;
    }
    for (i = 1 ; i <= maxRank ; ++i)
      count[i] += count[i - 1] ; 

    // Sorting difference cover using Manber-Myers algorithm 
    size_t *tmpSwap ;
    for (k = v ; k < n /*&& maxRank < _dcSize - 1*/ ; k <<= 1)
    {
      // Get the new SA, nextBuffer serves as the next SA 
      size_t ri ; // reverse i
      for (ri = 0 ; ri < _dcSize; ++ri)
      {
        i = _dcSize - 1 - ri ;

        //dci = _dc.CompactIndex( sa[i] ) ;
        if (sa[i] >= n - k)
          nextBuffer[i] = sa[i] ;
        if (sa[i] < k)
          continue ;
        size_t kbeforeDci = _dc.CompactIndex(sa[i] - k) ;
        nextBuffer[ count[ rank[kbeforeDci] ] - 1 ] = sa[i] - k ;
        --count[ rank[kbeforeDci] ] ;
      }
      //memcpy(sa, nextBuffer, sizeof(sa[0]) * _dcSize) ;
      tmpSwap = sa ;
      sa = nextBuffer ;
      nextBuffer = tmpSwap ;

      // Update the rank, nextBuffer serves as the next rank.
      maxRank = 0 ;
      count[0] = 0 ;
      size_t prevDci = 0 ;
      for (i = 0 ; i < _dcSize ; ++i)
      {
        dci = _dc.CompactIndex(sa[i]) ;
        if (i > 0 && 
            (rank[dci] != rank[prevDci] || 
             sa[i - 1] + k >= n || sa[i] + k >= n ||
             rank[ _dc.CompactIndex(sa[i - 1] + k)] != rank[_dc.CompactIndex(sa[i] + k)]))
        {
          ++maxRank ;
          count[maxRank] = 0 ;
        }
        nextBuffer[dci] = maxRank ; 
        ++count[maxRank] ;
        prevDci = dci ;
      }
      tmpSwap = rank ;
      rank = nextBuffer ;
      nextBuffer = tmpSwap ; 

      if (maxRank >= _dcSize - 1)
        break ;

      for (i = 1 ; i <= maxRank ; ++i)
        count[i] += count[i - 1] ;
      //Utils::PrintLog("SA sort MM iteration done: %llu %llu", k, _dcSize) ;
    }

    free(nextBuffer) ;
    free(count) ;
    
    //Utils::PrintLog("SA sort MM sort finishes") ;
    _dcISA = rank ;
    return sa ;
  }

  // Sort the suffixes in the difference cover using Larsson-Sadakane algorithm
  // @return: the suffix array of the difference cover
  size_t *SortSuffixInDCWithLS(const FixedSizeElemArray &T, size_t n)
  {
    size_t i, j, k ;
    size_t *sa ;
    size_t *rank ; // rank of a suffix consider the prefix of size k, allowing ties. Different from MM, LS stores the highest possible rank for a group. This allows efficient rank updates that can keep the other portion's rank unchanged in doubling
    size_t *nextBuffer ; // buffer for next iteration (double expanded) information
    size_t *L ; // length for single-ton runs (negative), and length for a h-group (k-group here)
    size_t maxRank ; // distinct ranks
    size_t sizeL ; // length of L, number of runs
    int v = _dc.GetV() ;

    sa = (size_t *)malloc(sizeof(size_t) * _dcSize) ;
    rank = (size_t *)malloc(sizeof(size_t) * _dcSize) ;
    nextBuffer = (size_t *)malloc(sizeof(size_t) * _dcSize) ;
    L = (size_t *)malloc(sizeof(size_t) * _dcSize) ;

    // Sort by their first v characters
    Utils::PrintLog("Start to sort |dc|-prefix") ;
    _dc.GetDiffCoverList(n, sa) ;
    size_t *alphabetCounts = (size_t *)malloc(sizeof(size_t) * (_alphabetSize + 1)) ;
    MultikeyQSort(T, n, sa, _dcSize, 0, _dcSize - 1, 0, /*dcStrategy=*/2, alphabetCounts) ;
    free(alphabetCounts) ; 
    Utils::PrintLog("Finish sorting |dc|-prefix") ;

    // Initialization
    size_t count = 0 ;
    sizeL = 0 ;
    for (i = 0 ; i <= _dcSize ; ++i)
    {
      if (i == _dcSize || (i > 0 && T.SubrangeCompare( sa[i - 1], sa[i - 1] + v - 1, T, sa[i], sa[i] + v - 1)))
      {
        if (count == 1) // previous element is a singleton
        {
          if (sizeL == 0 || (int64_t)L[sizeL - 1] > 0)
          {
            L[sizeL] = (size_t)(-1ll) ;
            ++sizeL ;
          }
          else
            L[sizeL - 1] = (size_t)((int64_t)L[sizeL - 1] - 1ll);
        }
        else
        {
          L[ sizeL ] = count ;
          ++sizeL ;
        }
        count = 0 ;

        if (i == _dcSize)
          break ;
      }
      ++count ;
    }
    
    maxRank = 0 ;
    size_t offset = 0 ;
    for (i = 0 ; i < sizeL ; ++i)
    {
      int64_t liValue = (int64_t)L[i] ;
      if (liValue < 0)
      {
        for (j = offset ; j < offset + (size_t)(-liValue) ; ++j)
        {
          size_t dcj = _dc.CompactIndex( sa[j] ) ; 
          rank[dcj] = maxRank ;
          ++maxRank ;
        }
        offset = j ;
      }
      else
      {
        maxRank += L[i] ;
        for (j = offset ; j < offset + L[i] ; ++j)
        {
          size_t dcj = _dc.CompactIndex( sa[j] ) ; 
          rank[dcj] = maxRank - 1 ;
        }
        offset = j ; 
      }
    }

    /*if (maxRank == _dcSize)
      {
      free(nextBuffer) ;
      free(L) ;
      _dcISA = rank ;
      return sa ;
      }*/

    // Sorting difference cover using Larsson-Sadakane algorithm 
    size_t *tmpSwap ;
    for (k = v ; k < n /*&& maxRank < _dcSize - 1*/ ; k <<= 1)
    {
      offset = 0 ;
      for (i = 0 ; i < sizeL ; ++i)
      {
        int64_t liValue = (int64_t)L[i] ;
        if (liValue < 0)
        {
          offset += (size_t)(-liValue) ;
        }
        else
        {
          MultikeyQSortForLSandDC(sa, rank, n, offset, offset + L[i] - 1,
              k, nextBuffer) ;

          offset += L[i] ;
        }
      }

      // Copy the updated rank back from the buffer
      offset = 0 ;
      for (i = 0 ; i < sizeL ; ++i)
      {
        int64_t liValue = (int64_t)L[i] ;
        if (liValue < 0)
        {
          offset += (size_t)(-liValue) ;
        }
        else
        {
          for (j = offset ; j < offset + L[i] ; ++j)
          {
            size_t dcj = _dc.CompactIndex(sa[j]) ;
            rank[dcj] = nextBuffer[dcj] ;
          }
          offset += L[i] ;
        }
      }

      // Update L to nextbuffer, and then swap the points
      size_t newSizeL = 0 ;
      offset = 0 ;
      for (i = 0 ; i < sizeL ; ++i) 
      {
        int64_t liValue = (int64_t)L[i] ;
        if (liValue < 0)
        {
          offset += (size_t)(-liValue) ;
          if (newSizeL == 0 || (int64_t)nextBuffer[ newSizeL - 1] > 0)
          {
            nextBuffer[ newSizeL ] = (size_t)(liValue) ;
            ++newSizeL ;
          }
          else // This happens when the last run in the previous 2k-group is singleton run
          {
            nextBuffer[ newSizeL - 1] = (size_t)((int64_t)nextBuffer[newSizeL-1] + liValue) ;
          }
        }
        else
        {
          size_t count = 1 ;
          for (j = offset + 1 ; j <= offset + liValue ; ++j)
          {
            if (j != offset + liValue 
                && rank[ _dc.CompactIndex(sa[j]) ] == rank[ _dc.CompactIndex(sa[j - 1]) ])
            {
              ++count ;
            }
            else // entering a new 2k-group, or end of the current k-group 
            {
              if (count > 1)
              {
                nextBuffer[ newSizeL ] = count ;
                ++newSizeL ;
                count = 1 ;
              }
              else // singleton
              {
                if (newSizeL == 0 || (int64_t)nextBuffer[ newSizeL - 1] > 0)
                {
                  nextBuffer[ newSizeL ] = (size_t)(-1ll) ;
                  ++newSizeL ;
                }
                else // This happens when the last run in the previous 2k-group is singleton run
                {
                  nextBuffer[ newSizeL - 1] = (size_t)((int64_t)nextBuffer[newSizeL-1] - 1ll) ;
                }
                //count = 1, not need to reset it again
              }
            }
          }

          offset += liValue ;
        }
      }

      if (newSizeL == 1 && (int64_t)nextBuffer[0] == -(int64_t)_dcSize)
        break ;

      tmpSwap = nextBuffer ;
      nextBuffer = L ;
      L = tmpSwap ;
      sizeL = newSizeL ;
    } // for-loop of k

    free(nextBuffer) ;
    free(L) ;

    _dcISA = rank ;
    return sa ;
  }
public:
  SuffixArrayGenerator() 
  {
    _b = 1<<24 ;  // 2^24, 16MB block size by default 
    _n = _space = 0 ;
    _cuts = NULL ;
    _dcISA = NULL ;
  }

  ~SuffixArrayGenerator() 
  {
    Free() ;
  }

  void Free() 
  {
    _space = 0 ;
  
    if (_cuts != NULL)
    {
      free(_cuts) ;
      if (_cutLCP != NULL)
      {
        size_t i ;
        for (i = 0 ; i < _cutCnt ; ++i)
          free(_cutLCP[i]) ;
        free(_cutLCP) ;
      }
    }
    if (_dcISA != NULL)
      free(_dcISA) ;
  }

  size_t GetSpace()
  {
    return _space + sizeof(*this) ;
  }

  // Initialize the generator to obtain the _cuts
  // _dcv: difference cover period
  // @return: the number of _cuts
  size_t Init(const FixedSizeElemArray &T, size_t n, size_t b, int dcv, int alphabetSize)
  {
    this->_n = n ;
    if (b > 0)
      this->_b = b ;
    this->_alphabetSize = alphabetSize ;
    _dc.Init(dcv) ;
    _dcSize = _dc.GetSize(n) ;
		Utils::PrintLog("Start to sort %lu difference cover points.", _dcSize) ;
    size_t *dcSA = SortSuffixInDCWithLS(T, n) ; 
    /*for (size_t i = 0 ; i < _dcSize ; ++i)
      printf("dcSA[%lu]=%lu\n", i, dcSA[i]) ;
    for (size_t i = 0 ; i < _dcSize ; ++i)
      printf("dcISA[%lu]=%lu\n", i, _dcISA[i]) ;*/
    GenerateCuts(dcSA) ; 
		Utils::PrintLog("Start to calculate LCP for %lu cuts.", _cutCnt) ;
    ComputeCutLCP(T, n, dcv) ;
    free(dcSA) ;
		Utils::PrintLog("Finish initializing suffix array generator.") ;
    return _cutCnt ;
  }

  size_t GetChunkCount()
  {
    return _cutCnt ;
  }

  // Estimate how many chunks would be given b and dcv
  // b: user-specific rough block size
  static size_t EstimateChunkCount(size_t n, size_t b, int dcv)
  {
    size_t blockCnt = DIV_CEIL(n, b) ;
    size_t eDcSize = DifferenceCover::EstimateCoverSize(dcv) * DIV_CEIL(n, dcv) ;
    size_t stride = DIV_CEIL(eDcSize, blockCnt) ;
    return DIV_CEIL(eDcSize, stride) ;
  }

  // Generate the from-th chunk to to-th chunk for T[s..e], both are inclusive
  // Each chunk is left close, right open for the cut.
  // The procedure utilized _cutLCP to expediate the search.
  void  GetChunksPositions(const FixedSizeElemArray &T, size_t n, size_t from, size_t to, size_t s, size_t e, std::vector< std::vector<size_t> > &pos)
  {
    size_t i, j ;
    if (to >= _cutCnt)
      to = _cutCnt - 1 ;
    if (e >= n)
      e = n - 1 ;
    std::vector< std::vector<size_t> >().swap(pos) ;
    for (j = from ; j <= to ; ++j)
      pos.push_back( std::vector<size_t>() ) ;
    
    size_t *rightmosti = (size_t *)calloc(to - from + 2, sizeof(size_t)) ;
    size_t *rightmostj = (size_t *)calloc(to - from + 2, sizeof(size_t)) ;
    for (i = s ; i <= e ; ++i)
    {
      if ((from == 0 || CompareCutUsingCutLCP(T, n, i, from, rightmosti[from - from], rightmostj[from - from]) >= 0)
          && (to == _cutCnt - 1 
              || CompareCutUsingCutLCP(T, n, i, to + 1, rightmosti[to + 1 - from], rightmostj[to + 1 - from]) < 0))
      {
        for (j = from + 1 ; j <= to ; ++j)
        {
          if (CompareCutUsingCutLCP(T, n, i, j, rightmosti[j-from], rightmostj[j-from]) < 0)
            break ;
        }
        pos[(j-1) - from].push_back(i) ;
      }
    }
    free(rightmosti) ;
    free(rightmostj) ;
  }

  void SortSuffixByPos(const FixedSizeElemArray &T, size_t n, size_t *pos, size_t m, size_t *sa)
  {
    if (m == 0)
      return ;
    size_t i ;
    if (sa != pos)
      for (i = 0 ; i < m ; ++i)
        sa[i] = pos[i] ;
    size_t *alphabetCounts = (size_t *)malloc(sizeof(size_t) * (_alphabetSize + 1)) ;
    MultikeyQSort(T, n, sa, m, 0, m - 1, 0, /*dcStrategy=*/1, alphabetCounts) ;
    free(alphabetCounts) ;
  }

  // TODO: Functions relating to use disk to hold chunks
  // Output each chunk to prefix_{xxx}.chunk file
  void OutputChunksToFiles(char *prefix)
  {
    int i ;
    char filename[1024] ;
    FILE **fps ;
    fps = (FILE **)malloc(sizeof(FILE*) * _cutCnt) ;
    for (i = 0 ; i < (int)_cutCnt ; ++i)
    {
      sprintf(filename, "%s_%d.chunk", prefix, i) ;
      fps[i] = fopen(filename, "w") ;
    }

    for (i = 0 ; i < (int)_cutCnt ; ++i)
      fclose(fps[i]) ;
    free(fps) ;
  }

  // TODO: Read in the i-th chunk file  
  void ReadChunkFile(char *prefix, int i, std::vector<size_t> &pos)
  {
    char filename[1024] ;
    sprintf(filename, "%s_%d.chunk", prefix, i) ;
    FILE *fp = fopen(filename, "r") ;
    
    fclose(fp) ;
  }

  // TODO: Remove all the temporary chunk files
  void CleanChunkFiles(char *prefix)
  {
    char filename[1024] ;
    for (int i = 0 ; i < (int)_cutCnt ; ++i)
    {
      sprintf(filename, "%s_%d.chunk", prefix, i) ;
    }
  }

  // Use the Theorem 2 from "Fast Lightweight Suffix Array Construction and Checking"
  // The simpler implementation requiring creating ISA, so it is memory intensive
  bool ValidateSA(const FixedSizeElemArray &T, size_t n, size_t *sa)
  {
    size_t i ;
    size_t *isa = (size_t *)malloc(sizeof(size_t) * n) ;
    for (i = 0 ; i < n ; ++i)
    {
      if (sa[i] >= n)
        return false ;
      if (i > 0)
      {
        if (sa[i - 1] == sa[i])
          return false ;
        if (T.Read(sa[i - 1]) > T.Read(sa[i]))
          return false ;
      }
    }
    for (i = 0 ; i < n ; ++i)
      isa[sa[i]] = i ;   
    
    for (i = 1 ; i < n ; ++i)
    {
      if (T.Read(sa[i - 1]) == T.Read(sa[i]))
      {
        if (sa[i-1] + 1 < n && sa[i] + 1 < n)
        {
          if (isa[ sa[i - 1] + 1] > isa[ sa[i] + 1])
            return false ;
        }
        else if (sa[i] + 1 == n) // only the previous can be followed by the end of the string 
          return false ;
      }
    }

    free(isa) ;
    return true ;
  }
} ;
}

#endif

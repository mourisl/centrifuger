#ifndef _MOURISL_COMPACTDS_SEQUENCE_PERMUTATION
#define _MOURISL_COMPACTDS_SEQUENCE_PERMUTATION

// This is for large, uncompact alphabet set
// UNFINISHED WORK

#include "Utils.hpp"
#include "Alphabet.hpp"
#include "FixedSizeElemArray.hpp"
#include "Sequence.hpp"
#include "CompactMapper.hpp"
#include "Permutation.hpp"
#include "Bitvector_Plain.hpp"

namespace compactds {
class Sequence_Permutation: public Sequence
{
private:
  bool _alphabetSparse ;
  Bitvector_Plain _A ; // concatenated alphabet chunk count
  //size_t *_alphabetCount ; // the number of times each alphabet show up in the data
  Bitvector_Plain _alphabetStart ; // record the start position for the i-th alphabet in _A
  Bitvector_Plain _D ; // Concatenated marker of the start position for each number/alphabet in the concatendated permutation list in each chunk. Each chunk use |chunk_size|+|\Sigma| bits. The difference to inverted index is that here there may be unused alphabet
  Permutation *_pi ; // permutation in each block/chunk 

  size_t _m ; // |\Sigma|, which is also the block size.
  CompactMapper _map ; // map alphabet set to consecutive 0..m, where _m is |\Sigma|
public:
  Sequence_Permutation()
  {
    _alphabetSparse = false ;
  }

  ~Sequence_Permutation()
  {
    Free() ;
  }

  Free()
  {
    if (_n > 0)
    {
      free(alphabetCount) ;
      delete[] _D ;
      delete[] _pi ;
      _n = 0 ;
    }
  }

  void SetAlphabetSparsity(bool sparse)
  {
    _alphabetSparse = sparse ;
  }
  
  void Init(const FixedSizeElemArray &S, size_t sequenceLength, const ALPHABET *alphabetMap = NULL) 
  {
    size_t i, j ;
    _n = sequenceLength ;
    
    _map.Init(S, _n, _alphabetSparse) ;
    
    _m = _map.GetCompactSize() ;
    size_t blockCnt = DIV_CEIL(_n, _m) ;
    size_t *alphabetPSum = (size_t *)calloc(m + 1, sizeof(*alphabetCount)) ; 
    for (i = 0 ; i < _n ; ++i)
    {
      size_t a = _map.Map(S.Read(i)) ;
      ++alphabetPSum[a + 1] ;
    }
    
    // Fill in _A
    _A.Malloc(n + blockCnt * _m) ;
    localCount = (size_t *)calloc(m + 1, sizeof(*localCount)) ;

    for (i = 0 ; i < blockCnt ; ++i)
    {
      for (j = i * _m ; j < _n && j < i * (m + 1) ; ++j)   
      {
        size_t a = _map.Map(S.Read(i)) ;
        ++localCount[a] ;
      }
      
      for (j = 0 ; j < _m ; ++j)
      {
        size_t k ;
        size_t offset = alphabetPSum[j] + j * blockCnt + i ; //each block is set as 1^k0, so the +i is for the 0s set before, and j*blockCnt is the 0s from other blocks.
        for (k = 0 ; k < localCount[j] ; ++k)
          _A.BitSet(offset + k) ;
      }

      for (j = 0 ; j < _m ; ++j)
        alphabetPSum[j] += localCount[j] ;
    }
    _A.Init() ;
    
    // Now alphabetPSum looks like count[0], count[0,1], count[0,1,2],...
    for (i = _m + 1 ; i > 0 ; --i)  
      alphabetPSum[i] = alphabetPsum[i - 1] ;
    alphabetPSum[0] = 0 ;
    _alphabetStart.Malloc(n + blockCnt * _m) ;
    for (i = 0 ; i < _m ; ++i)
      _alphabetStart.BitSet(alphabetPSum[i] + i * blockCnt) ;
    _alphabetStart.Init() ;
    
    // Fill in _D and _pi
    _D.Malloc(n + blockCnt * _m) ;
    _pi = new Permutation[blockCnt] ; // can we concatenate it?
    size_t *pi = (size_t *)malloc(sizeof(*pi) * _m) ;
    for (i = 0 ; i < blockCnt ; ++i)
    {
      for (j = 0 ; j < _m ; ++j)
        localCount[j] = 0 ;
      for (j = i * _m ; j < _n && j < i * (m + 1) ; ++j)
        ++localCount[ _map.Map(S.Read(j)) ] ;        
      
      size_t psum = 0 ;
      for (j = 0 ; j < _m ; ++j)
      {
        psum += localCount[j] ;
        _D.BitSet(i * (2 * _m) + psum + j) ; // each block use 2*m bits in _D.
      }
      
      // Make the permutation
      for (j = _m ; j > 0 ; --j)
        localCount[j] = localCount[j - 1] ;
      localCount[0] = 0 ;
      for (j = i * _m ; j < _n && j < i * (m + 1) ; ++j)
      {
        size_t l = _map.Map(S.Read(j)) ;
        pi[ localCount[l] ] = j - i * _m ;
        ++localCount[l] ;
      }
      _pi[i].Init(pi, localCount[m]) ;
    }
    _D.Init() ;

    free(pi) ;
    free(localCount) ;
  }

  void Free()
  {
    if (_n > 0)
    {
      _n = 0 ;
      size_t _m = _map.GetCompactSize() ;
      size_t blockCnt = DIV_CEIL(_n, _m) ;
      for (i = 0 ; i < blockCnt ; ++i)
        _pi[i].Free() ;
      _A.Free() ;
      _alphabetStart.Free() ;
      _D.Free() ;
      _mapper.Free() ;
    }
  }
  
  size_t GetSpace() 
  {
  }
  
  ALPHABET Access(size_t i) const 
  {
    return AccessLong(i) ;
  }

  size_t Rank(ALPHABET c, size_t i, int inclusive = 1) const
  {
    return RankLong(c, i, inclusive) ;
  }

  size_t Select(ALPHABET c, size_t i) const
  {
    return SelectLong(c, i) ;
  }

  // These .Long function goes beyond the limitation of ALPHABET
  size_t AccessLong(size_t i) const 
  {
    size_t k = i / _m ;
    size_t p = _pi[k].Next(i % _m) ;
    // No need to use rank agian, the excess is succificient enough
    return _map.MapBack( _D.Select(0, p + k * _m) - p - k * 2 * _m) ;
  }
  
  size_t RankLong(size_t c, size_t i, int inclusive = 1) const
  {
		if (!inclusive)
		{
			if (i == 0)
				return 0 ;
			--i ;
		}

    size_t k = i / _m ;
    //size_t j = i % _m ;
    c = _map.Map(c) ;

    size_t m, l, r ; // left, right for the binary search. This is why we need to split into blocks, otherwise the binary search range will be huge.
		if (!(c == 0 && k == 0))
      l = _D.Select(1, c + k * _m - 1) ;
		else 
			l = 0 ;
		r = _D.Select(1, c + k * _m) - 1 ;
    size_t lstart = l ; 

		while (l <= r)
		{
			m = (l + r) / 2 ;
			if (_pi[k].Next(m) <= i)
				l = m + 1 ;
			else
			{
				if (r == 0)
					break ;
				r = m - 1 ;
			}
		}

		// Use _A
		size_t aOffset = _alphabetStart.Select(1, c + 1) ;
		if (c == 0 && k == 0)
			return r ; 
		size_t cntBefore = _A.Select(0, _m * c + k) - aOffset - k ; // number of "c" show up before the current block
		return cntBefore + r - lstart ;
  }

  size_t SelectLong(size_t c, size_t i) const
  {
		if (i < 1 || i > _n)
			return POSITIVE_INF ;
		
    c = _map.Map(c) ;
		
		size_t aOffset = _alphabetStart.Select(1, c + 1) ; 
		size_t s = _A.Select(1, i + _A.Rank(1, aOffset) - 1) ;
		if (s >= _alphabetStart.Select(1, c + 2))
			return POSITIVE_INF ;
		size_t j = s - _A.Pred0(_A, s) ;
		size_t l = _D.Select(1, c) - c ;  
		
		return (s - j) * _m + _pi[s - j].Next(l + j) ;
  }

  void Save(FILE *fp)
  {
    Sequence::Save(fp) ;
  }

  void Load(FILE *fp)
  {
    Sequence::Load(fp) ;
  }
  
  void PrintStats() 
  {
  }
} ;
}

#endif

#ifndef _MOURISL_COMPACTDS_BITVECTOR_PLAIN
#define _MOURISL_COMPACTDS_BITVECTOR_PLAIN

#include "Utils.hpp"
#include "Bitvector.hpp"

#include "DS_Rank.hpp"
#include "DS_Select.hpp"

// The bitvector with  
namespace compactds {
class Bitvector_Plain : public Bitvector
{
private:
  size_t _n ; // the total raw length of the bits  
  
  // Variables for the bit vector
  WORD *_B ; // bitvector packed in WORD array

  // Variables for _ranking query
  DS_Rank9 _rank ;
  int _rb ; 

  // Variables for _selection
  DS_Select _select ;
  int _sb ;

  int _selectSpeed ;
  int _selectTypeSupport ;

public:
  Bitvector_Plain() 
  {
    _n = _rb = _sb = 0 ;
    _B = NULL ;
    _selectSpeed = BITVECTOR_DEFAULT_SELECT_SPEED ;
    _selectTypeSupport = 3 ;
  }
  ~Bitvector_Plain() {Free();}
  
  void SetRankBlockLength(int rBlockSize)
  {
    _rb = rBlockSize ;
  }

  void SetSelectBlockLength(int sBlockSize)
  {
    _sb = sBlockSize ;
  }

  void SetSelectSpeed(int selectSpeed)
  {
    this->_selectSpeed = selectSpeed ;
  }

  void SetSelectTypeSupport(int selectTypeSupport)
  {
    this->_selectTypeSupport = selectTypeSupport ;
  }

  
  void Malloc(const size_t &n)
  {
    this->_n = n ;
    _B = Utils::MallocByBits(n) ;
    
    _space = Utils::BitsToWordBytes(n) ;
  }
  
  void Free()
  {
    if (_B != NULL)
    {
      free(_B) ;
      _B = NULL ;
    }
    _rank.Free() ;
    _select.Free() ;
    _n = 0 ;
  }

  // Use with caution that the _rank and 
  void BitSet(size_t i)
  {
    Utils::BitSet(_B, i) ;
  }

  void BitClear(size_t i)
  {
    Utils::BitClear(_B, i) ;
  }

  // W is the plain bit vector
  void Init(const WORD *W, const size_t n)
  {
    _n = n ;
    Malloc(n) ;
    memcpy(_B, W, Utils::BitsToWordBytes(n)) ;

    Init() ;
  }
 
  // This is for when _B is already allocated
  void Init()
  {
    _space = Utils::BitsToWordBytes(_n) ;
    _rank.Free() ;
    _select.Free() ;
    //_rank.Init(_rb, _B, _n) ;
    _rank.Init(_B, _n) ;
    _space += _rank.GetSpace() - sizeof(_rank) ;
    _select.Init(_sb, _B, _n, _selectSpeed, _selectTypeSupport) ;
    _space += _select.GetSpace() - sizeof(_select) ;
  }

  void InitFromOnes(const uint64_t *S, const size_t onecnt, const size_t n)
  {
    Malloc(n) ;

    size_t i ;
    for (i = 0 ; i < onecnt ; ++i)
      BitSet(S[i]) ;

    Init() ;
  }

  // Return the ith bits (0-based)
  int Access(size_t i) const
  {
    return Utils::BitRead(_B, i) ;
  }

  // Return the number of 1s before i
  size_t Rank1(size_t i, int inclusive = 1) const
  {
    return _rank.Query(i, _B, _n, inclusive) ;
  }

  // Return the index of th i-th (this i is 1-based, so rank and select are inversible) 1
  size_t Select(size_t i) const
  {
    return _select.Query(i, _rank, _B, _n) ;
  }

  size_t Select(int type, size_t i) const
  {
    if (type == 1)
      return _select.Query(i, _rank, _B, _n) ;
    else
      return _select.Query0(i, _rank, _B, _n) ;
  }

  // Pred/successor on bit 0
  size_t Pred0(size_t i) const
  {
    return Select(0, Rank(0, i)) ;
  }

  size_t Succ0(size_t i) const
  {
    return Select(0, Rank(0, i, 0) + 1) ; 
  }

  size_t GetSpace() 
  {
    return _space + sizeof(*this) ;
  }

  const WORD *GetData() const
  {
    return _B ;
  }

  void Print(FILE *fp)
  {
    size_t i ;
    for (i = 0 ; i < _n ; ++i)
      fprintf(fp, "%d", Access(i)) ;
    fprintf(fp, "\n") ;
  }

  void Save(FILE *fp)
  {
    Bitvector::Save(fp) ;
    SAVE_VAR(fp, _n) ;
    SAVE_VAR(fp, _rb) ;
    SAVE_VAR(fp, _sb) ;
    SAVE_VAR(fp, _selectSpeed) ;
    SAVE_VAR(fp, _selectTypeSupport) ;
    if (_n > 0)
    {
      fwrite(_B, sizeof(*_B), Utils::BitsToWords(_n), fp) ;
      _rank.Save(fp) ;
      _select.Save(fp) ;
    }
  } 

  void Load(FILE *fp)
  {
    Free() ;
    Bitvector::Load(fp) ;
    LOAD_VAR(fp, _n) ;
    LOAD_VAR(fp, _rb) ;
    LOAD_VAR(fp, _sb) ;
    LOAD_VAR(fp, _selectSpeed) ;
    LOAD_VAR(fp, _selectTypeSupport) ;
    
    if (_n > 0)
    {
      _B = Utils::MallocByBits(_n) ;
      fread(_B, sizeof(*_B), Utils::BitsToWords(_n), fp) ;
      _rank.Load(fp) ;
      _select.Load(fp) ;
    }
    else
    {
      _B = NULL ;
      //_rank.Free() ;
      //_select.Free() ;
    }
  }
} ;
}

#endif

#ifndef _MOURISL_COMPACTDS_DUSTMASKER
#define _MOURISL_COMPACTDS_DUSTMASKER

// SDust algorithm
// Mask the low-complexity regions in the input sequence. Return the low-complexity intervals (0-based, inclusive on both end) 
// The algorithm is based on the DUST algorithm described in Morgulis et al. "A Fast and Symmetric DUST Implementation to Mask Low-Complexity DNA Sequences". J Comput Biol. 2006 Apr;13(5):1028-40. doi: 10.1089/cmb.2006.13.1028. PMID: 16646928.

#include <vector>

#include <algorithm>

struct _dustmasker_perfect_interval
{
  size_t start ;
  size_t end ;
  int score ; // kind of rv, not the final score. score / (end - start - 2) * 10 is the final score. We store the score numerator to avoid the precision loss from the division.
} ;


class Dustmasker
{
private: 
  int _w ; // window size
  int _T ; // threshold for masking
  int _l ; // linker
  
  int *_alphabetMap ;
  int _alphabetSize ;
  int _alphabetBit ;

  // queue data structure just for the purpose of this function.
  // We assume the capacity is fixed
  class Dustmasker_Queue
  {
    private:
      int _head ;
      int _tail ;
      int _mask ;
      int *_s ;
    public:
      Dustmasker_Queue()
      {
        _head = _tail = _mask = 0 ;
        _s = NULL ;
      }

      Dustmasker_Queue(int sz)
      {
        _head = _tail = 0 ;
        int capacityBits = 0 ;
        while ((1 << capacityBits) <= sz)
          ++capacityBits ;
        _s = (int *)malloc(sizeof(int) * (1 << capacityBits)) ;
        _mask = (1 << capacityBits) - 1 ;
      }

      ~Dustmasker_Queue()
      {
        if (_s != NULL)
          free(_s) ;
      }

      int Size()
      {
        return (_tail - _head) & _mask ;
      }

      void PushBack(int t)
      {
        _s[_tail] = t ;
        _tail = (_tail + 1) & _mask ;
      }

      int PopFront()
      {
        int t = _s[_head] ;
        _head = (_head + 1) & _mask ;
        return t ;
      }

      int Front()
      {
        return _s[_head] ;
      }

      int operator[](int i)
      {
        return _s[(_head + i) & _mask] ;
      }
  } ;

  // r is the score without normalize the length. It increase count[tripletCode] when adding a triplet.
  void AddTripletInfo(int tripletCode, int *count, int &r) 
  {
    r += count[tripletCode] ;
    ++count[tripletCode] ;
  }

  void RemoveTripletInfo(int tripletCode, int *count, int &r) 
  {
    --count[tripletCode] ;
    r -= count[tripletCode] ;
  }

  // rv, cv are for the longest suffix of the current window that satisfy max{c(v)}<=2T. lv tracks the length of this suffix
  void ShiftWindow(int t, Dustmasker_Queue &window, int &lv, int &rw, int &rv, int *cw, int *cv) 
  {
    if ((int)window.Size() >= _w - 2) //window stores triplet, so the actual bases cover by the window is window.size() + 2 
    {
      int oldTripletCode = window.Front();
      RemoveTripletInfo(oldTripletCode, cw, rw) ;
      window.PopFront();
      if (lv > (int)window.Size())
      {
        RemoveTripletInfo(oldTripletCode, cv, rv) ;
        --lv ;
      }
    }

    window.PushBack(t);
    ++lv ;
    AddTripletInfo(t, cw, rw) ;
    AddTripletInfo(t, cv, rv) ;

    if (cv[t] * 10 > 2 * _T) // Now the suffix does not satisfy the condition 1, we need to shrink the suffix
    {
      while (1)
      {
        int s = window[ window.Size() - lv] ;
        RemoveTripletInfo(s, cv, rv) ;
        --lv ;
        if (s == t) // The cv[t] will drop to 2*_T at this moment, so the suffix will satisfy the condition again. We can stop here.
          break ;
      }
    }
  }

  // We are processing the window from windowStart, so the perfect intervals that start before windowStart need to be saved to he result.
  void SaveMaskedRegions(std::vector<struct _dustmasker_perfect_interval> &result, std::vector<struct _dustmasker_perfect_interval> &P, size_t windowStart)
  {
    // P is mainteined in sorted order, where it is sorted first by descending order of start and then by ascending order of end. 
    // Based on how we maintain P, the last P is the one from windowStart - 1, and is the longest one including all the other perfect intervals starting before windowStart.
    if (P.size() > 0 && P.back().start < windowStart) 
    {
      struct _dustmasker_perfect_interval lastP = P.back() ;
      size_t l = result.size() ;
      if (l > 0)
      {
        if (lastP.start <= result[l - 1].end + 1) // merge the intervals. The interval is [start,end], closed.
        {
          result[l - 1].end = std::max(result[l - 1].end, lastP.end) ;
        }
        else
        {
          result.push_back(lastP) ;
        }
      }
      else
      {
        result.push_back(lastP) ;
      }

      while (P.size() > 0 && P.back().start < windowStart) 
      {
        P.pop_back() ;
      }
    }
  }

  // window: triplets of the window
  // windowStart: the start position of the window in the sequence.
  // lv is the length of the longest suffix of the current window that satisfy max{c(v)}<=2T. rv is the score of this suffix. cv is the count of triplets in this suffix.
  void FindPerfect(std::vector<struct _dustmasker_perfect_interval> &P, Dustmasker_Queue &window, size_t windowStart, int lv, int rv, int *cv)
  {
    int i ;
    int maxScore = 0 ;
    int maxScoreTripletCount = 1 ;
    //int *tmpc = (int *)calloc(64, sizeof(int)) ;
    //int tmpc[512] ;
    //memcpy(tmpc, cv, sizeof(int) * 512) ;
    //int oldrv = rv ;

    for (i = window.Size() - lv - 1 ; i >= 0 ; --i)
    {
      int t = window[i] ;
      AddTripletInfo(t, cv, rv) ;
      std::vector<struct _dustmasker_perfect_interval>::iterator it = P.begin() ;
      //int newScore =  rv * 10 / (window.size() - i - 1) ; // the score of the suffix starting form position i 
      /*{
        int n = window.size() - i ;
        if (n * (n - 1) / 2 != rv)
        {
          printf("Error: the score is not correct for the suffix starting from position %d. The score is %d, but it should be %d\n", i, rv, (n * (n - 1) / 2)) ;
        }
        else
          printf("good\n") ;
      }*/

      //printf("%d %d: %d %d %d %d\n", newScore, _T, windowStart, i, rv, window.size()) ;
      if (rv * 10 > _T * (window.Size() - i - 1))
      {
        // When the prefect interval is inside of the current suffix
        while (it != P.end() && it->start >= i + windowStart)
        {
          if ((size_t)it->score * maxScoreTripletCount > maxScore * (it->end - it->start - 2)) 
          {
            maxScore = it->score ;
            maxScoreTripletCount = it->end - it->start - 2 ;
          }
          ++it ;
        }

        if (rv * maxScoreTripletCount >= maxScore * (window.Size() - i - 1)) // Perfect interval requires that no subinterval has higher score.
        {
          maxScore = rv ;
          maxScoreTripletCount = window.Size() - i - 1 ;

          struct _dustmasker_perfect_interval newPerfectInterval ;
          newPerfectInterval.start = i + windowStart ;
          newPerfectInterval.end = windowStart + window.Size() + 1 ; //+1 is for the triplet
          newPerfectInterval.score = rv ;
          //printf("%lu %lu %d\n", newPerfectInterval.start, newPerfectInterval.end, newPerfectInterval.score) ;
          P.insert(it, newPerfectInterval) ; // insert the new interval before the iterator position.
        }
      }
    }

    // Resume the cv to the original state
    for (i = window.Size() - lv - 1 ; i >= 0 ; --i)
    {
      int t = window[i] ;
      RemoveTripletInfo(t, cv, rv) ;
    }
    //memcpy(cv, tmpc, sizeof(int) * 512) ;
    //rv = oldrv ;

    /*for (i = 0 ; i < 64 ; ++i)
    {
      if (cv[i] != tmpc[i])
        printf("Error: cv[%d] is changed from %d to %d\n", i, tmpc[i], cv[i]) ;
    }*/
  }

public:
  Dustmasker() 
  {
    _w = 64 ; // default window size
    _T = 20 ; // Based on the paper, the default threshold is 2, but the dustmasker program multipied the _T and S(a) by 10.
    _l = 1 ;
    _alphabetMap = new int[256] ;
    for (int i = 0 ; i < 256 ; ++i)
    {
      _alphabetMap[i] = -1 ;
    }
  }

  ~Dustmasker() 
  {
    delete[] _alphabetMap ;
  }

  void SetWindowSize(int w)
  {
    _w = w;
  }

  int GetWindowSize()
  {
    return _w ;
  }

  void SetThreshold(int T)
  {
    _T = T ;
  }

  void SetLinkerSize(int l)
  {
    _l = l ;
  }

  void Init(const char *alphabetMap)
  {
    int i, n ;
    
    for (i = 0 ; i < 256 ; ++i)
    {
      _alphabetMap[i] = -1 ;
    }

    for (i = 0 ; alphabetMap[i] != 0 ; ++i)
    {
      _alphabetMap[(int)alphabetMap[i]] = i ;
    }

    n = i ;

    for (i = 0 ; i < 256 ; ++i)
    {
      if (_alphabetMap[i] < 0)
        _alphabetMap[i] = n ; // all special character mapped to the same code
    }
    ++n ;

    _alphabetSize = n ;
    _alphabetBit = 0 ;
    while ((1 << _alphabetBit) < n)
      ++_alphabetBit ;
  }

  // The sdust algorithm
  void SDust(const char *S, size_t n, std::vector<struct _dustmasker_perfect_interval> &result)
  {
    if (n < 3)
      return ;

    size_t wstart, wfinish ;

    int triplet = 0 ;
    int tripletMask = (1<<(3 * _alphabetBit)) - 1 ; // mask for the triplet code 
    if (n < 3)
      return ;
    
    int *countV = (int *)calloc(tripletMask + 1, sizeof(int) ); // cv: count for the suffix v that satisfy max{c(v)}<=2T
    int *countW = (int *)calloc(tripletMask + 1, sizeof(int) ); // cw: count for the current window 
    int rv = 0, rw = 0, lv = 0 ;

    Dustmasker_Queue window(_w) ; // store the triplet code in the current window. The actual bases covered by the window is window.size() + 2.
    std::vector<struct _dustmasker_perfect_interval> P ; 
    
    triplet = (_alphabetMap[(int)S[0]] << _alphabetBit) + _alphabetMap[(int)S[1]] ;
    for (wfinish = 2 ; wfinish < n ; ++wfinish)
    {
      size_t wstart = 0 ;
      if (wfinish + 1 > (size_t)_w)
        wstart = wfinish + 1 - _w ;
      SaveMaskedRegions(result, P, wstart) ;
      
      triplet = ((triplet << _alphabetBit) & tripletMask) + _alphabetMap[(int)S[wfinish]] ;
      ShiftWindow(triplet, window, lv, rw, rv, countW, countV) ;
      //printf("%d %d %d. %d %d. %d. %d\n", rw, lv ,_T, wstart, wfinish, triplet, P.size()) ;
      if (rw * 10 > lv * _T) // The current window does not satisfy the condition 2. So it can have perfect interval inside.
        FindPerfect(P, window, wstart, lv, rv, countV) ;
    }
    wstart = 0 ;
    if (wfinish + 1 > (size_t)_w)
      wstart = wfinish + 1 - _w ;
    while (P.size() > 0)
    {
      SaveMaskedRegions(result, P, wstart) ;
      ++wstart ;
    }
    
    free(countV) ;
    free(countW) ;
  }
  
  // The main function to do dustmasking. Handling non-specific characters, like Ns, and conduct merging nearby low-complex intervals based on the _linker function.
  void MaskWithBuffer(const char *S, size_t n, std::vector<struct _dustmasker_perfect_interval> &windowResult, std::vector<struct _dustmasker_perfect_interval> &result)
  {
    size_t i, j ;
    
    result.clear() ;
    if (n < 3)
      return ;

    // Skip the Ns at the beginning
    for (i = 0 ; i < n && _alphabetMap[(int)S[i]] == _alphabetSize - 1 ; ++i)
      ;

    for (; i < n ; )
    {
      size_t nCount = 0 ;
      size_t lastValidPos = i ;
      for (j = i ; j < n ; ++j)
      {
        if (_alphabetMap[(int)S[j]] == _alphabetSize - 1)
        {
          ++nCount ;
        }
        else
        {
          if (nCount > (size_t)_w) // Get out of a very long run of Ns.
            break ;
          lastValidPos = j ;
          nCount = 0 ;
        }
      }

      if (lastValidPos > i)
      {
        windowResult.clear() ;
        SDust(S + i, lastValidPos - i + 1, windowResult) ;
        for (size_t k = 0 ; k < windowResult.size() ; ++k)
        {
          windowResult[k].start += i ;
          windowResult[k].end += i ;
          result.push_back(windowResult[k]) ;
        }
      }

      i = j ;
    }

    // Merge nearby low-complex intervals.
    if (_l > 1 && result.size() > 0)
    {
      size_t size = result.size() ;
      size_t k = 0 ;
      for (i = 1 ; i < size ; ++i)
      {
        if (result[i].start <= result[k].end + _l) // merge the intervals. The interval is [start,end], closed.
        {
          result[k].end = std::max(result[k].end, result[i].end) ;
        }
        else
        {
          ++k ;
          result[k] = result[i] ;
        }
      }
    }
  }

  // Thread-safe wrapper 
  void Mask(const char *S, size_t n, std::vector<struct _dustmasker_perfect_interval> &result)
  {
    std::vector<struct _dustmasker_perfect_interval> windowResult ;
    MaskWithBuffer(S, n, windowResult, result) ;
  }
} ;

#endif

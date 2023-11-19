#ifndef _MOURISL_COMPACTDS_VARIABLESIZEELEM_ARRAY
#define _MOURISL_COMPACTDS_VARIABLESIZEELEM_ARRAY

#include <vector>

#include "Utils.hpp"
#include "FixedSizeElemArray.hpp"

/*
 * The class for the array where each element has variable size
 */
namespace compactds {
class VariableSizeElemArray
{
public:
  VariableSizeElemArray() {}

  ~VariableSizeElemArray() {}

  virtual void Free() = 0;
  
  // Create the variable size element array  
  // b - block size
  // in - input array
  // n - the length of input array
  // 
  virtual void InitFromArray(int b, const unsigned int *in, const size_t &n) = 0 ;

  virtual unsigned int Read(size_t i) = 0 ;
} ;
}

#endif

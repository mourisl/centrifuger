#ifndef _LSONG_DEFS_HEADER
#define _LSONG_DEFS_HEADER

#include <stdint.h>

//#define DEBUG

#define CENTRIFUGER_VERSION "1.0.12-r273"

#define SAMPLE_SHEET_SEPARATOR_READ_ID "__centrifuger_sample_sheet_separator__"

extern char nucToNum[26] ; 
extern char numToNuc[26] ;

//#define MAX(x,y) (((x)>(y))?(x):(y))
//#define MIN(x,y) (((x)<(y))?(x):(y))
#define ABS(x) (((x)>(0))?(x):(-(x)))

#define MAX_SEG_COUNT 127
struct _pair
{
	int a, b ;
} ;

struct _pair64
{
	int64_t a, b ;
} ;

struct _triple
{
	int a, b, c ;
} ;
/*struct _pair_b64
{
	int a ;
	uint64_t b ;
} ;*/
#endif

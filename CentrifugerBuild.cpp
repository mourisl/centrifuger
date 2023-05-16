#include <stdio.h>
#include <time.h>
#include <getopt.h>

#include "argvdefs.h"
#include "ReadFiles.hpp"
#include "compactds/Sequence_Hybrid.hpp"
#include "compactds/FMBuilder.hpp"
#include "Taxonomy.hpp"

char nucToNum[26] = { 0, -1, 1, -1, -1, -1, 2, 
	-1, -1, -1, -1, -1, -1, 0,
	-1, -1, -1, -1, -1, 3,
	-1, -1, -1, -1, -1, -1 } ;

static const char *short_options = "f:o:t:" ;
static struct option long_options[] = {
			{ "bmax", required_argument, 0, ARGV_BMAX},
			{ "dcv", required_argument, 0, ARGV_DCV},
      { "offrate", required_argument, 0, ARGV_OFFRATE},
      { "convertion-table", required_argument, 0, ARGV_CONVERTION_TABLE},
			{ "name-table", required_argument, 0, ARGV_NAME_TABLE},
			{ (char *)0, 0, 0, 0} 
			} ;


int main(int argc, char *argv[])
{
	if ( argc <= 1 )
	{
		//fprintf( stderr, "%s", usage ) ;
		return 0 ;
  }
	int c, option_index ;
	option_index = 0 ;
  FMBuilder fmBuilder ;
  while (1)
  {
		c = getopt_long( argc, argv, short_options, long_options, &option_index ) ;
		
		if (c == -1)
			break ;
  
    if (c == 'f')
    {

    }
    else if (c == 'o')
    {

    }
    else if (c == 't')
    {

    }
  }
  FixedSizeElemArray genomes ;
  
  //fmBuilder.Build(genome, genomeLength,)
  
  return 0 ;
}

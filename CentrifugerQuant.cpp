#include <stdio.h>
#include <time.h>
#include <getopt.h>

#include "argvdefs.h"
#include "Taxonomy.hpp"
#include "Quantifier.hpp"

char usage[] = "./centrifuger-quant [OPTIONS]:\n"
	"Required:\n"
  "\t-x FILE: index prefix\n"
  "\t-c FILE: classification result file\n"
  //"When not giving "-x"\n"
  //"\t--taxonomy-tree FILE: taxonomy tree, i.e., nodes.dmp file\n"
  //"\t--name-table FILE: name table, i.e., names.dmp file\n"
  //"\t--size-table FILE: size table, table of contig (or genome) sizes\n"
  //"\t--conversion-table FILE: a table that converts any id to a taxonomy id\n"
	;

static const char *short_options = "x:c:" ;
static struct option long_options[] = {
  { (char *)0, 0, 0, 0} 
} ;

int main(int argc, char *argv[])
{
  if ( argc <= 1 )
  {
    fprintf( stderr, "%s", usage ) ;
    return 0 ;
  }

  int c, option_index ;
  option_index = 0 ;

  char *idxPrefix = NULL ; 
  char *classificationFile = NULL ;
	
  while (1)
  {
    c = getopt_long( argc, argv, short_options, long_options, &option_index ) ;

    if (c == -1)
      break ;

    if (c == 'x')
    {
      idxPrefix = strdup(optarg) ;
    }
    else if (c == 'c')
    {
      classificationFile = strdup(optarg) ;
    }
    else
    {
      fprintf(stderr, "Unknown parameter found.\n%s", usage) ;
      return EXIT_FAILURE ;
    }
  }
  
	Utils::PrintLog("Centrifuger-quant v" CENTRIFUGER_VERSION " starts." ) ;
  if (idxPrefix == NULL)
  {
    Utils::PrintLog("Need to use -x to specify index prefix.") ;
    return EXIT_FAILURE ;
  }
  
  Quantifier quantifier ;
  quantifier.Init(idxPrefix) ;
  quantifier.LoadReadAssignments(classificationFile, 0) ;

  quantifier.Quantification() ;

  quantifier.Output(stdout, 0) ;

  free(idxPrefix) ;
  free(classificationFile) ;
  
	Utils::PrintLog("Centrifuger-quant finishes." ) ;
  return 0 ;
}

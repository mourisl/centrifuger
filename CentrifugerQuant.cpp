#include <stdio.h>
#include <time.h>
#include <getopt.h>

#include "argvdefs.h"
#include "Taxonomy.hpp"
#include "Quantifier.hpp"

char usage[] = "./centrifuger-quant [OPTIONS]:\n"
	"Required:\n"
  "\t-c FILE: classification result file\n"
  "\t-x FILE: index prefix\n"
  "\tWhen not giving -x\n"
  "\t\t--taxonomy-tree FILE: taxonomy tree, i.e., nodes.dmp file\n"
  "\t\t--name-table FILE: name table, i.e., names.dmp file\n"
  "\t\t--size-table FILE: size table (optional), table of contig (or genome) sizes.\n"
  "Optional:\n"
  "\t--min-score INT: only consider reads with score at least <int> \n"
  "\t--min-length INT: only consider reads with classified length at least <int>\n"
  ""
	;

static const char *short_options = "x:c:" ;
static struct option long_options[] = {
  { "taxonomy-tree", required_argument, 0, ARGV_TAXONOMY_TREE},
  { "name-table", required_argument, 0, ARGV_NAME_TABLE},
  { "size-table", required_argument, 0, ARGV_SIZE_TABLE},
  { "min-score", required_argument, 0, ARGV_QUANT_MINSCORE},
  { "min-length", required_argument, 0, ARGV_QUANT_MINLENGTH},
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
	size_t classificationMinScore = 0 ;
  int classificationMinLength = 0 ;
  
  char *taxonomyFile = NULL ; // taxonomy tree file
  char *nameTable = NULL ;
  char *sizeTable = NULL ;

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
    else if (c == ARGV_QUANT_MINSCORE)
    {
      classificationMinScore = atoi(optarg) ;
    }
    else if (c == ARGV_QUANT_MINLENGTH)
    {
      classificationMinLength = atoi(optarg) ;
    }
    else if (c == ARGV_TAXONOMY_TREE)
    {
      taxonomyFile = strdup(optarg) ;  
    }
    else if (c == ARGV_NAME_TABLE)
    {
      nameTable = strdup(optarg) ;
    }
    else if (c == ARGV_SIZE_TABLE)
    {
      sizeTable = strdup(optarg) ;
    }
    else
    {
      fprintf(stderr, "Unknown parameter found.\n%s", usage) ;
      return EXIT_FAILURE ;
    }
  }
  
	Utils::PrintLog("Centrifuger-quant v" CENTRIFUGER_VERSION " starts." ) ;
  if (idxPrefix == NULL && 
      (taxonomyFile == NULL || nameTable == NULL))
  {
    Utils::PrintLog("Need to use -x to specify index prefix, or --taxonomy-tree/--name-table/(--size-table) to specify the taxonomy structure.") ;
    return EXIT_FAILURE ;
  }
  
  Quantifier quantifier ;
  if (idxPrefix)
    quantifier.Init(idxPrefix) ;
  else
    quantifier.Init(taxonomyFile, nameTable, sizeTable) ;
  quantifier.LoadReadAssignments(classificationFile, classificationMinScore, classificationMinLength, 0) ;
  Utils::PrintLog("Finish loading the read classification result.") ;

  quantifier.Quantification() ;

  quantifier.Output(stdout, 0) ;

  free(classificationFile) ;
 
  if (idxPrefix != NULL)
    free(idxPrefix) ;
  else
  {
    free(taxonomyFile) ;
    free(nameTable) ;
    free(sizeTable) ;
  }

	Utils::PrintLog("Centrifuger-quant finishes." ) ;
  return 0 ;
}

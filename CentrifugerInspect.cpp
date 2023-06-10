#include <stdio.h>
#include <getopt.h>

#include "Taxonomy.hpp"
#include "argvdefs.h"

char usage[] = "./centrifuger-inspect [OPTIONS]:\n"
  "Required:\n"
  "\t-x STRING: index prefix"
  "One of:\n"
  "\t--summary: \n"
  "\t--name: \n"
  "\t--conversion-table: \n"
  "\t--taxonomy-tree: \n"
  "\t--name-table: \n"
  ""
  ;

static const char *short_options = "x:" ;
static struct option long_options[] = {
  {"summary", no_argument, 0, ARGV_INSPECT_SUMMARY},
  {"name", no_argument, 0, ARGV_INSPECT_SEQNAME},
  {"conversion-table", no_argument, 0, ARGV_CONVERSION_TABLE},
  {"taxonomy-tree", no_argument, 0, ARGV_TAXONOMY_TREE},
  {"name-table", no_argument, 0, ARGV_NAME_TABLE},
  { (char *)0, 0, 0, 0} 
} ;

int main(int argc, char *argv[])
{
  char buffer[1024] ;
  char *idxPrefix = NULL ; 
  int c, option_index ;
	option_index = 0 ;

  Taxonomy taxonomy ;
  int inspectItem = -1 ; 
  while (1)
  {
		c = getopt_long( argc, argv, short_options, long_options, &option_index ) ;
		
		if (c == -1)
			break ;
    
    if (c == 'x') // reference genome file
    {
      idxPrefix = strdup(optarg) ;
    }
    else
      inspectItem = c ; 
  }

  if (idxPrefix == NULL)
  {
    fprintf(stderr, "Need -x to specify index.\n%s", usage) ;
    return EXIT_FAILURE ;
  }
  
  sprintf(buffer, "%s.2.cfr", idxPrefix) ;
  FILE *fp = fopen(buffer, "r") ;
  taxonomy.Load(fp) ;
  fclose(fp) ;

  if (inspectItem == ARGV_INSPECT_SEQNAME)
  {
  }
  else if (inspectItem == ARGV_INSPECT_SUMMARY)
  {
  }
  else if (inspectItem == ARGV_CONVERSION_TABLE)
  {
    taxonomy.PrintConversionTable(stdout) ;
  }
  else if (inspectItem == ARGV_TAXONOMY_TREE)
  {
    taxonomy.PrintTaxonomyTree(stdout) ;
  }
  else if (inspectItem == ARGV_NAME_TABLE)
  {
    taxonomy.PrintNameTable(stdout) ;
  }
  else
  {
    fprintf(stderr, "Use inspect options from %s", usage) ;
    return EXIT_FAILURE ;
  }

  free(idxPrefix) ;

  return 0 ;
}

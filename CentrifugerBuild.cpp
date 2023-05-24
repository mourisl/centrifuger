#include <stdio.h>
#include <time.h>
#include <getopt.h>

#include "argvdefs.h"
#include "ReadFiles.hpp"
#include "compactds/Sequence_Hybrid.hpp"
#include "compactds/FMBuilder.hpp"
#include "Taxonomy.hpp"

char usage[] = "./centrifuger-build [OPTIONS]:\n"
  "Required:\n"
  "\t-r FILE: reference sequence file\n"
  "\t--taxonomy-tree FILE: \n"
  "\t--name-table FILE: \n"
  "\t--conversion-table FILE: \n"
  "Optional:\n"
  "\t-o STRING: output prefix [centrifuger]\n"
  "\t-t INT: number of threads [1]\n"
  "\t--bmax INT: [1024]\n"
  "\t--offrate INT: [5]\n"
  ""
  ;

static const char *short_options = "r:o:t:" ;
static struct option long_options[] = {
			{ "bmax", required_argument, 0, ARGV_BMAX},
			{ "dcv", required_argument, 0, ARGV_DCV},
      { "offrate", required_argument, 0, ARGV_OFFRATE},
      { "taxonomy-tree", required_argument, 0, ARGV_TAXONOMY_TREE},
      { "conversion-table", required_argument, 0, ARGV_CONVERSION_TABLE},
			{ "name-table", required_argument, 0, ARGV_NAME_TABLE},
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
  char outputPrefix[1024] = "centrifuger" ;
  char outputFileName[1024] ;
  char *taxonomyFile = NULL ; // taxonomy tree file
  char *nameTable = NULL ;
  char *conversionTable = NULL ;
  int threadCnt = 1 ;
  ReadFiles refGenomeFile ;

  FMBuilder fmBuilder ;
  Taxonomy taxonomy ;

  while (1)
  {
		c = getopt_long( argc, argv, short_options, long_options, &option_index ) ;
		
		if (c == -1)
			break ;
  
    if (c == 'r') // reference genome file
    {
      refGenomeFile.AddReadFile(optarg, false) ;
    }
    else if (c == 'o')
    {
      strcpy(outputPrefix, optarg) ;
    }
    else if (c == 't')
    {
      threadCnt = atoi(optarg) ;
    }
    else if (c == ARGV_TAXONOMY_TREE)
    {
      taxonomyFile = strdup(optarg) ;  
    }
    else if (c == ARGV_NAME_TABLE)
    {
      nameTable = strdup(optarg) ;
    }
    else if (c == ARGV_CONVERSION_TABLE)
    {
      conversionTable = strdup(optarg) ;
    }
		else
		{
			fprintf( stderr, "Unknown parameter found\n%s", usage ) ;
			return EXIT_FAILURE ;
		}
  }

  if (!taxonomyFile)
  {
    fprintf(stderr, "Need to use --taxonomy-tree to specify the taxonomy tree.\n") ;
    return EXIT_FAILURE ;
  }
  if (!nameTable)
  {
    fprintf(stderr, "Need to use --name-table to specify taxonomy names.\n") ;
    return EXIT_FAILURE ;
  }
  if (!conversionTable)
  {
    fprintf(stderr, "Need to use --conversion-table to specify sequence id to taxonomy id mapping.\n") ;
    return EXIT_FAILURE ;
  }

  taxonomy.Init(taxonomyFile, nameTable, conversionTable) ;
  //FixedSizeElemArray genome ;
  //fmBuilder.Build(genome, genomeLength,)
  
  FILE *fpOutput = NULL ;
  sprintf(outputFileName, "%s.2.cfr", outputPrefix) ;
  fpOutput = fopen(outputFileName, "wb") ;
  taxonomy.Save(fpOutput) ;
  fclose(fpOutput) ;
  
  free(taxonomyFile) ;
  free(nameTable) ;
  free(conversionTable) ;

  return 0 ;
}

#include <stdio.h>
#include <time.h>
#include <getopt.h>

#include "argvdefs.h"
#include "Builder.hpp"

char usage[] = "./centrifuger-build [OPTIONS]:\n"
  "Required:\n"
  "\t-r FILE: reference sequence file (can use multiple -r to specify more than one input file)\n"
  "\t\tor\n"
  "\t-l FILE: list of reference sequence file stored in <file>, one file per row\n"
  "\t--taxonomy-tree FILE: taxonomy tree, i.e., nodes.dmp file\n"
  "\t--name-table FILE: name table, i.e., names.dmp file\n"
  "\t--conversion-table FILE: seqID to taxID conversion file\n"
  "Optional:\n"
  "\t-o STRING: output prefix [centrifuger]\n"
  "\t-t INT: number of threads [1]\n"
  "\t--build-mem STR: automatic infer bmax and dcv to match memory constraints, can use T,G,M,K to specify the memory size [not used]\n"
  "\t--bmax INT: block size for blockwise suffix array sorting [16777216]\n"
  "\t--offrate INT: SA/offset is sampled every (2^<int>) BWT chars [4]\n"
  "\t--dcv INT: difference cover period [4096]\n"
  "\t--subset-tax INT: only consider the subset of input genomes under taxonomy node INT [0]\n"
  ""
  ;

static const char *short_options = "r:l:o:t:" ;
static struct option long_options[] = {
      { "bmax", required_argument, 0, ARGV_BMAX},
			{ "dcv", required_argument, 0, ARGV_DCV},
      { "build-mem", required_argument, 0, ARGV_BUILD_MEMORY},
      { "offrate", required_argument, 0, ARGV_OFFRATE},
      { "taxonomy-tree", required_argument, 0, ARGV_TAXONOMY_TREE},
      { "conversion-table", required_argument, 0, ARGV_CONVERSION_TABLE},
			{ "name-table", required_argument, 0, ARGV_NAME_TABLE},
      { "subset-tax", required_argument, 0, ARGV_SUBSET_TAXONOMY}, 
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
  char *taxonomyFile = NULL ; // taxonomy tree file
  char *nameTable = NULL ;
  char *conversionTable = NULL ;
  uint64_t subsetTax = 0 ; 
  size_t buildMemoryConstraint = 0 ;
  ReadFiles refGenomeFile ;

  struct _FMBuilderParam fmBuilderParam ;
  fmBuilderParam.sampleRate = 16 ;

  while (1)
  {
		c = getopt_long( argc, argv, short_options, long_options, &option_index ) ;
		
		if (c == -1)
			break ;
  
    if (c == 'r') // reference genome file
    {
      refGenomeFile.AddReadFile(optarg, false) ;
    }
    else if (c == 'l')
    {
      char *fileName = (char *)malloc(sizeof(char) * 4096) ;
      FILE *fpList = fopen(optarg, "r") ;
      while (fscanf(fpList, "%s", fileName) != EOF)
        refGenomeFile.AddReadFile(fileName, false) ;
      fclose(fpList) ;
      free(fileName) ;
    }
    else if (c == 'o')
    {
      strcpy(outputPrefix, optarg) ;
    }
    else if (c == 't')
    {
      fmBuilderParam.threadCnt = atoi(optarg) ;
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
    else if (c == ARGV_DCV)
    {
      fmBuilderParam.saDcv = atoi(optarg) ;
    }
    else if (c == ARGV_BMAX)
    {
      fmBuilderParam.saBlockSize = atoi(optarg) ;
    }
    else if (c == ARGV_BUILD_MEMORY)
    {
      buildMemoryConstraint = Utils::SpaceStringToBytes(optarg) ;
    }
    else if (c == ARGV_OFFRATE)
    {
      fmBuilderParam.sampleRate = (1<<atoi(optarg)) ;
    }
    else if (c == ARGV_SUBSET_TAXONOMY)
    {
      sscanf(optarg, "%lu", &subsetTax) ;
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

  const char alphabetList[] = "ACGT" ;

  Builder builder ;
  builder.Build(refGenomeFile, taxonomyFile, nameTable, conversionTable, subsetTax, buildMemoryConstraint, fmBuilderParam, alphabetList) ;
  builder.Save(outputPrefix) ;

  free(taxonomyFile) ;
  free(nameTable) ;
  free(conversionTable) ;

  return 0 ;
}

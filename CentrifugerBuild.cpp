#include <stdio.h>
#include <time.h>
#include <getopt.h>

#include "argvdefs.h"
#include "ReadFiles.hpp"
#include "compactds/Sequence_Hybrid.hpp"
#include "compactds/FMBuilder.hpp"
#include "compactds/FMIndex.hpp"
#include "compactds/Alphabet.hpp"
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
  "\t--block INT: [16777216]\n"
  "\t--offrate INT: [5]\n"
  "\t--dcv INT: [4096]\n"
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
  ReadFiles refGenomeFile ;

  Alphabet alphabets ;
  Taxonomy taxonomy ;
  struct _FMBuilderParam fmBuilderParam ;

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
  
  const char alphabetList[] = "ACGT" ;
  const int alphabetSize = strlen(alphabetList) ;
  int alphabetCodeLen = alphabets.InitFromList(alphabetList, alphabetSize) ;

  FixedSizeElemArray genomes ;
  genomes.Malloc(alphabetCodeLen, 1000000) ;
  genomes.SetSize(0) ;
  std::map<size_t, size_t> seqLength ; // we use map here is for the case that a seq show up in the conversion table but not in the actual genome file.

  while (refGenomeFile.Next())
  {
    size_t seqid = taxonomy.SeqNameToId(refGenomeFile.id) ;
    if (seqid >= taxonomy.GetSeqCount())
    {
      fprintf(stderr, "WARNING: taxonomy id doesn't exist for %s!\n", refGenomeFile.id) ;
      seqid = taxonomy.AddExtraSeqName(refGenomeFile.id) ;
    }
    
    // Remove the Ns and convert lower-case sequences to upper-case
    size_t i, k ;
    char *s = refGenomeFile.seq ;
    k = 0 ;
    for (i = 0 ; s[i] ; ++i)
    {
      if (s[i] >= 'a' && s[i] <= 'z')
        s[i] = s[i] - 'a' + 'A' ;
      if (alphabets.IsIn(s[i]))
      {
        s[k] = s[i] ;
        ++k ;
      }
    }
    s[k] = '\0' ;
    for (i = 0 ; i < k ; ++i)
    {
      genomes.PushBack( alphabets.Encode(s[i])) ;   
    }
    seqLength[seqid] = k ;
  }
  
  FixedSizeElemArray BWT ;
  size_t firstISA ;
  struct _FMIndexAuxData fmAuxData ;
  FMBuilder::MallocAuxiliaryData(fmAuxData, alphabetCodeLen, genomes.GetSize(), fmBuilderParam) ;
  FMBuilder::Build(genomes, genomes.GetSize(), alphabetSize, BWT, firstISA, fmAuxData, fmBuilderParam) ;
  FMIndex<Sequence_Hybrid> fmIndex ;
  fmIndex.Init(BWT, genomes.GetSize(), 
      firstISA, fmAuxData, alphabetList, alphabetSize) ;
  
  // Convert the sampled point to seqID.
  FILE *fpOutput ;
  // .1.cfr file is for the index
  sprintf(outputFileName, "%s.1.cfr", outputPrefix) ;
  fpOutput = fopen(outputFileName, "w") ;
  fmIndex.Save(fpOutput) ;
  fclose(fpOutput) ;

  // .2.cfr file is for taxonomy structure
  sprintf(outputFileName, "%s.2.cfr", outputPrefix) ;
  fpOutput = fopen(outputFileName, "w") ;
  taxonomy.Save(fpOutput) ;
  fclose(fpOutput) ;

  // .3.cfr file is for sequence length
  sprintf(outputFileName, "%s.3.cfr", outputPrefix) ;
  fpOutput = fopen(outputFileName, "w") ;
  for (std::map<size_t, size_t>::iterator iter = seqLength.begin() ; 
      iter != seqLength.end() ; ++iter)
  {
    size_t tmp[2] = {iter->first, iter->second} ;
    fwrite(tmp, sizeof(tmp[0]), 2, fpOutput) ;
  }
  fclose(fpOutput) ;
  
  free(taxonomyFile) ;
  free(nameTable) ;
  free(conversionTable) ;

  return 0 ;
}

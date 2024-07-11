#include <stdio.h>
#include <getopt.h>

#include "defs.h"
#include "argvdefs.h"
#include "Taxonomy.hpp"
#include "compactds/FMIndex.hpp"
#include "compactds/Sequence_RunBlock.hpp"

char usage[] = "./centrifuger-inspect [OPTIONS]:\n"
  "Required:\n"
  "\t-x STRING: index prefix\n"
  "One of:\n"
  "\t--summary: print the summary information for each strain in the database\n"
  //"\t--name: print the scientific name \n"
  "\t--conversion-table: print the seqID to taxonomy ID translation information\n"
  "\t--taxonomy-tree: print the taxonomy tree\n"
  "\t--name-table: print the scientific name for each strain in the database\n"
  "\t--size-table: print the lengths of the sequences belonging to the same taxonomic ID"
  "\t--index-size: print the index information\n"
  ""
  ;

static const char *short_options = "x:" ;
static struct option long_options[] = {
  {"summary", no_argument, 0, ARGV_INSPECT_SUMMARY},
  {"seq-name", no_argument, 0, ARGV_INSPECT_SEQNAME},
  {"conversion-table", no_argument, 0, ARGV_CONVERSION_TABLE},
  {"taxonomy-tree", no_argument, 0, ARGV_TAXONOMY_TREE},
  {"name-table", no_argument, 0, ARGV_NAME_TABLE},
  {"size-table", no_argument, 0, ARGV_INSPECT_SIZE_TABLE},
  {"index-size", no_argument, 0, ARGV_INSPECT_INDEXSIZE},
  { (char *)0, 0, 0, 0} 
} ;

using namespace compactds ;

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

  std::map<size_t, size_t> seqLength ;
  sprintf(buffer, "%s.3.cfr", idxPrefix) ;
  fp = fopen(buffer, "r") ;
  size_t tmp[2] ;
  while (fread(tmp, sizeof(tmp[0]), 2, fp))
    seqLength[tmp[0]] = tmp[1] ;
  fclose(fp) ;

  if (inspectItem == ARGV_INSPECT_SEQNAME)
  {
    
  }
  else if (inspectItem == ARGV_INSPECT_SUMMARY)
  {
    for (std::map<size_t, size_t>::iterator iter = seqLength.begin() ;
        iter != seqLength.end() ; ++iter)
    {
      size_t ctid = taxonomy.SeqIdToTaxId(iter->first) ; 
      fprintf(stdout, "%s\t%lu\t%lu\t%s\n", 
          taxonomy.SeqIdToName(iter->first).c_str(), // sequence name
          taxonomy.GetOrigTaxId(ctid),
          iter->second, // sequence length
          taxonomy.GetTaxIdName(ctid).c_str()
          ) ;
    }
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
  else if (inspectItem == ARGV_INSPECT_SIZE_TABLE)
  {
    size_t i ;
    size_t n = taxonomy.GetNodeCount() ;
    size_t *taxidLength = (size_t *)malloc(sizeof(size_t) * n) ;
    taxonomy.ConvertSeqLengthToTaxLength(seqLength, taxidLength) ;
    for (i = 0 ; i < n ; ++i)
    {
      if (taxidLength[i] == 0)//|| i == taxonomy.GetRoot() )
        continue ;
      fprintf(stdout, "%lu\t%lu\n", taxonomy.GetOrigTaxId(i), taxidLength[i]) ;
    }
    free(taxidLength) ;
  }
  else if (inspectItem == ARGV_INSPECT_INDEXSIZE)
  {
    sprintf(buffer, "%s.1.cfr", idxPrefix) ; 
    FMIndex<Sequence_RunBlock> fm ;
    FILE *fp = fopen(buffer, "r") ;
    fm.Load(fp) ;
    fclose(fp) ;

    fm.PrintSpace() ;
  }
  else
  {
    fprintf(stderr, "Use inspect options from %s", usage) ;
    return EXIT_FAILURE ;
  }

  free(idxPrefix) ;

  return 0 ;
}

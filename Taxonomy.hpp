#ifndef _MOURISL_TAXONOMY_HEADER
#define _MOURISL_TAXONOMY_HEADER

// Partly based on the taxonomy.h file implemented by Florian Breitwieser in Centrifuge.

#include<map>
#include<utility>
#include<string>
#include <fstream>

#include <stdio.h> 

enum 
{
    RANK_UNKNOWN = 0,
    RANK_STRAIN,
    RANK_SPECIES,
    RANK_GENUS,
    RANK_FAMILY,
    RANK_ORDER,
    RANK_CLASS,
    RANK_PHYLUM,
    RANK_KINGDOM,
    RANK_DOMAIN,
    RANK_FORMA,
    RANK_INFRA_CLASS,
    RANK_INFRA_ORDER,
    RANK_PARV_ORDER,
    RANK_SUB_CLASS,
    RANK_SUB_FAMILY,
    RANK_SUB_GENUS,
    RANK_SUB_KINGDOM,
    RANK_SUB_ORDER,
    RANK_SUB_PHYLUM,
    RANK_SUB_SPECIES,
    RANK_SUB_TRIBE,
    RANK_SUPER_CLASS,
    RANK_SUPER_FAMILY,
    RANK_SUPER_KINGDOM,
    RANK_SUPER_ORDER,
    RANK_SUPER_PHYLUM,
    RANK_TRIBE,
    RANK_VARIETAS,
    RANK_LIFE,
    RANK_MAX
};

struct TaxonomyNode 
{
    uint64_t parentTid;
    uint8_t  rank;
    uint8_t  leaf;

    TaxonomyNode(uint64_t _parent_tid, uint8_t  _rank, uint8_t _leaf):
    	parentTid(_parent_tid), rank(_rank), leaf(_leaf) {};

    TaxonomyNode(): parentTid(0), rank(RANK_UNKNOWN), leaf(false) {};
};


class Taxonomy
{
private:
  std::map<uint64_t, uint64_t> _compactTaxId ;
  std::map<uint64_t, TaxonomyNode> _taxonomyTree;
  uint8_t _taxRankNum[RANK_MAX] ;

  void InitTaxRankNum()
  {
    uint8_t rank = 0;

    _taxRankNum[RANK_SUB_SPECIES] = rank;
    _taxRankNum[RANK_STRAIN] = rank++;

    _taxRankNum[RANK_SPECIES] = rank++;

    _taxRankNum[RANK_SUB_GENUS] = rank;
    _taxRankNum[RANK_GENUS] = rank++;

    _taxRankNum[RANK_SUB_FAMILY] = rank;
    _taxRankNum[RANK_FAMILY] = rank;
    _taxRankNum[RANK_SUPER_FAMILY] = rank++;

    _taxRankNum[RANK_SUB_ORDER] = rank;
    _taxRankNum[RANK_INFRA_ORDER] = rank;
    _taxRankNum[RANK_PARV_ORDER] = rank;
    _taxRankNum[RANK_ORDER] = rank;
    _taxRankNum[RANK_SUPER_ORDER] = rank++;

    _taxRankNum[RANK_INFRA_CLASS] = rank;
    _taxRankNum[RANK_SUB_CLASS] = rank;
    _taxRankNum[RANK_CLASS] = rank;
    _taxRankNum[RANK_SUPER_CLASS] = rank++;

    _taxRankNum[RANK_SUB_PHYLUM] = rank;
    _taxRankNum[RANK_PHYLUM] = rank;
    _taxRankNum[RANK_SUPER_PHYLUM] = rank++;

    _taxRankNum[RANK_SUB_KINGDOM] = rank;
    _taxRankNum[RANK_KINGDOM] = rank;
    _taxRankNum[RANK_SUPER_KINGDOM] = rank++;

    _taxRankNum[RANK_DOMAIN] = rank;
    _taxRankNum[RANK_FORMA] = rank;
    _taxRankNum[RANK_SUB_TRIBE] = rank;
    _taxRankNum[RANK_TRIBE] = rank;
    _taxRankNum[RANK_VARIETAS] = rank;
    _taxRankNum[RANK_UNKNOWN] = rank;
  }

  
  void ReadTaxonomyTree(string taxonomy_fname) 
  {
    ifstream taxonomy_file(taxonomy_fname.c_str(), ios::in);
    if(taxonomy_file.is_open()) {
      char line[1024];
      while(!taxonomy_file.eof()) {
        line[0] = 0;
        taxonomy_file.getline(line, sizeof(line));
        if(line[0] == 0 || line[0] == '#') continue;
        istringstream cline(line);
        uint64_t tid, parent_tid;
        char dummy; string rank_string;
        cline >> tid >> dummy >> parent_tid >> dummy >> rank_string;
        if(_taxonomyTree.find(tid) != _taxonomyTree.end()) {
          cerr << "Warning: " << tid << " already has a parent!" << endl;
          continue;
        }

        _taxonomyTree[tid] = TaxonomyNode(parent_tid, GetTaxRankId(rank_string.c_str()), false);
      }
      taxonomy_file.close();
    } else {
      cerr << "Error: " << taxonomy_fname << " doesn't exist!" << endl;
      throw 1;
    }
  }

  void ReadTaxonomyName(string fname)
  {
    
  }

public:
  Taxonomy() {}
  ~Taxonomy() {}
  
  void Init(char *nodesFile, char *namesFile)
  {
    InitTaxRankNum() ;

    const int bufferSize = 4096 ;
    char *bufferSize = (char *)malloc(sizeof(char) * bufferSize) ;
    
    FILE *fp ;

    fp = fopen(nodesFile, "r") ;
    fclose(fp) ;
    
    fp = fopen(nodesFile, "r") ;
    fclose(fp) ;
  }

  const char *GetTaxRankString(uint8_t rank)
  {
    switch(rank) {
      case RANK_STRAIN:        return "strain";
      case RANK_SPECIES:       return "species";
      case RANK_GENUS:         return "genus";
      case RANK_FAMILY:        return "family";
      case RANK_ORDER:         return "order";
      case RANK_CLASS:         return "class";
      case RANK_PHYLUM:        return "phylum";
      case RANK_KINGDOM:       return "kingdom";
      case RANK_FORMA:         return "forma";
      case RANK_INFRA_CLASS:   return "infraclass";
      case RANK_INFRA_ORDER:   return "infraorder";
      case RANK_PARV_ORDER:    return "parvorder";
      case RANK_SUB_CLASS:     return "subclass";
      case RANK_SUB_FAMILY:    return "subfamily";
      case RANK_SUB_GENUS:     return "subgenus";
      case RANK_SUB_KINGDOM:   return "subkingdom";
      case RANK_SUB_ORDER:     return "suborder";
      case RANK_SUB_PHYLUM:    return "subphylum";
      case RANK_SUB_SPECIES:   return "subspecies";
      case RANK_SUB_TRIBE:     return "subtribe";
      case RANK_SUPER_CLASS:   return "superclass";
      case RANK_SUPER_FAMILY:  return "superfamily";
      case RANK_SUPER_KINGDOM: return "superkingdom";
      case RANK_SUPER_ORDER:   return "superorder";
      case RANK_SUPER_PHYLUM:  return "superphylum";
      case RANK_TRIBE:         return "tribe";
      case RANK_VARIETAS:      return "varietas";
      case RANK_LIFE:          return "life";
      default:                 return "no rank";
    };
  }

  uint8_t GetTaxRankId(const char* rank) 
  {
    if(strcmp(rank, "strain") == 0) {
      return RANK_STRAIN;
    } else if(strcmp(rank, "species") == 0) {
      return RANK_SPECIES;
    } else if(strcmp(rank, "genus") == 0) {
      return RANK_GENUS;
    } else if(strcmp(rank, "family") == 0) {
      return RANK_FAMILY;
    } else if(strcmp(rank, "order") == 0) {
      return RANK_ORDER;
    } else if(strcmp(rank, "class") == 0) {
      return RANK_CLASS;
    } else if(strcmp(rank, "phylum") == 0) {
      return RANK_PHYLUM;
    } else if(strcmp(rank, "kingdom") == 0) {
      return RANK_KINGDOM;
    } else if(strcmp(rank, "forma") == 0) {
      return RANK_FORMA;
    } else if(strcmp(rank, "infraclass") == 0) {
      return RANK_INFRA_CLASS;
    } else if(strcmp(rank, "infraorder") == 0) {
      return RANK_INFRA_ORDER;
    } else if(strcmp(rank, "parvorder") == 0) {
      return RANK_PARV_ORDER;
    } else if(strcmp(rank, "subclass") == 0) {
      return RANK_SUB_CLASS;
    } else if(strcmp(rank, "subfamily") == 0) {
      return RANK_SUB_FAMILY;
    } else if(strcmp(rank, "subgenus") == 0) {
      return RANK_SUB_GENUS;
    } else if(strcmp(rank, "subkingdom") == 0) {
      return RANK_SUB_KINGDOM;
    } else if(strcmp(rank, "suborder") == 0) {
      return RANK_SUB_ORDER;
    } else if(strcmp(rank, "subphylum") == 0) {
      return RANK_SUB_PHYLUM;
    } else if(strcmp(rank, "subspecies") == 0) {
      return RANK_SUB_SPECIES;
    } else if(strcmp(rank, "subtribe") == 0) {
      return RANK_SUB_TRIBE;
    } else if(strcmp(rank, "superclass") == 0) {
      return RANK_SUPER_CLASS;
    } else if(strcmp(rank, "superfamily") == 0) {
      return RANK_SUPER_FAMILY;
    } else if(strcmp(rank, "superkingdom") == 0) {
      return RANK_SUPER_KINGDOM;
    } else if(strcmp(rank, "superorder") == 0) {
      return RANK_SUPER_ORDER;
    } else if(strcmp(rank, "superphylum") == 0) {
      return RANK_SUPER_PHYLUM;
    } else if(strcmp(rank, "tribe") == 0) {
      return RANK_TRIBE;
    } else if(strcmp(rank, "varietas") == 0) {
      return RANK_VARIETAS;
    } else if(strcmp(rank, "life") == 0) {
      return RANK_LIFE;
    } else {
      return RANK_UNKNOWN;
    }
  }

  uint64_t GetTaxIdAtParentRank(const TaxonomyTree& tree, uint64_t taxid, uint8_t at_rank) 
  {
    while (true) {
      TaxonomyTree::const_iterator itr = tree.find(taxid);
      if(itr == tree.end()) {
        break;
      }
      const TaxonomyNode& node = itr->second;

      if (node.rank == at_rank) {
        return taxid;
      } else if (node.rank > at_rank || node.parent_tid == taxid) {
        return 0;
      }

      taxid = node.parent_tid;
    }
    return 0;
  }


  void SaveBinary(FILE *fp)
  {
  }

  void LoadBinary(FILE *fp)
  {
  }
} ;

#endif

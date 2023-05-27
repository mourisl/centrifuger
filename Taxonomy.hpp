#ifndef _MOURISL_TAXONOMY_HEADER
#define _MOURISL_TAXONOMY_HEADER

// Partly based on the taxonomy.h file implemented by Florian Breitwieser in Centrifuge.
// This class handles the taxonomy-related information, including taxonomy tree and taxonomy id mappings 

#include<map>
#include<utility>
#include<string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <fstream>
#include <fstream>

#include <stdio.h> 

#include "MapID.hpp"

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
  MapID<uint64_t> _taxIdMap ;
  struct TaxonomyNode *_taxonomyTree; // Use arrays to hold the taxonomy information, more efficient to access
  std::string *_taxonomyName ;
  MapID<std::string> _seqStrNameMap ; 
  uint64_t *_seqIdToTaxId ;

  size_t _nodeCnt ;
  size_t _seqCnt ; // the sequences with taxonomy information
  size_t _extraSeqCnt ; // the number of sequences 
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
  
  void ReadTaxonomyTree(std::string taxonomy_fname, std::map<uint64_t, int> &presentTax) 
  {
    std::ifstream taxonomy_file(taxonomy_fname.c_str(), std::ios::in);
    std::map<uint64_t, struct TaxonomyNode> tree ;
    if(taxonomy_file.is_open()) {
      char line[1024];
      while(!taxonomy_file.eof()) {
        line[0] = 0;
        taxonomy_file.getline(line, sizeof(line));
        if(line[0] == 0 || line[0] == '#') continue;
        std::istringstream cline(line);
        uint64_t tid, parent_tid;
        char dummy; std::string rank_string;
        cline >> tid >> dummy >> parent_tid >> dummy >> rank_string;
        if(tree.find(tid) != tree.end()) {
          std::cerr << "Warning: " << tid << " already has a parent!" << std::endl;
          continue;
        }
        tree[tid] = TaxonomyNode(parent_tid, GetTaxRankId(rank_string.c_str()), true);
      }
      taxonomy_file.close();
    } else {
      std::cerr << "Error: " << taxonomy_fname << " doesn't exist!" << std::endl;
      throw 1;
    }

    // Get the parent nodes related to the present leaf tax nodes
    std::map<uint64_t, int> selectedTax ;  
    for (std::map<uint64_t, int>::iterator iter = presentTax.begin() ; iter != presentTax.end(); ++iter)
    {
      if (tree.find(iter->first) ==tree.end()) {
        std::cerr << "Warning: " << iter->first << " is not in the taxonomy tree" << std::endl;
        continue;
      }
      uint64_t p = iter->first ;
      while (1)
      {
        if (selectedTax.find(p) != selectedTax.end())
          break ;
        selectedTax[p] = 1;
        p = tree[p].parentTid;  
      }
    }
    presentTax = selectedTax;

    // Clean up the tree 
    std::map<uint64_t, struct TaxonomyNode> cleanTree;
    for (std::map<uint64_t, struct TaxonomyNode>::iterator iter = tree.begin(); iter != tree.end(); ++iter)
    {
      if (selectedTax.find(iter->first) == selectedTax.end())
        continue ;
      cleanTree[iter->first] = tree[iter->first] ;
      _taxIdMap.Add(iter->first) ;
    }

    // Flatten the taxonomy tree to the array
    _nodeCnt = _taxIdMap.GetSize() ;
    _taxonomyTree = new struct TaxonomyNode[_nodeCnt] ;
    memset(_taxonomyTree, 0, sizeof(TaxonomyNode) * _nodeCnt) ;
    for (std::map<uint64_t, struct TaxonomyNode>::iterator it = cleanTree.begin() ;
        it != cleanTree.end() ; ++it)
    {
      uint64_t i = _taxIdMap.Map(it->first) ; 
      _taxonomyTree[i] = it->second ;
    }
    
    // We need to split the update parent node here because the order is random.
    for (uint64_t i = 0 ; i < _nodeCnt ; ++i)
    {
      _taxonomyTree[i].parentTid = _taxIdMap.Map(_taxonomyTree[i].parentTid) ;
      _taxonomyTree[ _taxonomyTree[i].parentTid ].leaf = false ;
    }

  }

  void ReadTaxonomyName(std::string fname, std::map<uint64_t, int> &presentTax)
  {
    std::ifstream taxname_file(fname.c_str(), std::ios::in);
    _taxonomyName = new std::string[_nodeCnt] ;
    if(taxname_file.is_open()) {
      char line[1024];
      while(!taxname_file.eof()) {
        line[0] = 0;
        taxname_file.getline(line, sizeof(line));
        if(line[0] == 0 || line[0] == '#') continue;
        if(!strstr(line, "scientific name")) continue;  
        std::istringstream cline(line);
        uint64_t tid;
        char dummy; std::string scientific_name;
        cline >> tid >> dummy >> scientific_name ;
        if (presentTax.find(tid) == presentTax.end())
          continue ;
        uint64_t ctid = _taxIdMap.Map(tid) ;// compact tid
        std::string temp;
        
        while(true) {
          cline >> temp;
          if(temp == "|") break;
          scientific_name.push_back('_');
          scientific_name += temp;
        }
        
        _taxonomyName[ctid] = scientific_name ;
      }
      taxname_file.close();
    } else {
      std::cerr << "Error: " << fname << " doesn't exist!" << std::endl;
      throw 1;
    }
  }
  
  void ReadPresentTaxonomyLeafs(std::string fname, std::map<uint64_t, int> &presentLeafs)
  {
    std::ifstream seqmap_file(fname.c_str(), std::ios::in);
    std::map<std::string, uint64_t> rawSeqNameMap ;
    if(seqmap_file.is_open()) {
      char line[1024];
      while(!seqmap_file.eof()) {
        line[0] = 0;
        seqmap_file.getline(line, sizeof(line));
        if(line[0] == 0 || line[0] == '#') continue;
        std::istringstream cline(line);
        uint64_t tid;
        std::string seqId;
        cline >> seqId >> tid ;
        presentLeafs[tid] = 0;
      }
      seqmap_file.close();
    } else {
      std::cerr << "Error: " << fname << " doesn't exist!" << std::endl;
      throw 1;
    }
  }

  void ReadSeqNameFile(std::string fname)
  {
    std::ifstream seqmap_file(fname.c_str(), std::ios::in);
    std::map<std::string, uint64_t> rawSeqNameMap ;
    if(seqmap_file.is_open()) {
      char line[1024];
      while(!seqmap_file.eof()) {
        line[0] = 0;
        seqmap_file.getline(line, sizeof(line));
        if(line[0] == 0 || line[0] == '#') continue;
        std::istringstream cline(line);
        uint64_t tid;
        std::string seqIdStr;
        cline >> seqIdStr >> tid ;
        _seqStrNameMap.Add(seqIdStr) ;
        rawSeqNameMap[seqIdStr] = tid ;
      }
      seqmap_file.close();
    } else {
      std::cerr << "Error: " << fname << " doesn't exist!" << std::endl;
      throw 1;
    }
    
    // Map sequence string identifier to compact taxonomy id
    _seqIdToTaxId = new uint64_t[ _seqStrNameMap.GetSize() ] ;
    for (std::map<std::string, uint64_t>::iterator iter = rawSeqNameMap.begin() ;
        iter != rawSeqNameMap.end() ; ++iter)
    {
     _seqIdToTaxId[ _seqStrNameMap.Map(iter->first) ] = _taxIdMap.Map(iter->second) ; 
    }
    _seqCnt = _seqStrNameMap.GetSize() ;
  }
  
public:
  Taxonomy() 
  {
    _taxonomyTree = NULL ;
    _taxonomyName = NULL ;
    _seqIdToTaxId = NULL ;
    _nodeCnt = 0 ;
    _seqCnt = 0 ;
    _extraSeqCnt = 0 ;
    InitTaxRankNum() ;
  }

  ~Taxonomy() 
  {
    Free() ;
  }

  void Free()
  {
    if (_nodeCnt > 0)
    {
      if (_taxonomyTree != NULL)
        delete[] _taxonomyTree ;
      if (_taxonomyName != NULL)
        delete[] _taxonomyName ;
      if (_seqIdToTaxId != NULL)
        delete[] _seqIdToTaxId ;
      _nodeCnt = 0 ;
    }
  }

  void Init(const char *nodesFile, const char *namesFile, const char *seqIdFile)
  {
    std::map<uint64_t, int> presentTax;
    ReadPresentTaxonomyLeafs(std::string(seqIdFile), presentTax) ;
    ReadTaxonomyTree(std::string(nodesFile), presentTax) ;
    ReadTaxonomyName(std::string(namesFile), presentTax) ;
    ReadSeqNameFile(std::string(seqIdFile)) ;
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

  // Also returns compact tax id
  uint64_t GetTaxIdAtParentRank(uint64_t taxid, uint8_t at_rank) 
  {
    while (true) {
      const TaxonomyNode& node = _taxonomyTree[taxid] ;

      if (node.rank == at_rank) {
        return taxid;
      } else if (node.rank > at_rank || node.parentTid == taxid) {
        return 0;
      }

      taxid = node.parentTid;
    }
    return 0;
  }

  uint64_t GetNodeCount()
  {
    return _nodeCnt ;
  }

  uint64_t GetSeqCount()
  {
    return _seqCnt ;
  }
  
  uint64_t GetOrigTaxId(uint64_t taxid)
  {
    return _taxIdMap.Inverse(taxid) ;
  }

  const char *GetName(uint64_t ctid) // compact taxtonomy id
  {
    return _taxonomyName[ctid].c_str() ;
  }

  uint64_t SeqNameToId(std::string &s)
  {
    if (!_seqStrNameMap.IsIn(s))
      return _seqStrNameMap.GetSize() ; 
    else
      return _seqStrNameMap.Map(s) ;
  }
  
  uint64_t SeqNameToId(const char *s)
  {
    std::string tmps(s) ;
    return SeqNameToId(tmps) ;
  }
 
  // Directly add a seqId(string)
  void AddExtraSeqName(char *s)
  {
    _seqStrNameMap.Add(s) ;
    ++_extraSeqCnt ;
  }

  uint64_t SeqIdToTaxId(uint64_t seqId)
  {
    if (seqId < _seqCnt)
      return _seqIdToTaxId[seqId] ;
    else
      return 0 ;
  }

  void Save(std::ofstream &ofs)
  {
    ofs << _nodeCnt << _seqCnt << _extraSeqCnt ;
    // Save the taxnomoy information
    ofs.write((char *)_taxonomyTree, sizeof(_taxonomyTree[0]) * _nodeCnt) ;
    _taxIdMap.Save(ofs) ;
    size_t i ;
    for (i = 0 ; i < _nodeCnt ; ++i)
      ofs << _taxonomyName[i] ;
    
    // Save the seqID information
    ofs.write((char *)_seqIdToTaxId, sizeof(_seqIdToTaxId[0]) * _seqCnt) ;
    for (i = 0 ; i < _seqCnt + _extraSeqCnt ; ++i)
    {
      ofs << _seqStrNameMap.Inverse(i) ;
    }
  }

  void Load(std::ifstream &ifs)
  {
    Free() ;
    InitTaxRankNum() ;
    
    ifs >> _nodeCnt >> _seqCnt >> _extraSeqCnt ;
    
    // Load the taxnomoy information
    _taxonomyTree = new struct TaxonomyNode[_nodeCnt] ;
    ifs.read((char *)_taxonomyTree, sizeof(_taxonomyTree[0]) * _nodeCnt) ;
    _taxIdMap.Load(ifs) ;
    size_t i ;
    _taxonomyName = new std::string[_nodeCnt] ;
    for (i = 0 ; i < _nodeCnt ; ++i)
      ifs >> _taxonomyName[i] ;

    // Load the seqID information
    _seqIdToTaxId = new uint64_t[_seqCnt] ;
    ifs.read((char *)_seqIdToTaxId, sizeof(_seqIdToTaxId[0]) *_seqCnt) ;
    _seqStrNameMap.Clear() ;
    for (i = 0 ; i < _seqCnt + _extraSeqCnt ; ++i)
    {
      std::string seqStr ;
      ifs >> seqStr ;
      _seqStrNameMap.Add(seqStr) ;
    }
  }
} ;

#endif

#ifndef _MOURISL_TAXONOMY_HEADER
#define _MOURISL_TAXONOMY_HEADER

// Partly based on the taxonomy.h file implemented by Florian Breitwieser in Centrifuge.
// This class handles the taxonomy-related information, including taxonomy tree and taxonomy id mappings 

#include <map>
#include<utility>
#include<string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <stdio.h> 
#include <string.h>

#include "MapID.hpp"
#include "compactds/SimpleVector.hpp"
#include "compactds/Utils.hpp"
#include "compactds/Tree_Plain.hpp"

using namespace compactds ;

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
  size_t _rootCTaxId ;// the compact tax Id for the root

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
    _taxRankNum[RANK_LIFE] = rank;
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
    // The set 0 should be handled by the constructor now.
    //memset(_taxonomyTree, 0, sizeof(TaxonomyNode) * _nodeCnt) ;
    for (std::map<uint64_t, struct TaxonomyNode>::iterator it = cleanTree.begin() ;
        it != cleanTree.end() ; ++it)
    {
      uint64_t i = _taxIdMap.Map(it->first) ; 
      _taxonomyTree[i] = it->second ;
    }
    
    // We need to split the update parent node here because the order is random.
    for (uint64_t i = 0 ; i < _nodeCnt ; ++i)
    {
      if (_taxIdMap.IsIn(_taxonomyTree[i].parentTid))
      {
        _taxonomyTree[i].parentTid = _taxIdMap.Map(_taxonomyTree[i].parentTid) ;
        _taxonomyTree[ _taxonomyTree[i].parentTid ].leaf = false ;
      }
      else
      {
        Utils::PrintLog("WARNING: parent tax ID of %lu does not exist. Set its parent to itself.", GetOrigTaxId(i)) ;
        _taxonomyTree[i].parentTid = i ;
      }
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

  void ReadSeqNameFile(std::string fname, bool conversionTableAtFileLevel)
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
        if (conversionTableAtFileLevel)
        {
          char buffer[1024] ;
          Utils::GetFileBaseName(seqIdStr.c_str(), "fna|fa|fasta|faa", buffer) ;
          seqIdStr = buffer ;
        }
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

  // Whether b is next to a in accession id
  bool IsNextSeqName(const char *a, const char *b)
  {
    int i, j ;
    uint64_t id[2] ;
    for (i = 0 ; i < 2 ; ++i)
    {
      const char *s = a ;
      if (i == 1)  
        s = b ;
      id[i] = 0 ;
      for (j = 0 ; s[j] ; ++j)
        if (s[j] >= '0' && s[j] <= '9')
          break ;
      
      //if (j > 2) // It's not something like GCFXXXX numbering
      //  return false ;

      for (; s[j] ; ++j)
      {
        if (s[j] >= '0' && s[j] <= '9')
        {
          id[i] = id[i] * 10 + s[j] - '0' ;
        }
        else
          break ;
      }
    }
    if (id[1] == id[0] + 1)
      return true ;
    else
      return false ;
  }
 
  void SaveString(FILE *fp, std::string &s)
  {
    size_t len = s.length() ;
    fwrite(&len, sizeof(len), 1, fp) ;
    fwrite(s.c_str(), sizeof(char), len, fp) ;
  }

  void LoadString(FILE *fp, std::string &s)
  {
    size_t len ;
    fread(&len, sizeof(len), 1, fp) ;
    char *buffer = (char *)malloc(sizeof(char) * (len + 1)) ;
    fread(buffer, sizeof(char), len, fp) ;
    buffer[len] = '\0' ;
    s = buffer ;
    free(buffer) ;
  }
  
  size_t FindRoot()
  {
    size_t i ;
    for (i = 0 ; i < _nodeCnt ; ++i)
      if (_taxonomyTree[i].parentTid == i)
        return i ;
    return _nodeCnt ;
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
    _rootCTaxId = 0 ;
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
  
  void Init(const char *nodesFile, const char *namesFile, const char *seqIdFile, bool conversionTableAtFileLevel)
  {
    std::map<uint64_t, int> presentTax;
    ReadPresentTaxonomyLeafs(std::string(seqIdFile), presentTax) ;
    ReadTaxonomyTree(std::string(nodesFile), presentTax) ;
    ReadTaxonomyName(std::string(namesFile), presentTax) ;
    ReadSeqNameFile(std::string(seqIdFile), conversionTableAtFileLevel) ;
  
    _rootCTaxId = FindRoot() ;
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
  size_t GetTaxIdAtParentRank(uint64_t taxid, uint8_t at_rank) 
  {
    while (true) {
      const TaxonomyNode& node = _taxonomyTree[taxid] ;

      if (node.rank == at_rank) {
        return taxid;
      } else if (node.rank > at_rank || node.parentTid == taxid) {
        return _nodeCnt ;
      }

      taxid = node.parentTid;
    }
    return _nodeCnt;
  }

  size_t GetNodeCount()
  {
    return _nodeCnt ;
  }

  size_t GetSeqCount()
  {
    return _seqCnt ;
  }

  size_t GetAllSeqCount()
  {
    return _seqCnt + _extraSeqCnt ;
  }
  
  uint64_t GetOrigTaxId(size_t taxid)
  {
    if (taxid >= _nodeCnt)
      return _taxIdMap.Inverse(_rootCTaxId) ;
    else
      return _taxIdMap.Inverse(taxid) ;
  }

  uint64_t GetRoot()
  {
    return _rootCTaxId ;
  }

  size_t CompactTaxId(uint64_t taxid)
  {
    if (_taxIdMap.IsIn(taxid))
      return _taxIdMap.Map(taxid) ;
    else
      return _nodeCnt ;
  }
  
  uint64_t GetParentTid(uint64_t ctid)
  {
    return _taxonomyTree[ctid].parentTid ;
  }

  uint8_t GetTaxIdRank(size_t ctid)
  {
    if (ctid >= _nodeCnt)
      return RANK_UNKNOWN ;
    else
      return _taxonomyTree[ctid].rank ;
  }

  std::string GetTaxIdName(size_t ctid)
  {
    if (ctid < _nodeCnt)
      return _taxonomyName[ctid] ;
    else
    {
      std::string tmp("Unknown") ;
      return tmp ;
    }
  }

  size_t SeqNameToId(std::string &s)
  {
    if (!_seqStrNameMap.IsIn(s))
      return _seqStrNameMap.GetSize() ; 
    else
      return _seqStrNameMap.Map(s) ;
  }
  
  size_t SeqNameToId(const char *s)
  {
    std::string tmps(s) ;
    return SeqNameToId(tmps) ;
  }

  std::string SeqIdToName(size_t seqid)
  {
    return _seqStrNameMap.Inverse(seqid) ;
  }
 
  // Directly add a seqId(string)
  // @return: id ;
  size_t AddExtraSeqName(char *s)
  {
    size_t ret = _seqStrNameMap.Add(s) ;
    ++_extraSeqCnt ;
    return ret ; 
  }

  uint64_t SeqIdToTaxId(size_t seqId)
  {
    if (seqId < _seqCnt)
      return _seqIdToTaxId[seqId] ;
    else
      return _nodeCnt ;
  }

  // Get the seq names
  void GetSeqNames(std::vector<std::string> &seqNames)
  {
    _seqStrNameMap.GetElemList(seqNames) ;
  }

  // Promote the tax id to higher level until number of taxids <= k, or reach LCA 
  void ReduceTaxIds(const SimpleVector<size_t> &taxIds, SimpleVector<size_t> &promotedTaxIds, int k)
  {
    int i ;
    int taxCnt = taxIds.Size() ;
    std::vector< SimpleVector<size_t> > taxPaths ;
    promotedTaxIds.Clear() ;
    if (taxIds.Size() <= k)
    {
      promotedTaxIds = taxIds ;
      return ;
    }

    // If there is a tax id not in the tree, we 
    //   give it no rank directly.
    for (i = 0 ; i < taxCnt ; ++i)
      if (taxIds[i] >= _nodeCnt)
      {
        promotedTaxIds.PushBack(_nodeCnt) ;
        return ;
      }
    // For each tax level, collect the found tax id on this level 
    std::map<size_t, int> taxIdsInRankNum[RANK_MAX] ;
    for (i = 0 ; i < taxCnt ; ++i)
    {
      size_t t = taxIds[i]; 
      uint8_t prevRankNum = 0 ;
      uint8_t ri ;// rank index

      taxIdsInRankNum[prevRankNum][t] = 1 ; // the input is at the base level
      do
      {
        uint8_t rankNum = _taxRankNum[_taxonomyTree[t].rank] ;
        if (rankNum != _taxRankNum[RANK_UNKNOWN] && rankNum > prevRankNum)
        {
          // Handle the case of missing taxonomy level in between
          for (ri = rankNum - 1 ; ri > prevRankNum ; --ri)
            taxIdsInRankNum[ri][t] = 1 ;

          if (taxIdsInRankNum[rankNum].find(t) == taxIdsInRankNum[rankNum].end())
            taxIdsInRankNum[rankNum][t] = 1 ;
          else
            break ; // the upper tax id has already been added, so no need to process anymore
          prevRankNum = rankNum ;
        }
        t = _taxonomyTree[t].parentTid ; 
      } while (t != _taxonomyTree[t].parentTid) ;
    }

    // Go through the levels until the tax ids <= k
    uint8_t ri ;
    for (ri = 0 ; ri < _taxRankNum[RANK_UNKNOWN] ; ++ri)
      if ((int)taxIdsInRankNum[ri].size() <= k)
        break ;
    
    for (std::map<size_t, int>::iterator iter = taxIdsInRankNum[ri].begin() ;
        iter != taxIdsInRankNum[ri].end() ; ++iter)
      promotedTaxIds.PushBack(iter->first) ;
    if (promotedTaxIds.Size() == 0)
      promotedTaxIds.PushBack(_rootCTaxId) ;
  }

  // @return: Number of children tax ids. childrenTax: compact tax ids below or equal to ctid.
  size_t GetChildrenTax(size_t ctid, std::map<size_t, int> &childrenTax)
  {
    childrenTax.clear() ;
    if (ctid >= _nodeCnt)
      return 0 ;

    size_t i, j ;
    size_t t ;
    int *visited ; // -1: not visited, 0: not a children, 1: is a children
    visited = (int *)malloc(sizeof(int) * _nodeCnt) ;
    memset(visited, -1, sizeof(visited[0]) * _nodeCnt) ;
    
    visited[ctid] = 1 ;
    SimpleVector<size_t> path ;
    for (i = 0 ; i < _nodeCnt ; ++i)
    {
      t = i ;
      path.Clear() ;
      
      while (t != _taxonomyTree[t].parentTid) // It's fine to put the while before
                          // we only need to add root if ctid is the root
      {
        if (visited[t] != -1)
          break ;
        path.PushBack(t) ;
        t = _taxonomyTree[t].parentTid ;
      }
      size_t pathSize = path.Size() ;
      int res = visited[t] ;
      if (res == -1)
        res = 0 ;
      
      for (j = 0 ; j < pathSize ; ++j)
        visited[ path[j] ] = res ;       
    }

    for (i = 0 ; i < _nodeCnt ; ++i)
    {
      if (visited[i] == 1)
        childrenTax[i] = 1 ;
    }
    free(visited) ;
  
    return childrenTax.size() ;
  }

  // Convert the taxonomy tree structure to a general tree that supports children operations.
  //   also the converted tree make sure every non-root node has a parent
  size_t ConvertToGeneralTree(Tree_Plain &tree)
  {
    size_t i ;
    tree.SetRoot(_rootCTaxId) ;
    tree.Init(_nodeCnt) ;
    for (i = 0 ; i < _nodeCnt ; ++i) 
    {
      if (i != GetParentTid(i))
        tree.AddEdge(i, GetParentTid(i)) ;
    }

    // Connect the disjoint trees to the root 
    std::vector<size_t> rootChildrenList = tree.GetChildren( tree.Root() ) ;
    std::map<size_t, int> rootChildrenMap ;
    size_t rootChildrenListSize = rootChildrenList.size() ;
    for (i = 0 ; i < rootChildrenListSize ; ++i)
      rootChildrenMap[ rootChildrenList[i] ] = 1 ;
    for (i = 0 ; i < _nodeCnt ; ++i)
      if (tree.Parent(i) == tree.Root() && rootChildrenMap.find(i) == rootChildrenMap.end())
        tree.AddEdge(i, tree.Root()) ;
  
    return _nodeCnt ;
  }

  // Assume taxIdLength is allocated of size _nodeCnt
  void ConvertSeqLengthToTaxLength(std::map<size_t, size_t> seqLength, size_t *taxidLength)
  {
    size_t i, j ;
    std::vector<std::string> seqNames ;
    GetSeqNames(seqNames) ;

    std::sort(seqNames.begin(), seqNames.end()) ;

    for (i = 0 ; i < _nodeCnt ; ++i)
      taxidLength[i] = 0 ;

    size_t seqNameCnt = seqNames.size() ;
    for (i = 0 ; i < seqNameCnt ; )
    {
      size_t seqId = SeqNameToId(seqNames[i]) ;
      size_t len = seqLength[seqId] ;
      uint64_t taxid = SeqIdToTaxId(seqId) ;
      for (j = i + 1 ; j < seqNameCnt ; ++j)
      {
        size_t nextSeqId = SeqNameToId(seqNames[j]) ;
        if (SeqIdToTaxId(nextSeqId) != taxid
            || !IsNextSeqName(seqNames[j - 1].c_str(),
              seqNames[j].c_str()))  
          break ;
        len += seqLength[nextSeqId] ;
      }

      if (taxid < _nodeCnt)
      {
        if (len > taxidLength[taxid])
          taxidLength[taxid] = len ;
      }
      //else // ignore the case the sequence is not in the tree.
      //  taxidLength[taxid] += len ;

      i = j ;
    }
  }

  void Save(FILE *fp)
  {
    SAVE_VAR(fp, _nodeCnt) ;
    SAVE_VAR(fp, _seqCnt) ;
    SAVE_VAR(fp, _extraSeqCnt) ;
    // Save the taxnomoy information
    fwrite(_taxonomyTree, sizeof(_taxonomyTree[0]), _nodeCnt, fp) ;
    _taxIdMap.Save(fp) ;
    size_t i ;
    for (i = 0 ; i < _nodeCnt ; ++i)
      SaveString(fp, _taxonomyName[i]) ;
    
    // Save the seqID information
    fwrite(_seqIdToTaxId, sizeof(_seqIdToTaxId[0]), _seqCnt, fp ) ;
    for (i = 0 ; i < _seqCnt + _extraSeqCnt ; ++i)
    {
      std::string s = _seqStrNameMap.Inverse(i) ;
      SaveString(fp, s) ;
    }
  }

  void Load(FILE *fp)
  {
    Free() ;
    InitTaxRankNum() ;

    LOAD_VAR(fp, _nodeCnt) ;
    LOAD_VAR(fp, _seqCnt) ;
    LOAD_VAR(fp, _extraSeqCnt) ;

    // Load the taxnomoy information
    _taxonomyTree = new struct TaxonomyNode[_nodeCnt] ;
    fread(_taxonomyTree, sizeof(_taxonomyTree[0]), _nodeCnt, fp) ;
    _taxIdMap.Load(fp) ;
    size_t i ;
    _taxonomyName = new std::string[_nodeCnt] ;
    for (i = 0 ; i < _nodeCnt ; ++i)
      LoadString(fp, _taxonomyName[i]) ;

    // Load the seqID information
    _seqIdToTaxId = new uint64_t[_seqCnt] ;
    fread(_seqIdToTaxId, sizeof(_seqIdToTaxId[0]), _seqCnt, fp ) ;
    for (i = 0 ; i < _seqCnt + _extraSeqCnt ; ++i)
    {
      std::string seqStr ;
      LoadString(fp, seqStr) ;
      _seqStrNameMap.Add(seqStr) ;
    }
    _rootCTaxId = FindRoot() ;
  }

  void PrintTaxonomyTree(FILE *fp)
  {
    size_t i ;
    for (i = 0 ; i < _nodeCnt ; ++i)
      printf("%lu\t|\t%lu\t|\t%s\n", 
          GetOrigTaxId(i), GetOrigTaxId( _taxonomyTree[i].parentTid ), 
          GetTaxRankString(_taxonomyTree[i].rank)) ;
  }

  void PrintNameTable(FILE *fp)
  {
    size_t i ;
    for (i = 0 ; i < _nodeCnt ; ++i)
    {
      printf("%lu\t%s\n", GetOrigTaxId(i), _taxonomyName[i].c_str()) ;
    }
  }

  void PrintConversionTable(FILE *fp)
  {
    size_t i ;
    for (i = 0 ; i < _seqCnt + _extraSeqCnt ; ++i)
      printf("%s\t%lu\n",_seqStrNameMap.Inverse(i).c_str(),
          GetOrigTaxId( SeqIdToTaxId(i) ) ) ;
  }
} ;

#endif

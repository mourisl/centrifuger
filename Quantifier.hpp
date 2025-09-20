#ifndef _MOURISL_CLASSIFIER_QUANTIFIER
#define _MOURISL_CLASSIFIER_QUANTIFIER

#include "Taxonomy.hpp"
#include "compactds/Tree_Plain.hpp"
#include "BufferManager.hpp"
#include "Classifier.hpp"

#include "defs.h"

#include <vector> 
#include <algorithm>
#include <zlib.h>

using namespace compactds ;

enum
{
  QUANTIFIER_OUTPUT_FORMAT_CENTRIFUGER,
  QUANTIFIER_OUTPUT_FORMAT_METAPHLAN,
  QUANTIFIER_OUTPUT_FORMAT_CAMI
} ;

struct _readAssignment
{
  std::vector<uint64_t> targets ;
  double weight ; // weighted number of assignment, considering the score
  double count ; // number of assignment
  double uniqCount ; // The number of assignment that it is unique

  _readAssignment() : targets() {}

  _readAssignment(const struct _readAssignment &a): targets(a.targets) 
  {
    weight = a.weight ;
    count = a.count ;
    uniqCount = a.uniqCount ;
  }

  struct _readAssignment& operator=(const struct _readAssignment &a)
  {
    targets = a.targets ;
    weight = a.weight ;
    count = a.count ;
    uniqCount = a.uniqCount ;
    return *this ;
  }

  bool operator <(const struct _readAssignment &b) const
  {
    if (targets.size() != b.targets.size())
      return targets.size() < b.targets.size() ;
    else
    {
      int i, size ;
      size = b.targets.size() ;
      for (i = 0 ; i < size ; ++i)
        if (b.targets[i] != targets[i])
          return targets[i] < b.targets[i] ;
    }
    return false ;
  }
  
  bool operator ==(const struct _readAssignment &b) const
  {
    if (targets.size() != b.targets.size())
      return false ;
    else
    {
      int i, size ;
      size = b.targets.size() ;
      for (i = 0 ; i < size ; ++i)
        if (b.targets[i] != targets[i])
          return false ;
    }
    return true ;
  }
} ;

class Quantifier
{
private:
  BufferManager<char> _buffers ;

  Taxonomy _taxonomy ;

  std::map<size_t, size_t> _seqLength ;
  size_t *_taxidLength ; // (average) genome length for each taxonomy ID.

  std::vector<struct _readAssignment> _assignments ;
  double *_abund ;
  double *_readCount ; // number of reads assigned to this tax ID and its subtree, taking the probability distribution into account.
  double *_uniqReadCount ; // number of reads uniquely assigned to this tax ID. Unique is at the strain/sequence level in its subtree.
  size_t uncountedReadCount ; // number of reads 

  bool _hasExpandedTaxIds ;
  std::map< std::pair<size_t, size_t>, double > _childReadCount ; // pair<a,b>: a: parent ctid, b: child ctid. double: count sum (each count is fractioned by the expanded tax id size)

  // NOT USED now. Original implementation of taxonomy ID genome length is taking the max, now is taking average.
  size_t GenerateTreeInternalNodeLength(size_t tag, const Tree_Plain &tree, size_t *taxidLen)
  {
    size_t i ;
    if (tree.IsLeaf(tag))
    {
      return taxidLen[tag] ; // it should be set outside
    }

    std::vector<size_t> children = tree.GetChildren(tag) ;
    size_t childrenCnt = children.size() ;
    size_t len = 0 ;
    for (i = 0 ; i < childrenCnt ; ++i)
    {
      size_t tmp = GenerateTreeInternalNodeLength(children[i], tree, taxidLen) ;
      if (tmp > len)
        len = tmp ;
    }
    return taxidLen[tag] = len ;
  }

  // Cumulative the abundance across the tree
  // The initial abund is estimated for each node separately 
  double GenerateTreeAbundance(uint64_t tag, double *abund, const Tree_Plain &tree)
  {
    size_t i ;
    std::vector<size_t> children = tree.GetChildren(tag) ;
    size_t csize = children.size() ;
    double sum = abund[tag] ;
    for (i = 0 ; i < csize ; ++i)
      sum += GenerateTreeAbundance(children[i], abund, tree) ;

    return abund[tag] = sum ;
  } 

  // Redistribute the parent node's abundance to the children.
  void RedistributeAbundToChildren(uint64_t tag, double *abund, const Tree_Plain &tree, size_t *taxidLen, double *treeEdgeWeight)
  {
    size_t i ;
    std::vector<size_t> children = tree.GetChildren(tag) ;
    size_t csize = children.size() ;
    
    double childrenSum = 0 ;
    double weightedChildrenSum = 0 ;
    for (i = 0 ; i < csize ; ++i)
    {
      childrenSum += abund[children[i]] ;
      // Since the abund represents the fraction of cells, there is no need to normalize of the genome length
      //weightedChildrenSum += abund[children[i]] / (taxidLen ? taxidLen[children[i]] : 1);
    }
    double excess = abund[tag] - childrenSum ;
    if (excess < 0)
      excess = 0 ;
    if (childrenSum == 0)
      return ;

    double expandedChildSum = 0 ;
    if (treeEdgeWeight != NULL)
    {
      for (i = 0 ; i < csize ; ++i)
        expandedChildSum += treeEdgeWeight[ children[i] ] ;
    }

    for (i = 0 ; i < csize ; ++i)
    {
      weightedChildrenSum += abund[children[i]] / (taxidLen ? taxidLen[children[i]] : 1) *
        ((excess - expandedChildSum) / csize + 
         (expandedChildSum == 0 ? 0 : treeEdgeWeight[children[i]] / expandedChildSum)) ;
    }
    
    if (weightedChildrenSum == 0) // excess == 0
      weightedChildrenSum = 1 ;

    for (i = 0 ; i < csize ; ++i)
    {
      abund[children[i]] += excess * 
        (abund[children[i]] / (taxidLen ? taxidLen[children[i]] : 1)
          * ((excess - expandedChildSum) / csize 
            + (expandedChildSum == 0 ? 0 : treeEdgeWeight[children[i]] / expandedChildSum)))
        / weightedChildrenSum ;
      RedistributeAbundToChildren(children[i], abund, tree, taxidLen, treeEdgeWeight) ;
    }
  }
  
  // Update the abund0 to abund1 using one iteration of EM
  // return: |abund1-abund0|
  double EMupdate(double *abund0, double *abund1, double *readCount, const std::vector< struct _readAssignment> &assignments, const Tree_Plain &tree, size_t *taxidLen, double *treeEdgeWeight)
  {
    size_t i, j ;
    size_t size = assignments.size() ;
    size_t treeSize = tree.GetSize() ;
    double sum ;

    memset(readCount, 0, sizeof(double) * treeSize) ;
    
    // E-step
    for (i = 0 ; i < size ; ++i)
    {
      sum = 0 ;
      const std::vector<uint64_t> &targets = assignments[i].targets ;
      size_t targetCnt = targets.size() ;

      for (j = 0 ; j < targetCnt ; ++j)
        sum += abund0[ targets[j] ] ;
      for (j = 0 ; j < targetCnt ; ++j)
      {
        readCount[ targets[j] ] += assignments[i].weight * abund0[ targets[j] ] / sum ;
      }
    }
    sum = 0 ;
    //GenerateTreeAbundance(0, readCount, tree) ;
    
    // M-step: renormalize the read count 
    double diffSum = 0 ;
    sum = 0 ;
    for (i = 0 ; i < treeSize ; ++i)
    {
      sum += readCount[i] / (double)taxidLen[i] ;
      //printf("%d: %lf %lu %lf %lf\n", i, readCount[i], taxidLen[i], readCount[i] / (double)taxidLen[i], sum) ;
    }
    for (i = 0 ; i < treeSize ; ++i)
    {
      double tmp = readCount[i] / (double)taxidLen[i] / sum ;
      abund1[i] = tmp ;
    }
    
    GenerateTreeAbundance(0, abund1, tree) ;
    // When redistrubitng using abundance, we don't need to reconsider the genome length
    RedistributeAbundToChildren(0, abund1, tree, /*taxidLen*/NULL, treeEdgeWeight);
    for (i = 0 ; i < treeSize ; ++i)
    {
      diffSum += ABS(abund0[i] - abund1[i]) ;
    }
    return diffSum ;
  }

  void EstimateAbundanceWithEM(const std::vector< struct _readAssignment > &assignments, const Tree_Plain &tree, size_t *taxidLen, double *treeEdgeWeight, double *readCount, double *abund)
  {
    size_t i, j ;

    // Initalize the abundance
    size_t assignCnt = assignments.size() ;
    double totalWeight = 0 ;
    for (i = 0 ; i < assignCnt ; ++i)
    {
      size_t targetCnt = assignments[i].targets.size() ;
      for (j = 0 ; j < targetCnt ; ++j)
        readCount[assignments[i].targets[j]] += assignments[i].weight / (double)targetCnt ;
      totalWeight += assignments[i].weight ;
    }
    
    double tmp = 0 ;
    for (i = 0 ; i < tree.GetSize() ; ++i)
      tmp += readCount[i] ;
    GenerateTreeAbundance(tree.Root(), readCount, tree) ;
    RedistributeAbundToChildren(tree.Root(), readCount, tree, taxidLen, treeEdgeWeight);
    size_t treeSize = tree.GetSize() ;
    double factor = readCount[tree.Root()] ;
    for (i = 0 ; i < treeSize ; ++i)
    {
      abund[i] = readCount[i] / factor ; 
    }
    
    // EM algorithm
    double *nextAbund = (double *)malloc(sizeof(nextAbund[0]) * treeSize) ;
    double delta = 0 ;

    const int maxIterCnt = 1000 ;
    int t ;
    for (t = 0 ; t < maxIterCnt ; ++t)
    {
      delta = EMupdate(abund, nextAbund, readCount, assignments, tree, taxidLen, treeEdgeWeight) ; 
      memcpy(abund, nextAbund, sizeof(double) * treeSize) ;
      //printf("delta: %lf\n", delta) ;
      if (delta < 1e-6 && delta < 0.1 / (double)treeSize)
        break ;
    }
    
    GenerateTreeAbundance(0, readCount, tree) ;
    RedistributeAbundToChildren(tree.Root(), readCount, tree, taxidLen, treeEdgeWeight);
    free(nextAbund) ;
  }

  double CalculateAssignmentWeight(size_t score, size_t hitLength, size_t readLength)
  {
    int diff = readLength - hitLength ;
    if (diff < int(readLength * 0.01))
      return 1 ;
    else
      diff -= int(readLength * 0.01) ;
    if (diff > 10)
      diff = 11 ;
    return 1.0 / (double)(1 << (2 * diff)) ; // Every difference decrease the probability by 1/4
  }

  // Get the taxonomy lineage string for the taxId
  // style: from quantifier_output_format
  // useName: the path string uses scientific name. (false is to use ID)
  // canonicalRankOnly: for intermediate node, only use those in the canonical rank
  // @return: length of the selected nodes. 
  int GetTaxLineagePathString(size_t ctid, int style, bool useName, bool canonicalRankOnly, std::string &pathString)
  {
    int i ;

    SimpleVector<size_t> path ;
    _taxonomy.GetTaxLineagePath(ctid, path) ;
    path.Reverse() ;

    int pathSize = path.Size() ;
    pathString = "" ;
    for (i = 0 ; i < pathSize ; ++i)
    {
      if (canonicalRankOnly && !_taxonomy.IsInCanonicalRank(path[i]))
        continue ;
      
      if (style == QUANTIFIER_OUTPUT_FORMAT_METAPHLAN
          && useName) // Add some prefix to the taxonomy name
      {
        char r[4] ; //rank descriptor
        if (_taxonomy.IsInCanonicalRank(path[i]))
        {
          if (_taxonomy.GetTaxIdRank(path[i]) == RANK_SUPER_KINGDOM || _taxonomy.GetTaxIdRank(path[i]) == RANK_ACELLULAR_ROOT)
            r[0] = 'd' ;
          else
            r[0] = _taxonomy.GetTaxRankString( _taxonomy.GetTaxIdRank(path[i]))[0] ;
          r[1] = '_' ;
          r[2] = '_' ;
          r[3] = '\0' ;
        }
        else
        {
          r[0] = r[1] = '_' ;
          r[2] = '\0' ;
        }
        pathString += r ;
      }
      
      if (useName)
        pathString += _taxonomy.GetTaxIdName(path[i]) ;
      else
      {
        char buffer[30] ;
        sprintf(buffer, "%lu", _taxonomy.GetOrigTaxId(path[i])) ;
        pathString += buffer ; 
      }
      if (i < pathSize - 1) 
        pathString += "|" ;
    }

    return pathSize ;
  }

public:
  Quantifier()
  {
    _buffers.Init(4) ;

    _buffers.Get(0, 65536) ;
    _buffers.Get(1, 65536) ;
    _buffers.Get(2, 65536) ;
    _buffers.Get(3, 65536) ;
  
    _abund = NULL ;
    _readCount = NULL ;
    _uniqReadCount = NULL ;
    _taxidLength = NULL ;
  }

  ~Quantifier()
  {
    if (_abund != NULL)
    {
      free(_taxidLength) ;
      free(_uniqReadCount) ;
      free(_readCount) ;
      free(_abund) ;
    }
  }

  // file: classificaiont output file
  // format: 0: centrifuger. Future: 1-kraken, 2-kmcp/ganon
  void Init(char *indexPrefix)
  {
    char fileName[1024] ;

    // read in the index
    _taxonomy.Free() ;
    sprintf(fileName, "%s.2.cfr", indexPrefix) ;
    FILE *fp = fopen(fileName, "r") ;
    _taxonomy.Load(fp) ;
    fclose(fp) ;

    _seqLength.clear() ;
    sprintf(fileName, "%s.3.cfr", indexPrefix) ;
    fp = fopen(fileName, "r") ;
    size_t tmp[2] ;
    while (fread(tmp, sizeof(tmp[0]), 2, fp))
      _seqLength[tmp[0]] = tmp[1] ;
    fclose(fp) ;
  
    _abund = (double *)calloc(_taxonomy.GetNodeCount() + 1, sizeof(_abund[0])) ;
    _readCount = (double *)calloc(_taxonomy.GetNodeCount() + 1, sizeof(_readCount)) ;
    _uniqReadCount = (double *)calloc(_taxonomy.GetNodeCount() + 1, sizeof(_uniqReadCount)) ;
    _taxidLength = (size_t *)calloc(_taxonomy.GetNodeCount() + 1, sizeof(size_t)) ;
    
    // Initialize genome length. It stores the average genome size if there are multiple genomes, such as internal nodes
    _taxonomy.ConvertSeqLengthToTaxLength(_seqLength, _taxidLength) ; 
  }

  void Init(char *taxonomyTree, char *nameTable, char *sizeTable)
  {
    _taxonomy.Init(taxonomyTree, nameTable) ;
    
    _abund = (double *)calloc(_taxonomy.GetNodeCount() + 1, sizeof(_abund[0])) ;
    _readCount = (double *)calloc(_taxonomy.GetNodeCount() + 1, sizeof(_readCount)) ;
    _uniqReadCount = (double *)calloc(_taxonomy.GetNodeCount() + 1, sizeof(_uniqReadCount)) ;
    _taxidLength = (size_t *)calloc(_taxonomy.GetNodeCount() + 1, sizeof(size_t)) ;

    if (sizeTable)
    {
      FILE *fp = fopen(sizeTable, "r") ;
      size_t taxid ;
      size_t length ;
      while (fscanf(fp, "%lu %lu", &taxid, &length) != EOF)
      {
        _taxidLength[ _taxonomy.CompactTaxId(taxid) ] = length ;
      }
      _taxonomy.InferAllTaxLength(_taxidLength, false) ;
      fclose(fp) ;
    }
    else
    {
      size_t i ;
      for (i = 0 ; i < _taxonomy.GetNodeCount() ; ++i)
        _taxidLength[i] = 1000000 ;
    }
  }
  
  // Coalsce the assignment that mapped to the same set of target
  size_t CoalesceAssignments()
  {
    size_t i, k ;
    size_t size = _assignments.size() ;
    std::sort(_assignments.begin(), _assignments.end()) ;
    
    k = 1 ;
    for (i = 1 ; i < size ; ++i)
    {
      if (_assignments[i] == _assignments[k - 1])
      {
        _assignments[k - 1].weight += _assignments[i].weight ;
        _assignments[k - 1].count += _assignments[i].count ;
        _assignments[k - 1].uniqCount += _assignments[i].uniqCount ;
      }
      else
      {
        _assignments[k] = _assignments[i] ;
        ++k ;
      }
    }
    _assignments.resize(k) ;
    return k ;
  }
  
  void LoadReadAssignments(char *file, uint64_t minScore, uint64_t minHitLength, int format)
  {
    _assignments.clear() ;
    
    // read in the classificaiton result
    size_t lineCnt = 0 ;
    gzFile gzfp = strcmp(file, "-") ? gzopen(file, "r") : gzdopen(fileno(stdin), "r");

    char *line =  _buffers.Get(0, 0) ;
    char *readId = _buffers.Get(2, 0) ;
    char *prevReadId = _buffers.Get(3, 0) ;

    struct _readAssignment assign ;
    prevReadId[0] = '\0' ;
    bool hasExpandedTaxIds = false ;
    std::vector<uint64_t> expandedTaxIds ;

    while (gzgets(gzfp, line, sizeof(char) * _buffers.GetBufferSize(0)))
    {
      if (lineCnt == 0) // header
      {
        //if (strstr(line, "expandedTaxIDs"))
        //  hasExpandedTaxIds = true ;
        ++lineCnt ; 
        continue ;
      }

      char *buffer = _buffers.Get(1, 0) ;
      uint64_t taxid, score, secondScore, hitLength, readLength ;
      sscanf(line, "%s\t%[^\t]\t%lu\t%lu\t%lu\t%lu\t%lu", readId, buffer, &taxid, &score, &secondScore, &hitLength, &readLength) ;
      if (hitLength < minHitLength || score < minScore || taxid == 0)
        continue ;

      if (hasExpandedTaxIds) // It is not set. The feature is disabled for now.
      {
        //TODO: find the actual column. Currently assume the last column is the expanded taxIdx
        int lineLength = strlen(line) ;
        int i = lineLength - 1 ;
        int j, k = 0 ;
        if (line[i] == '\n')
          --i ;
        for (; i >= 0 && line[i] != '\t' ; --i)
        {
          buffer[k] = line[i] ;
          ++k ;
        }
        k = '\0' ;

        for (i = 0, j = k - 1 ; i < j ; ++i, --j)
        {
          char tmp = buffer[i] ;
          buffer[i] = buffer[j] ;
          buffer[j] = tmp ;
        }

        uint64_t childTid = 0 ;
        expandedTaxIds.clear() ;
        for (i = 0 ; i < k ; ++i)
        {
          if (buffer[i] >= '0' && buffer[i] <= '9')
          {
            childTid = childTid * 10 + buffer[i] - '0' ; 
          }
          else
          {
            expandedTaxIds.push_back(childTid) ;
            childTid = 0 ;
          }
        }

        int expandedTaxIdSize = expandedTaxIds.size() ;
        for (i = 0 ; i < expandedTaxIdSize ; ++i)
        {
          std::pair<size_t, size_t> p ;
          p.first = _taxonomy.CompactTaxId(taxid) ;
          p.second = _taxonomy.CompactTaxId(expandedTaxIds[i]) ;
          _childReadCount[p] += CalculateAssignmentWeight(score, hitLength, readLength) / (double)expandedTaxIdSize ;
        }
      }

      if (strcmp(readId, prevReadId))
      {
        if (prevReadId[0] != '\0' && assign.targets.size() > 0)
          _assignments.push_back(assign) ;

        assign.targets.clear() ;
        assign.weight = CalculateAssignmentWeight(score, hitLength, readLength) ;
        assign.count = 1 ;
        assign.uniqCount = score > secondScore ? 1 : 0 ;

        strcpy(prevReadId, readId) ;
      }
      assign.targets.push_back(_taxonomy.CompactTaxId(taxid)) ;
      ++lineCnt ;
    
      // Reduce the size about every 10,000,000 assignments 
      if (lineCnt % 10000000 == 0)
        CoalesceAssignments() ;
    }
    if (assign.targets.size() > 0)
      _assignments.push_back(assign) ;
    gzclose(gzfp) ;
    CoalesceAssignments() ;
  }

  void AddReadAssignment(const struct _classifierResult &result)
  {
    int i ;
    struct _readAssignment assign ;
    int size = result.taxIds.size() ;
    for (i = 0 ; i < size ; ++i)
      assign.targets.push_back( _taxonomy.CompactTaxId(result.taxIds[i])) ; 
    assign.weight = CalculateAssignmentWeight(result.score, result.hitLength, 
        result.queryLength) ;
    assign.count = 1 ;
    assign.uniqCount = result.score > result.secondaryScore ? 1 : 0 ;

    _assignments.push_back(assign) ;
  }

  // Main function. Should be called after Init and set up the read assignment
  void Quantification()
  {
    size_t i, j ;

    CoalesceAssignments() ;
    size_t assignCnt = _assignments.size() ;
    std::vector< struct _readAssignment> subtreeAssignments ; // the assignment where the targe ID is with respect to the subtree 

    Tree_Plain allTree ;
    _taxonomy.ConvertToGeneralTree(allTree) ;

    // Reduce the _taxonomy tree to the nodes with some coverage. Also we need the tree to have pointers to the children
    MapID<size_t> coveredTaxIds ; 
    Tree_Plain subtree ; // the subtree that has some coverage information

    // Get the subtree's ids and also convert read assign to the subtree node
    size_t subtreeSize = 1 ;
    coveredTaxIds.Add( allTree.Root() ) ; // Map the root to 0.
    for (i = 0 ; i < assignCnt ; ++i)  
    {
      size_t targetCnt = _assignments[i].targets.size() ;
      subtreeAssignments.push_back( _assignments[i] ) ;
      
      for (j = 0 ; j < targetCnt ; ++j)
      {
        uint64_t ctid = subtreeAssignments[i].targets[j] ;
        // If a read hit a node not in the tree, we set it to the root
        if (ctid == _taxonomy.GetNodeCount())
        {
          subtreeAssignments[i].targets[j] = 0 ; // Subtree's root is 0. We allow duplicated 0 here, so the probability reflect the underlying read assignment
          _readCount[ allTree.Root() ] += _assignments[i].count / targetCnt;
          _uniqReadCount[ allTree.Root() ] += _assignments[i].uniqCount ; // targetCnt must be 1, otherwise uniqCount will be 0.
          continue ;
        }
        _readCount[ _assignments[i].targets[j] ] += _assignments[i].count / targetCnt;
        _uniqReadCount[ _assignments[i].targets[j] ] += _assignments[i].uniqCount ; // targetCnt must be 1, otherwise uniqCount will be 0.

        uint64_t p = ctid ;
        while (coveredTaxIds.Add(p) == subtreeSize)
        {
          ++subtreeSize ;
          p = _taxonomy.GetParentTid(p) ;
        }
        subtreeAssignments[i].targets[j] = coveredTaxIds.Map(ctid) ;
      }
      subtreeAssignments[i].targets.resize(targetCnt) ;
    }
    GenerateTreeAbundance(allTree.Root(), _readCount, allTree) ;
    GenerateTreeAbundance(allTree.Root(), _uniqReadCount, allTree) ;
    
    // Create the subtree's structure
    subtree.Init(subtreeSize) ;
    // 0 is the root based on the mapping, so we start from i=1
    for (i = 1 ; i < subtreeSize ; ++i)
      subtree.AddEdge(i, coveredTaxIds.Map(_taxonomy.GetParentTid( coveredTaxIds.Inverse(i)))) ;
        
    // Copy the tax Id length to subtree 
    size_t *subtreeTaxIdLen = (size_t *)calloc(subtreeSize, sizeof(*subtreeTaxIdLen)) ;
    size_t allNodeCnt = allTree.GetSize() ;
    for (i = 0 ; i < allNodeCnt ; ++i) 
    {
      if (coveredTaxIds.IsIn(i))
      {
        subtreeTaxIdLen[coveredTaxIds.Map(i)] = _taxidLength[i] + _taxidLength[_taxonomy.GetRoot()] / 10 ; // Adding a baseline length to avoid extremely short genome confounding the length
      }
    }
    
    // Initialize abundance
    double *subtreeAbund = (double *)calloc(subtreeSize, sizeof(*subtreeAbund)) ; 
    double *subtreeReadCount = (double *)calloc(subtreeSize, sizeof(*subtreeReadCount)) ; 
    double *subtreeEdgeWeight = (double *)calloc(subtreeSize, sizeof(*subtreeEdgeWeight)) ; // The edge from node i to its parent. The weight is the number of normalized read count 

    if (_hasExpandedTaxIds)
    {
      for (std::map< std::pair<size_t, size_t>, double>::iterator iter = _childReadCount.begin() ;
          iter != _childReadCount.end() ; ++iter)
      {
        std::pair<size_t, size_t> p = iter->first ;
        double weight = iter->second ;
        if (coveredTaxIds.IsIn(p.first) && coveredTaxIds.IsIn(p.second))
        {
          size_t mf, ms ; //mapped p.first, mapped p.second
          mf = coveredTaxIds.Map(p.first) ;
          ms = coveredTaxIds.Map(p.second) ;
          if (subtree.Parent(ms) == mf)
            subtreeEdgeWeight[ms] = weight ;
        }
      }
    }
    
    // Start the calculation using EM.
    EstimateAbundanceWithEM(subtreeAssignments, subtree, subtreeTaxIdLen, 
        _hasExpandedTaxIds ? subtreeEdgeWeight : NULL, 
        subtreeReadCount, subtreeAbund) ;

    for (i = 0 ; i < subtreeSize; ++i)
    {
      _abund[ coveredTaxIds.Inverse(i) ] = subtreeAbund[i] ;
    }
    free(subtreeAbund) ;
    free(subtreeReadCount) ;
    free(subtreeTaxIdLen) ;
    free(subtreeEdgeWeight) ;
  }

  // format: check the enum 
  void Output(FILE *fp, int format)
  {
    size_t i ;
    size_t nodeCnt = _taxonomy.GetNodeCount() ;
      
    if (format == QUANTIFIER_OUTPUT_FORMAT_METAPHLAN)
    {
      fprintf(fp, "#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species\n") ;      
      for (i = 0 ; i < nodeCnt ; ++i)
      {
        if (_readCount[i] < 1e-6)
          continue ;
        if (!_taxonomy.IsInCanonicalRank(i))
          continue ;
        //if (_abund[i] < 5e-8)
        //  continue ;

        std::string taxIdPathString ;
        std::string taxNamePathString ;

        GetTaxLineagePathString(i, format, false, true, taxIdPathString) ;
        GetTaxLineagePathString(i, format, true, true, taxNamePathString) ;

        fprintf(fp, "%s\t%s\t%.5lf\t\n", taxNamePathString.c_str(), taxIdPathString.c_str(), _abund[i] * 100.0 ) ; 
      }
    }
    else if (format == QUANTIFIER_OUTPUT_FORMAT_CAMI)
    {
      fprintf(fp, "@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n") ;
      for (i = 0 ; i < nodeCnt ; ++i)
      {
        if (_readCount[i] < 1e-6)
          continue ;
        if (!_taxonomy.IsInCanonicalRank(i))
          continue ;
        //if (_abund[i] < 5e-8)
        //  continue ;
        
        std::string taxIdPathString ;
        std::string taxNamePathString ;

        GetTaxLineagePathString(i, format, false, true, taxIdPathString) ;
        GetTaxLineagePathString(i, format, true, true, taxNamePathString) ;
        
        fprintf(fp, "%lu\t%s\t%s\t%s\t%.5lf\n", _taxonomy.GetOrigTaxId(i), _taxonomy.GetTaxRankString(_taxonomy.GetTaxIdRank(i)), taxIdPathString.c_str(), taxNamePathString.c_str(), _abund[i] * 100.0 ) ; 
      }
    }
    else
    {
      if (format != QUANTIFIER_OUTPUT_FORMAT_CENTRIFUGER)
        Utils::PrintLog("Warning: unknown output format, will output in Centrifuger format.\n") ;
      
      fprintf(fp, "name\ttaxID\ttaxRank\tgenomeSize\tnumReads\tnumUniqueReads\tabundance\n") ;
      for (i = 0 ; i < nodeCnt ; ++i)
      {
        if (_readCount[i] < 1e-6)
          continue ;

        fprintf(fp, "%s\t%lu\t%s\t%lu\t%d\t%d\t%.7lf\n",
            _taxonomy.GetTaxIdName(i).c_str(), 
            _taxonomy.GetOrigTaxId(i),
            _taxonomy.GetTaxRankString( _taxonomy.GetTaxIdRank(i)),
            _taxidLength[i], 
            (int)(_readCount[i] + 1e-3), (int)(_uniqReadCount[i] + 1e-3), _abund[i]) ;
      }
    }
  }
} ;

#endif

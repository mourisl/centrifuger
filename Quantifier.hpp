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

struct _readAssignment
{
  std::vector<uint64_t> targets ;
  double weight ;
  double uniqWeight ; // The number of assignment that it is unique
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
  double *_readCount ; // number of reads assigned to this tax ID and its subtree
  double *_uniqReadCount ; // number of reads uniquely assigned to this tax ID. Unique is at the strain/sequence level in its subtree.

  // NOT USED now. Original implementation of taxonomy ID genome length is taking the max, now is taking average.
  size_t GenerateTreeInternalNodeLength(size_t tag, const Tree_Plain &tree, size_t *taxidLen)
  {
    size_t i ;
    if (tree.IsLeaf(tag))
    {
      return taxidLen[tag] ; // it is set through the all-tree
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
  void RedistributeAbundToChildren(uint64_t tag, double *abund, const Tree_Plain &tree, size_t *taxidLen)
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
      weightedChildrenSum += abund[children[i]] ; //* taxidLen[children[i]] ;
    }
    double excess = abund[tag] - childrenSum ;
    if (excess < 0)
      excess = 0 ;
    for (i = 0 ; i < csize ; ++i)
    {
      abund[children[i]] += excess * (abund[children[i]]) / weightedChildrenSum ;
      RedistributeAbundToChildren(children[i], abund, tree, taxidLen) ;
    }
  }
  
  // Update the abund0 to abund1 using one iteration of EM
  // return: |abund1-abund0|
  double EMupdate(double *abund0, double *abund1, double *readCount, const std::vector< struct _readAssignment> &assignments, const Tree_Plain &tree, size_t *taxidLen)
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
      sum += readCount[i] / (double)taxidLen[i] ;
    for (i = 0 ; i < treeSize ; ++i)
    {
      double tmp = readCount[i] / (double)taxidLen[i] / sum ;
      abund1[i] = tmp ;
    }
    
    GenerateTreeAbundance(0, abund1, tree) ;
    RedistributeAbundToChildren(0, abund1, tree, taxidLen);
    for (i = 0 ; i < treeSize ; ++i)
    {
      diffSum += ABS(abund0[i] - abund1[i]) ;
    }
    return diffSum ;
  }


  void EstimateAbundanceWithEM(const std::vector< struct _readAssignment > &assignments, const Tree_Plain &tree, size_t *taxidLen, double *abund)
  {
    size_t i, j ;

    // Initalize the abundance
    size_t assignCnt = assignments.size() ;
    for (i = 0 ; i < assignCnt ; ++i)
    {
      size_t targetCnt = assignments[i].targets.size() ;
      for (j = 0 ; j < targetCnt ; ++j)
        abund[assignments[i].targets[j]] = 1.0 / (double)targetCnt ;
    }
    GenerateTreeAbundance(tree.Root(), abund, tree) ;
    RedistributeAbundToChildren(tree.Root(), abund, tree, taxidLen);
    
    // EM algorithm
    size_t treeSize = tree.GetSize() ;
    double *nextAbund = (double *)malloc(sizeof(nextAbund[0]) * treeSize) ;
    double *readCount = (double *)malloc(sizeof(readCount[0]) * treeSize) ;
    double delta = 0 ;

    const int maxIterCnt = 1000 ;
    int t ;
    for (t = 0 ; t < maxIterCnt ; ++t)
    {
      delta = EMupdate(abund, nextAbund, readCount, assignments, tree, taxidLen) ; 
      memcpy(abund, nextAbund, sizeof(double) * treeSize) ;
      //printf("delta: %lf\n", delta) ;
      if (delta < 1e-6 && delta < 0.1 / (double)treeSize)
        break ;
    }
    free(nextAbund) ;
    free(readCount) ;
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
  
    _abund = (double *)calloc(_taxonomy.GetNodeCount(), sizeof(_abund[0])) ;
    _readCount = (double *)calloc(_taxonomy.GetNodeCount(), sizeof(_readCount)) ;
    _uniqReadCount = (double *)calloc(_taxonomy.GetNodeCount(), sizeof(_uniqReadCount)) ;
    _taxidLength = (size_t *)calloc(_taxonomy.GetNodeCount(), sizeof(size_t)) ;
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
        _assignments[k - 1].uniqWeight += _assignments[i].uniqWeight ;
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
  
  void LoadReadAssignments(char *file, int format)
  {
    _assignments.clear() ;
    
    // read in the classificaiton result
    size_t lineCnt = 0 ;
    gzFile gzfp = gzopen(file, "r") ;

    char *line =  _buffers.Get(0, 0) ;
    char *readId = _buffers.Get(2, 0) ;
    char *prevReadId = _buffers.Get(3, 0) ;

    struct _readAssignment assign ;
    prevReadId[0] = '\0' ;
    while (gzgets(gzfp, line, sizeof(char) * _buffers.GetBufferSize(0)))
    {
      if (lineCnt == 0) // header
      {
        ++lineCnt ; 
        continue ;
      }

      char *buffer = _buffers.Get(1, 0) ;
      uint64_t taxid, score, secondScore ;
      sscanf(line, "%s\t%[^\t]\t%lu\t%lu\t%lu", readId, buffer, &taxid, &score, &secondScore) ;
      if (strcmp(readId, prevReadId))
      {
        if (prevReadId[0] != '\0')
          _assignments.push_back(assign) ;

        assign.targets.clear() ;
        assign.weight = 1 ;
        assign.uniqWeight = score > secondScore ? 1 : 0 ;

        strcpy(prevReadId, readId) ;
      }
      assign.targets.push_back(_taxonomy.CompactTaxId(taxid)) ;
      ++lineCnt ;
    
      // Reduce the size about every 1,000,000 assignments 
      if (lineCnt % 1000000 == 0)
        CoalesceAssignments() ;
    }
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
    assign.weight = 1 ;
    assign.uniqWeight = result.score > result.secondaryScore ? 1 : 0 ;

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
      size_t nohitCnt = 0 ;
      subtreeAssignments.push_back( _assignments[i] ) ;
      
      for (j = 0 ; j < targetCnt ; ++j)
      {
        if (subtreeAssignments[i].targets[j] == 0)
        {
          ++nohitCnt ;
          continue ;
        }
        
        uint64_t ctid = subtreeAssignments[i].targets[j] ; 
        // If a read hit a node not in the tree, we set it to the root
        if (ctid == _taxonomy.GetNodeCount())
        {
          subtreeAssignments[i].targets[j] = 0 ; // Subtree's root is 0.
          _readCount[ allTree.Root() ] += _assignments[i].weight / targetCnt ;
          _uniqReadCount[ allTree.Root() ] += _assignments[i].uniqWeight ; // targetCnt must be 1, otherwise uniqWeight will be 0.
          continue ;
        }
        _readCount[ _assignments[i].targets[j] ] += _assignments[i].weight / targetCnt ;
        _uniqReadCount[ _assignments[i].targets[j] ] += _assignments[i].uniqWeight ; // targetCnt must be 1, otherwise uniqWeight will be 0.


        uint64_t p = ctid ;
        while (coveredTaxIds.Add(p) == subtreeSize)
        {
          ++subtreeSize ;
          p = _taxonomy.GetParentTid(p) ;
        }
        subtreeAssignments[i].targets[j] = coveredTaxIds.Map(ctid) ;
      }
      subtreeAssignments[i].targets.resize(targetCnt - nohitCnt) ;
    }
    GenerateTreeAbundance(allTree.Root(), _readCount, allTree) ;
    GenerateTreeAbundance(allTree.Root(), _uniqReadCount, allTree) ;
    
    // Create the subtree's structure
    subtree.Init(subtreeSize) ;
    // 0 is the root based on the mapping, so we start from i=1
    for (i = 1 ; i < subtreeSize ; ++i)
      subtree.AddEdge(i, coveredTaxIds.Map(_taxonomy.GetParentTid( coveredTaxIds.Inverse(i)))) ;
        
    // Initialize genome length. It stores the average genome size if there are multiple genomes, such as internal nodes
    _taxonomy.ConvertSeqLengthToTaxLength(_seqLength, _taxidLength) ; 
   
    // Copy the tax Id length to subtree 
    size_t *subtreeTaxIdLen = (size_t *)calloc(subtreeSize, sizeof(*subtreeTaxIdLen)) ;
    size_t allNodeCnt = allTree.GetSize() ;
    for (i = 0 ; i < allNodeCnt ; ++i) 
    {
      if (coveredTaxIds.IsIn(i))
      {
        subtreeTaxIdLen[coveredTaxIds.Map(i)] = _taxidLength[i] ;
      }
    }
    
    // Initialize abundance
    double *subtreeAbund = (double *)calloc(subtreeSize, sizeof(*subtreeAbund)) ; 
    
    // Start the calculation using EM.
    EstimateAbundanceWithEM(subtreeAssignments, subtree, subtreeTaxIdLen, subtreeAbund) ;

    for (i = 0 ; i < subtreeSize; ++i)
      _abund[ coveredTaxIds.Inverse(i) ] = subtreeAbund[i] ;
    
    free(subtreeAbund) ;
    free(subtreeTaxIdLen) ;
  }

  // format: 0-centrifuge's report
  //        1 - kraken
  void Output(FILE *fp, int format)
  {
    size_t i ;

    // Get the read assignment information
    fprintf(fp, "name\ttaxID\ttaxRank\tgenomeSize\tnumReads\tnumUniqueReads\tabundance\n") ;
    size_t nodeCnt = _taxonomy.GetNodeCount() ;
    for (i = 0 ; i < nodeCnt ; ++i)
    {
      if (_readCount[i] < 1e-6)
        continue ;

      printf("%s\t%lu\t%s\t%lu\t%d\t%d\t%lf\n",
          _taxonomy.GetTaxIdName(i).c_str(), 
          _taxonomy.GetOrigTaxId(i),
          _taxonomy.GetTaxRankString( _taxonomy.GetTaxIdRank(i)),
          _taxidLength[i], 
          (int)(_readCount[i] + 1e-3), (int)(_uniqReadCount[i] + 1e-3), _abund[i]) ;
    }
  }
} ;

#endif

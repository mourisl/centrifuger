#!/usr/bin/env python3

# The script to extract information in taxonomy structures
import sys
import argparse

def ReadTaxonomyTree(tax_file):
  taxonomyTree = {}
  fp = open(tax_file)
  for line in fp:
    fields = line.strip().split('\t')
    #assert len(fields) == 5
    tax_id, parent_tax_id, rank = fields[0], fields[2], fields[4]
    assert tax_id not in taxonomyTree
    taxonomyTree[tax_id] = [parent_tax_id, rank]       
  fp.close()
  return taxonomyTree

# @return: taxonomy ids (represented by a set) in the subtree using taxid as the root
def GetSubTree(taxonomyTree, taxid):
  isInSubTree = {taxid:True}
  ret = set({taxid})
  # For each taxonomy id, search the path to the root or a node 
  #   with known subtree information.
  for tid in taxonomyTree:
    if (tid in isInSubTree):
      continue

    flag = False
    path = []
    while (True):
      path.append(tid)
      parentTid = taxonomyTree[tid][0]
      if (parentTid in isInSubTree):
        flag = isInSubTree[parentTid]
        break 

      if (tid == parentTid):
        break
      else:
        tid = taxonomyTree[tid][0]
    
    for tid in path:
      isInSubTree[tid] = flag
      if (flag):
        ret.add(tid)

  return ret

def PromoteTaxLevel(taxonomyTree, taxid, rank):
  tid = taxid
  if (tid not in taxonomyTree):
    return -1

  while (True):
    if (taxonomyTree[tid][1] == rank):
      return tid 
    else:
      parentTid = taxonomyTree[tid][0]
      if (parentTid == tid):
        break
      tid = parentTid ;
  return -1 

if (__name__ == "__main__"):
  parser = argparse.ArgumentParser(description="")
  parser.add_argument("--op", help="operation: subtree,promote", dest="op", required=True)
  parser.add_argument("--tree", help="tree structure, usuallay nodes.dmp", dest="tree", required=True)
  parser.add_argument("--taxid", help="taxonomy id for operation", dest="taxid", required=False)
  parser.add_argument("--taxid-list", help="taxonomy id list", dest="taxidListFile", required=False)
  parser.add_argument("--rank", help="taxonomy rank (species, rank,..)", dest="taxRank", required=False)

  args = parser.parse_args()
  
  taxonomyTree = ReadTaxonomyTree(args.tree)

  if (args.op == "subtree"):
    for taxid in GetSubTree(taxonomyTree, args.taxid):
      print(taxid)
  elif (args.op == "promote"):
    taxidList = []
    fp = open(args.taxidListFile, "r")
    for line in fp:
      taxidList.append(line.rstrip())
    fp.close()
    
    for taxid in taxidList:
      print(PromoteTaxLevel(taxonomyTree, taxid, args.taxRank))
    

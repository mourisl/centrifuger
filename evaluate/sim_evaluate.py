#!/usr/bin/env python3

# Evaluate classfication result at various level
#   Adapted from: https://github.com/DaehwanKimLab/centrifuge/blob/master/evaluation/test/centrifuge_evaluate_mason.py
# To make it more general, it assumes the truth table, prediction table are two-column tsv files, where the first column is the read names/ids, and the second column is the tax ids. Another required input is the taxonomoy tree (nodes.dmp).

import sys
import argparse

"""
"""
def read_taxonomy_tree(tax_file):
  taxonomy_tree = {}
  fp = open(tax_file)
  for line in fp:
    fields = line.strip().split('\t')
    #assert len(fields) == 5
    tax_id, parent_tax_id, rank = fields[0], fields[2], fields[4]
    assert tax_id not in taxonomy_tree
    taxonomy_tree[tax_id] = [parent_tax_id, rank]       
  fp.close()
  return taxonomy_tree
"""
"""

"""
"""
def compare_scm(class_out, true_out, taxonomy_tree, rank):
  higher_ranked = {}
        
  ancestors = set() # non-leaf nodes
  for tax_id in taxonomy_tree.keys():
    if tax_id in ancestors:
      continue
    while True:
      parent_tax_id, cur_rank = taxonomy_tree[tax_id]
      if parent_tax_id in ancestors:
        break
      if tax_id == parent_tax_id:
        break
      tax_id = parent_tax_id
      ancestors.add(tax_id)

  db_dic = {}
  first = True
  
  for read_name, tax_id in class_out:
    # Traverse up taxonomy tree to match the given rank parameter
    rank_tax_id = tax_id
    if rank != "strain":
      while True:
        if tax_id not in taxonomy_tree:
          rank_tax_id = ""
          break
        parent_tax_id, cur_rank = taxonomy_tree[tax_id]
        if cur_rank == rank:
          rank_tax_id = tax_id
          break
        if tax_id == parent_tax_id:
          rank_tax_id = ""
          break
        tax_id = parent_tax_id
    else:
      assert rank == "strain"
      if tax_id in ancestors:
        continue

    if rank_tax_id == "":
      # higher_ranked[read_name] = True  
        continue
    
    if read_name not in db_dic:
      db_dic[read_name] = set()
    db_dic[read_name].add(rank_tax_id)

  classified, unclassified, unique_classified = 0, 0, 0
  for read_name, tax_id in true_out: 
    #read_name, tax_id = fields[0], fields[1]
    # Traverse up taxonomy tree to match the given rank parameter
    rank_tax_id = tax_id
    if (rank != "strain"):
      while True:
        if tax_id not in taxonomy_tree:
          rank_tax_id = ""
          break
        parent_tax_id, cur_rank = taxonomy_tree[tax_id]
        if cur_rank == rank:
          rank_tax_id = tax_id
          break
        if tax_id == parent_tax_id:
          rank_tax_id = ""
          break
        
        tax_id = parent_tax_id
    if rank_tax_id == "":
      continue
    if read_name not in db_dic:
      unclassified += 1
      continue
    maps = db_dic[read_name]
    if rank_tax_id in maps:
      classified += 1
      if len(maps) == 1 and read_name not in higher_ranked:
        unique_classified += 1
    else:
      unclassified += 1
      # daehwan - for debugging purposes
      # print read_name

  raw_unique_classified = 0
  for read_name, maps in db_dic.items():
    if len(maps) == 1 and read_name not in higher_ranked:
      raw_unique_classified += 1
  return classified, unique_classified, unclassified, len(db_dic), raw_unique_classified

"""
"""

def evaluate(predictFile, truthFile, taxonomyFile, ranks, options):
  taxonomyTree = read_taxonomy_tree(taxonomyFile)

  truth = []
  fp = open(truthFile)
  for line in fp:
    cols = line.rstrip().split()
    truth.append((cols[0], cols[1]))
  fp.close()  
  
  predict = []
  fp = open(predictFile)
  for line in fp:
    if ("seqID" in line):
      continue
    read_name, tax_id, score = line.strip().split('\t')
    predict.append( (read_name, tax_id) )
  fp.close() 

  results = {"strain"  : [0, 0, 0],
               "species" : [0, 0, 0],
               "genus"   : [0, 0, 0],
               "family"  : [0, 0, 0],
               "order"   : [0, 0, 0],
               "class"   : [0, 0, 0],
               "phylum"  : [0, 0, 0]}
  for rank in ranks:
    classified, unique_classified, unclassified, raw_classified, raw_unique_classified = \
        compare_scm(predict, truth, taxonomyTree, rank)
    results[rank] = [classified, unique_classified, unclassified]
    num_cases = classified + unclassified
    # if rank == "strain":
    #    assert num_cases == num_fragment
    print("\t\t%s" % rank)
    print("\t\t\tsensitivity: {:,} / {:,} ({:.2%})".format(classified, num_cases, float(classified) / num_cases))
    print("\t\t\tprecision  : {:,} / {:,} ({:.2%})".format(classified, raw_classified, float(classified) / raw_classified))
    print("\n\t\t\tfor uniquely classified ")
    print("\t\t\t\t\tsensitivity: {:,} / {:,} ({:.2%})".format(unique_classified, num_cases, float(unique_classified) / num_cases))
    print("\t\t\t\t\tprecision  : {:,} / {:,} ({:.2%})".format(unique_classified, raw_unique_classified, float(unique_classified) / raw_unique_classified))

if (__name__ == "__main__"):
  rank_list_default = "strain,species,genus,family,order,class,phylum"
  
  parser = argparse.ArgumentParser(description="Evaluate the classification results given the truth tabel")
  parser.add_argument("--truth", help="truth data (read_id strain_tax_id)", dest="truth", required=True)
  parser.add_argument("-c", help="classification result", dest="classificationRes", required=True)
  parser.add_argument("--tree", help="taxonomy tree", dest="tree", required=True)
  parser.add_argument("--tool", help="classification method", dest="tool", type=str,
      default="centrifuge")
  #parser.add_argument("--level", help="evaluate at the specified taxonomy level", dest="level", default="strain") 
  parser.add_argument("--rank-list",
        dest="ranks",
        type=str,
        default=rank_list_default,
        help="A comma-separated list of ranks (default: %s)" % rank_list_default)
  args = parser.parse_args() 
  options = {}
  evaluate(args.classificationRes, args.truth, args.tree, args.ranks.split(','), options)
  
  

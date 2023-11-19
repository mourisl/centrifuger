#ifndef _MOURISL_COMPACTDS_TREE_CARDINAL
#define _MOURISL_COMPACTDS_TREE_CARDINAL

#include "Utils.hpp"

namespace compactds {
  
class Tree_Cardinal: public Tree 
{
protected:
  size_t _c ; // cardinality (number of max children)
public:
  Tree_Cardinal() 
  {
    _space = 0 ;
    _n = 0 ;
    _c = 0 ;
  }
  ~Tree_Cardinal() {}

  // Number of children with label l. 1: has such children. 0-don't 
  virtual size_t ChildrenLabeled(size_t v, size_t l) const = 0 ;
  // The child with label l.
  virtual size_t LabeledChild(size_t v, size_t l) const = 0 ;
  // The label of the edge that leads to node v.
  virtual size_t ChildLabel(size_t v) const = 0 ;

  virtual void Save(FILE *fp)
  {
    Tree::Save(fp) ;
    SAVE_VAR(fp, _c) ;
  }

  virtual void Load(FILE *fp)
  {
    Tree::Load(fp) ;
    LOAD_VAR(fp, _c) ;
  }
} ;

} // end of namespace

#endif

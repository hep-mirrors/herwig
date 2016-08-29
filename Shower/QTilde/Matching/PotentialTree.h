// -*- C++ -*-
#ifndef HERWIG_PotentialTree_H
#define HERWIG_PotentialTree_H
//
// This is the declaration of the PotentialTree class.
//

#include "CKKWTree.h"

namespace Herwig {

using namespace ThePEG;

/**
 *  Struct to store a potential CKKWTree 
 */  
struct PotentialTree {

  /**
   *  Constructor
   */
  PotentialTree() {}

  /**
   *  Constructor
   */
  PotentialTree(CKKWTreePtr itree, DiagPtr idiag, 
		Ptr<ColourLines>::transient_const_pointer icl) 
    : tree_(itree), diagram_(idiag), cl_(icl), wgt_(0.)
  {}

  /**
   *  Get the tree
   */
  CKKWTreePtr tree() const {return tree_;}

  /**
   *  set the tree
   */
  void tree(CKKWTreePtr in) {tree_=in;}

  /**
   *  Get the diagram
   */
  tcDiagPtr diagram() {return diagram_;}

  /**
   *  set the tree
   */
  void diagram(tcDiagPtr in) {diagram_=in;}

  /**
   *  All diagrams
   */
  MEBase::DiagramVector & diagrams() {return diagrams_;}

  /**
   *  Colour Structure
   */
  Ptr<ColourLines>::transient_const_pointer colourLines() {return cl_;}

  /**
   *  Colour Structure
   */
  void colourLines(Ptr<ColourLines>::transient_const_pointer in) {cl_=in;}

  /**
   *  Set the weight
   */
  void weight(double weight) {wgt_ = weight;}

  /**
   *  Get the weight
   */
  double weight() const {return wgt_;}

private:

  /**
   *  The tree
   */
  CKKWTreePtr tree_;

  /**
   *  All diagrams
   */
  MEBase::DiagramVector diagrams_;

  /**
   *  The diagram
   */
  tcDiagPtr diagram_;

  /**
   *  The colour structure
   */
  Ptr<ColourLines>::transient_const_pointer cl_;

  /**
   *  The weight
   */
  double wgt_;
};

}

#endif /* HERWIG_PotentialTree_H */

// -*- C++ -*-
//
// TwoBodyDecay.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_TwoBodyDecay_H
#define HERWIG_TwoBodyDecay_H
//
// This is the declaration of the TwoBodyDecay struct.
//

namespace Herwig {
using namespace ThePEG;
using Helicity::tVertexBasePtr;

/**
 * Struct for the prototype of a two-body decay mode
 */
struct TwoBodyDecay {
  
public:
  
  /**
   *  Constructor
   * @param pa Decaying particle
   * @param pb First  decay product
   * @param pc Second decay product
   */    
  TwoBodyDecay(tPDPtr pa, tPDPtr pb, tPDPtr pc,
	       tVertexBasePtr vertex) : parent_(pa), vertex_(vertex) {
    ParticleOrdering order;
    if( order(pb, pc) ) {
      children_.first = pb;
      children_.second = pc;
    }
    else {
      children_.first = pc;
      children_.second = pb;
    }
  }
  
  /**
   *  The parent
   */    
  tPDPtr parent_;
  
  /**
   *  The children
   */
  tPDPair children_;
  
  /**
   *  Vertex
   */
  tVertexBasePtr vertex_;
  
private:
  
  TwoBodyDecay();
};

/**
 * Test whether one TwoBodyDecay is less than another
 */
inline bool operator<(const TwoBodyDecay & x, const TwoBodyDecay & y) {
  if(x.parent_->id()!=y.parent_->id())
    return x.parent_->id()<y.parent_->id();
  if(x.children_.first->id()!=y.children_.first->id())
    return x.children_.first->id() < y.children_.first->id();
  if(x.children_.second->id()!=y.children_.second->id())
    return x.children_.second->id() < y.children_.second->id();
  return x.vertex_<y.vertex_;
}

}
#endif /* HERWIG_TwoBodyDecay_H */


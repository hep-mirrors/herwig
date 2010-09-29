// -*- C++ -*-
//
// DecayVertex.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DecayVertex_H
#define HERWIG_DecayVertex_H
//
// This is the declaration of the DecayVertex class.
//
#include <ThePEG/EventRecord/HelicityVertex.h>
#include "DecayMatrixElement.h"
#include "DecayVertex.fh"

namespace Herwig {

/** \ingroup Helicity
 *  \author Peter Richardson
 *
 *  The DecayVertex class is designed to implement the vertex
 *  for a decay for use with the spin correlation algorithm. 
 *  It inherits from HelicityVertex class of ThePEG and implements 
 *  the methods to calculate the \f$\rho\f$ and \f$D\f$ matrices.
 *
 *  It uses the DecayMatrixElement class to store the matrix element and
 *  it is this class which performs the calculations of the matrices.
 *
 *  @see HelicityVertex
 *  @see DecayMatrixElement
 */ 

class DecayVertex: public HelicityVertex {
  
public:
  
  /**
   *  Access to the matrix element
   */
  //@{
  /**
   * Get the matrix element
   */
  const DecayMatrixElement & ME() const { 
    return _matrixelement;
  }

  /**
   * Set the matrix element
   */
  void ME(const DecayMatrixElement & in) const {
    _matrixelement.reset(in);
  }
  //@}  

public:
  
  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();
  
public:
  
  /**
   * Method to calculate the \f$\rho\f$ matrix for one of the decay products
   * @param iprod The product we are calculating the \f$\rho\f$ matrix for.
   * @param recursive Whether or not to recursive calculate the matrix
   */
  virtual RhoDMatrix getRhoMatrix(int iprod,bool recursive) const;

  /**
   * Method to calculate the \f$D\f$ matrix for the decaying particle. It this
   * case the argument is a dummy.
   */
  virtual RhoDMatrix getDMatrix(int) const;
  
private:
  
  /**
   * Describe a concrete class without persistent data.
   */
  static NoPIOClassDescription<DecayVertex> initDecayVertex;
  
  /** 
   * Private and non-existent assignment operator.
   */
  DecayVertex & operator=(const DecayVertex &);
  
private:
  
  /**
   * Storage of the decay matrix element.
   */
  DecayMatrixElement _matrixelement;
  
};

}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */
  
/**
 * The following template specialization informs ThePEG about the
 * base class of DecayVertex.
 */
template <>
struct BaseClassTrait<Herwig::DecayVertex,1> {
  /** Typedef of the base class of DecayVertex. */
  typedef ThePEG::HelicityVertex NthBase;
};
  
/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::DecayVertex>
  : public ClassTraitsBase<Herwig::DecayVertex> {
  
  /**
   * Return the class name.
   */
  static string className() { return "Herwig::DecayVertex"; }
};

/** @endcond */
  
}

#endif /* HERWIG_DecayVertex_H */

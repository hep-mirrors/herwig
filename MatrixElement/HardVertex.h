// -*- C++ -*-
//
// HardVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_HardVertex_H
#define HERWIG_HardVertex_H
//
// This is the declaration of the HardVertex class.

#include "ThePEG/EventRecord/HelicityVertex.h"
#include "ProductionMatrixElement.h"
#include "HardVertex.fh"
// #include "HardVertex.xh"

namespace Herwig {

using namespace ThePEG;
    
/** \ingroup Helicity
 *  \author Peter Richardson
 *
 *  The HardVertex class is designed to implement the vertex for a 
 *  hard interaction for the Herwig spin correlation algorithm. 
 *  It inherits from the HelicityVertex class of ThePEG and implements 
 *  the methods to calculate the \f$\rho\f$ and \f$D\f$ matrices.
 * 
 *  The ProductionMatrixElement class is used to store the matrix element
 *  and this class performs the calculations of the matrices.
 *
 *  @see HelicityVertex
 *  @see ProductionMatrixElement
 */ 

class HardVertex: public HelicityVertex {
  
public:
  
  /**
   *  Access to the matrix element
   */
  //@{
  /**
   * Get the matrix element
   */
  const ProductionMatrixElement & ME() const {
    return _matrixelement;
  }

  /**
   * Set the matrix element
   */
  void ME(const ProductionMatrixElement & in) const {
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
   * Method to calculate the \f$\rho\f$ matrix for one of the outgoing particles
   * @param iout The outgoing particle we are calculating the \f$\rho\f$ matrix for.
   */
  virtual RhoDMatrix getRhoMatrix(int iout,bool) const;

  /**
   * Method to calculate the \f$D\f$ matrix for an incoming particle.
   * @param in The incoming particle we are calculating the \f$D\f$ matrix for.
   */
  virtual RhoDMatrix getDMatrix(int in) const;
  
private:
  
  /**
   * Describe a concrete class without persistent data.
   */
  static NoPIOClassDescription<HardVertex> initHardVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  HardVertex & operator=(const HardVertex &) = delete;
  
private:
  
  /**
   * Storage of the matrix element.
   */
  ProductionMatrixElement _matrixelement;
  
};
}



namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */
  
/**
 * The following template specialization informs ThePEG about the
 * base class of HardVertex.
 */
template <>
struct BaseClassTrait<Herwig::HardVertex,1> {
  /** Typedef of the base class of HardVertex. */
  typedef ThePEG::HelicityVertex NthBase;
};
  
/**  
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::HardVertex>
  : public ClassTraitsBase<Herwig::HardVertex> {
  
  /**
   * Return the class name.
   */
  static string className() { return "Herwig::HardVertex"; }
};

/** @endcond */
  
}

#endif /* HERWIG_HardVertex_H */

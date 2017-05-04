// -*- C++ -*-
//
// ShowerVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ShowerVertex_H
#define HERWIG_ShowerVertex_H
//
// This is the declaration of the ShowerVertex class.
//

#include "ThePEG/EventRecord/HelicityVertex.h"
#include "Herwig/Decay/DecayMatrixElement.h"
#include "ShowerVertex.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Helicity
 *  \author Peter Richardson
 *
 * The ShowerVertex class is designed to implement the vertex for a branching
 * in the shower for use with the spin correlation alogorithm.
 *  It inherits from HelicityVertex class of ThePEG and implements 
 *  the methods to calculate the \f$\rho\f$ and \f$D\f$ matrices.
 *
 *  It uses the DecayMatrixElement class to store the matrix element and
 *  it is this class which performs the calculations of the matrices.
 *
 *  @see HelicityVertex
 *  @see DecayMatrixElement
 *
 * @see \ref ShowerVertexInterfaces "The interfaces"
 * defined for ShowerVertex.
 */
class ShowerVertex: public HelicityVertex {

public:

  /**
   * Default constructor
   */
  ShowerVertex() : convertIn_(false), convertOut_(2,false), outMatrix_(2,RhoDMatrix())
  {}

public:

  /**
   *  Access to the matrix element
   */
  //@{
  /**
   * Get the matrix element
   */
  inline const DecayMEPtr ME() const {
    return matrixElement_;
  }

  /**
   * Set the matrix element
   */
  inline void ME(DecayMEPtr in) const {
    matrixElement_ = in;
  }
  //@}

  /**
   *  Set the change of basis for the incoming particle
   */
  inline void incomingBasisTransform(RhoDMatrix conv) {
    convertIn_ = true;
    inMatrix_ = conv;
  }

public:
  
  /**
   * Method to calculate the \f$\rho\f$ matrix for one of the decay products
   * @param iprod The product we are calculating the \f$\rho\f$ matrix for.
   */
  virtual RhoDMatrix getRhoMatrix(int iprod, bool ) const;

  /**
   * Method to calculate the \f$D\f$ matrix for the decaying particle. It this
   * case the argument is a dummy.
   */
  virtual RhoDMatrix getDMatrix(int) const;

  /**
   *  Get the transformed rho matrix
   */ 
  RhoDMatrix mapIncoming(RhoDMatrix rho) const;

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ShowerVertex & operator=(const ShowerVertex &);

private:
  
  /**
   * Storage of the decay matrix element.
   */
  mutable DecayMEPtr matrixElement_;

  /**
   *  Whether or not the incoming spin density matrices need to be converted
   */
  bool convertIn_;

  /**
   *  Whether or not the outgoing spin density matrices need to be converted
   */
  vector<bool> convertOut_;

  /**
   * Storage of conversion for incoming particle
   */
  RhoDMatrix inMatrix_;
  

  /**
   *  Storage of conversion for outgoing particles
   */
  vector<RhoDMatrix> outMatrix_;

};
}

#endif /* HERWIG_ShowerVertex_H */

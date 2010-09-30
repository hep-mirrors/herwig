// -*- C++ -*-
//
// ShowerVertex.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ShowerVertex_H
#define HERWIG_ShowerVertex_H
//
// This is the declaration of the ShowerVertex class.
//

#include "ThePEG/EventRecord/HelicityVertex.h"
#include "Herwig++/Decay/DecayMatrixElement.h"
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
   *  Access to the matrix element
   */
  //@{
  /**
   * Get the matrix element
   */
  inline const DecayMatrixElement & ME() const;

  /**
   * Set the matrix element
   */
  inline void ME(const DecayMatrixElement &) const;
  //@}  

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
   * The static object used to initialize the description of this class.
   * Indicates that this is an concrete class without persistent data.
   */
  static NoPIOClassDescription<ShowerVertex> initShowerVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ShowerVertex & operator=(const ShowerVertex &);

private:
  
  /**
   * Storage of the decay matrix element.
   */
  DecayMatrixElement _matrixelement;

};
}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ShowerVertex. */
template <>
struct BaseClassTrait<Herwig::ShowerVertex,1> {
  /** Typedef of the first base class of ShowerVertex. */
  typedef Herwig::HelicityVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ShowerVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ShowerVertex>
  : public ClassTraitsBase<Herwig::ShowerVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::ShowerVertex"; }
};

/** @endcond */

}

#include "ShowerVertex.icc"

#endif /* HERWIG_ShowerVertex_H */

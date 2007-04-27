// -*- C++ -*-
#ifndef HERWIG_DecayVertex_H
#define HERWIG_DecayVertex_H
//
// This is the declaration of the DecayVertex class.
//
#include <ThePEG/Helicity/HelicityVertex.h>
#include "DecayMatrixElement.h"
// #include "DecayVertex.fh"
// #include "DecayVertex.xh"

namespace Herwig {
namespace Helicity {

using ThePEG::Helicity::HelicityVertex;

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
  inline const DecayMatrixElement & ME() const;

  /**
   * Set the matrix element
   */
  inline void ME(const DecayMatrixElement &) const;
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
   */
  virtual RhoDMatrix getRhoMatrix(int iprod);

  /**
   * Method to calculate the \f$D\f$ matrix for the decaying particle. It this
   * case the argument is a dummy.
   */
  virtual RhoDMatrix getDMatrix(int);
  
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
}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */
  
/**
 * The following template specialization informs ThePEG about the
 * base class of DecayVertex.
 */
template <>
struct BaseClassTrait<Herwig::Helicity::DecayVertex,1> {
  /** Typedef of the base class of DecayVertex. */
  typedef Herwig::Helicity::HelicityVertex NthBase;
};
  
/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::Helicity::DecayVertex>
  : public ClassTraitsBase<Herwig::Helicity::DecayVertex> {
  
  /**
   * Return the class name.
   */
  static string className() { return "Herwig++::DecayVertex"; }
};

/** @endcond */
  
}

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

  ThePEG_DECLARE_CLASS_POINTERS(Herwig::Helicity::DecayVertex,DVertexPtr);

}
}

#include "DecayVertex.icc"

#endif /* HERWIG_DecayVertex_H */

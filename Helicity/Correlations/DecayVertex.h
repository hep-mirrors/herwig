// -*- C++ -*-
#ifndef HERWIG_DecayVertex_H
#define HERWIG_DecayVertex_H
//
// This is the declaration of the DecayVertex class.

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
 *  It inherits from elicityVertex class of ThePEG and implements 
 *  the methods to calculate the rho and D matrices.
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
   * Standard ctors and dtor.
   */
  inline DecayVertex();
  inline DecayVertex(const DecayVertex &);
  virtual ~DecayVertex();
  
public:
  
  /**
   * Set and get the matrix element.
   */
  inline const DecayMatrixElement & ME() const;
  inline void ME(const DecayMatrixElement &) const;
  
public:
  
  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();
  
public:
  
  /**
   * Methods to calculate the rho and D matrices.
   */
  virtual RhoDMatrix getRhoMatrix(int);

  /**
   * Mthod to get the rho matrix for a given outgoing particle
   * and to get the D matrix for an incoming particle.
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
  
  /**
   * The following template specialization informs ThePEG about the
   * base class of DecayVertex.
   */
  template <>
  struct BaseClassTrait<Herwig::Helicity::DecayVertex,1> {
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
    static string className() { return "/Herwig++/Helicity/DecayVertex"; }

    /**
     * Return the name of the shared library to be loaded to get
     * access to this class and every other class it uses
     * (except the base class).
     */
    static string library() { return "HwCorrelations.so"; }

  };
  
}

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

  ThePEG_DECLARE_CLASS_POINTERS(Herwig::Helicity::DecayVertex,DVertexPtr);

}
}

#include "DecayVertex.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DecayVertex.tcc"
#endif

#endif /* HERWIG_DecayVertex_H */

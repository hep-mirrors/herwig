// -*- C++ -*-
#ifndef HERWIG_HardVertex_H
#define HERWIG_HardVertex_H
//
// This is the declaration of the HardVertex class.

#include "ThePEG/Helicity/HelicityVertex.h"
#include "ProductionMatrixElement.h"
// #include "HardVertex.fh"
// #include "HardVertex.xh"

namespace Herwig {
using ThePEG::Helicity::HelicityVertex;

namespace Helicity {
using namespace ThePEG;
    
/** \ingroup Helicity
 *  \author Peter Richardson
 *
 *  The HardVertex class is designed to implement the vertex for a 
 *  hard interaction for the Herwig++ spin correlation algorithm. 
 *  It inherits from the HelicityVertex class of ThePEG and implements 
 *  the methods to calculate the rho and D matrices.
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
   * Standard ctors and dtor.
   */
  inline HardVertex();
  inline HardVertex(const HardVertex &);
  virtual ~HardVertex();
  
public:
  
  /**
   * Set and get the matrix element.
   */
  inline const ProductionMatrixElement & ME() const;
  inline void ME(const ProductionMatrixElement &) const;
  
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
   * Method to get the rho matrix for a given outgoing particle
   * and the D matrix for an incoming particle.
   */
  virtual RhoDMatrix getDMatrix(int);
  
private:
  
  /**
   * Describe a concrete class without persistent data.
   */
  static NoPIOClassDescription<HardVertex> initHardVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  HardVertex & operator=(const HardVertex &);
  
private:
  
  ProductionMatrixElement _matrixelement;
  vector<double> _rubbish;
  
};

}
}


namespace ThePEG {
  
  /**
   * The following template specialization informs ThePEG about the
   * base class of HardVertex.
   */
  template <>
  struct BaseClassTrait<Herwig::Helicity::HardVertex,1> {
    typedef ThePEG::Helicity::HelicityVertex NthBase;
  };
  
  /**  
   * The following template specialization informs ThePEG about the
   * name of this class and the shared object where it is defined.
   */
  template <>
  struct ClassTraits<Herwig::Helicity::HardVertex>
    : public ClassTraitsBase<Herwig::Helicity::HardVertex> {

    /**
     * Return the class name.
     */
    static string className() { return "/Herwig++/Helicity/HardVertex"; }

    /**
     * Return the name of the shared library to be loaded to get
     * access to this class and every other class it uses
     * (except the base class).
     */
    static string library() { return "HwCorrelations.so"; }

  };
  
}

#include "HardVertex.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "HardVertex.tcc"
#endif

#endif /* HERWIG_HardVertex_H */

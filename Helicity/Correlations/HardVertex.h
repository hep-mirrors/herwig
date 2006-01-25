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
  
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline HardVertex();

  /**
   * Copy-constructor.
   */
  inline HardVertex(const HardVertex &);

  /**
   * Destructor.
   */
  virtual ~HardVertex();
  //@}
  
public:
  
  /**
   *  Access to the matrix element
   */
  //@{
  /**
   * Get the matrix element
   */
  inline const ProductionMatrixElement & ME() const;

  /**
   * Set the matrix element
   */
  inline void ME(const ProductionMatrixElement &) const;
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
  virtual RhoDMatrix getRhoMatrix(int iout);

  /**
   * Method to calculate the \f$D\f$ matrix for an incoming particle.
   * @param in The incoming particle we are calculating the \f$D\f$ matrix for.
   */
  virtual RhoDMatrix getDMatrix(int in);
  
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
  
  /**
   * Storage of the matrix element.
   */
  ProductionMatrixElement _matrixelement;
  
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
    /** Typedef of the base class of HardVertex. */
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
    static string className() { return "Herwig++::HardVertex"; }

    /**
     * Return the name of the shared library to be loaded to get
     * access to this class and every other class it uses
     * (except the base class).
     */
    static string library() { return "libHwCorrelations.so"; }

  };
  
}

#include "HardVertex.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "HardVertex.tcc"
#endif

#endif /* HERWIG_HardVertex_H */

// -*- C++ -*-
#ifndef HERWIG_ProductionMatrixElement_H
#define HERWIG_ProductionMatrixElement_H
//
// This is the declaration of the ProductionMatrixElement class.

#include <ThePEG/Config/ThePEG.h>
#include <ThePEG/Utilities/ClassDescription.h>
#include <ThePEG/Helicity/RhoDMatrix.h>
// #include "ProductionMatrixElement.fh"
// #include "ProductionMatrixElement.xh"

namespace Herwig {
using ThePEG::Helicity::RhoDMatrix;
namespace Helicity {

using namespace ThePEG;

/** \ingroup Helicity
 *  \author Peter Richardson
 *
 *  The storage of the helicity amplitude expression for the matrix element 
 *  of a hard process. Two incoming particles and an arbitary number of 
 *  external particles are supported.
 *
 *  @see DecayMatrixElement
 *  @see RhoDMatrix
 *  @see HardVertex
 */

class ProductionMatrixElement: public Base {  
      
public:
  
  // constructor for 2-2 scattering	 
  inline ProductionMatrixElement(int,int,int,int);

  // constructor for 2-3 scattering			    
  inline ProductionMatrixElement(int,int,int,int,int);

  // constructor for 2-4 scattering			    
  inline ProductionMatrixElement(int,int,int,int,int,int);

  // constructor for 2-5 scattering		    
  inline ProductionMatrixElement(int,int,int,int,int,int,int);

  // constructor for 2-6 scattering
  inline ProductionMatrixElement(int,int,int,int,int,int,int,int);

  // constructor for 2-n body process
  inline ProductionMatrixElement(int,int,vector<int>);

  // Standard ctors and dtor.
  inline ProductionMatrixElement();
  inline ProductionMatrixElement(const ProductionMatrixElement &);
  virtual ~ProductionMatrixElement();

public:
      
  // get the spins of the incoming particles particle
  inline vector<int> inspin();

  // get the spins of the outgoing particles
  inline vector<int> outspin();
  
public:

  /**
   * Access to the individual helicity components.
   */

  /**
   * 2-2 scattering.
   */
  inline Complex   operator () (int,int,int,int) const;
  inline Complex & operator () (int,int,int,int);

  /**
   * 2-3 scattering.
   */
  inline Complex   operator () (int,int,int,int,int) const;
  inline Complex & operator () (int,int,int,int,int);

  /**
   * 2-4 scattering.
   */
  inline Complex   operator () (int,int,int,int,int,int) const;
  inline Complex & operator () (int,int,int,int,int,int);

  /**
   * 2-5 scattering.
   */
  inline Complex   operator () (int,int,int,int,int,int,int) const;
  inline Complex & operator () (int,int,int,int,int,int,int);

  /**
   * 2-6 scattering.
   */
  inline Complex   operator () (int,int,int,int,int,int,int,int) const;
  inline Complex & operator () (int,int,int,int,int,int,int,int);

  /**
   * 2-n scattering.
   */
  inline Complex   operator () (vector<int>) const;
  inline Complex & operator () (vector<int>);

public:

  /**
   * Calculate the decay matrix for an incoming particle.
   */
  RhoDMatrix calculateDMatrix(int,RhoDMatrix, vector<RhoDMatrix>);

  /**
   * Calculate the rho matrix for a given outgoing particle.
   */
  RhoDMatrix calculateRhoMatrix(int,RhoDMatrix,RhoDMatrix,
				vector<RhoDMatrix>);
  
public:

  /**
   * Reset the matrix element.
   */
  inline void reset(const ProductionMatrixElement &) const;

public:
  
  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();
  
private:
  
  /**
   * Describe a concrete class without persistent data.
   */
  static NoPIOClassDescription<ProductionMatrixElement> initProductionMatrixElement;
  
  /**
   * Private and non-existent assignment operator.
   */
  ProductionMatrixElement & operator=(const ProductionMatrixElement& );
  
private:
  
  /**
   * Set the size of the vector containing the matrix element.
   */
  inline void setMESize();
  
private:
  
  /**
   * Number of outgoing particles.
   */
  mutable int _nout;

  /**
   * Spin of the incoming particles as 2s+1.
   */
  mutable vector<int> _inspin;

  /**
   * Spins of the outgoing particles.
   */
  mutable vector<int> _outspin;

  /**
   * Storage of the matrix element, a vector is better for memory usage.
   */
  mutable vector<Complex> _matrixelement;

  /**
   * Constants needed to map the index of the vector to a helicity structure.
   */
  mutable vector<int> _constants;

};

}
}


namespace ThePEG {

  /**
   * The following template specialization informs ThePEG about the
   * base class of ProductionMatrixElement.
   */
  template <>
  struct BaseClassTrait<Herwig::Helicity::ProductionMatrixElement,1> {
    typedef Base NthBase;
  };
  
  /**
   * The following template specialization informs ThePEG about the
   * name of this class and the shared object where it is defined.
   */
  template <>
  struct ClassTraits<Herwig::Helicity::ProductionMatrixElement>
    : public ClassTraitsBase<Herwig::Helicity::ProductionMatrixElement> {

    /**
     * Return the class name.
     */
    static string className() { return "/Herwig++/Helicity/ProductionMatrixElement"; }

    /**
     * Return the name of the shared library to be loaded to get
     * access to this class and every other class it uses
     * (except the base class).
     */
    static string library() { return "hwCorrelations.so"; }

  };
  
}

#include "ProductionMatrixElement.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ProductionMatrixElement.tcc"
#endif

#endif /* HERWIG_ProductionMatrixElement_H */

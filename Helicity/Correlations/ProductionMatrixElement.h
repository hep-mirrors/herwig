// -*- C++ -*-
#ifndef HERWIG_ProductionMatrixElement_H
#define HERWIG_ProductionMatrixElement_H
//
// This is the declaration of the <!id>ProductionMatrixElement<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  The storage of the helicity amplitude expression for the matrix element of
//  a hard process. Two incoming particles and an arbitary number of external particles
//  are supported.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="DecayMatrixElement.html">DecayMatrixElement.h</a>,
// <a href="RhoDMatrix.html">RhoDMatrix.h</a>,
// <a href="HardVertex.html">HardVertex.h</a>.
// 
// Author: Peter Richardson
//

#include <ThePEG/Config/ThePEG.h>
#include <ThePEG/Utilities/ClassDescription.h>
#include <ThePEG/Helicity/RhoDMatrix.h>
// #include "ProductionMatrixElement.fh"
// #include "ProductionMatrixElement.xh"

namespace Herwig {
namespace Helicity {

using namespace ThePEG;
using ThePEG::Helicity::RhoDMatrix;

class ProductionMatrixElement: public Base {  
      
public:
  
  inline ProductionMatrixElement(int,int,int,int);
  // constructor for 2-2 scattering	 

  inline ProductionMatrixElement(int,int,int,int,int);
  // constructor for 2-3 scattering			    

  inline ProductionMatrixElement(int,int,int,int,int,int);
  // constructor for 2-4 scattering			    

  inline ProductionMatrixElement(int,int,int,int,int,int,int);
  // constructor for 2-5 scattering		    

  inline ProductionMatrixElement(int,int,int,int,int,int,int,int);
  // constructor for 2-6 scattering

  inline ProductionMatrixElement(int,int,vector<int>);
  // constructor for 2-n body process

  inline ProductionMatrixElement();
  inline ProductionMatrixElement(const ProductionMatrixElement &);
  virtual ~ProductionMatrixElement();
  // Standard ctors and dtor.

public:
      
  inline vector<int> inspin();
  // get the spins of the incoming particles particle

  inline vector<int> outspin();
  // get the spins of the outgoing particles
  
public:

  // access to the individual helicity components
  inline Complex   operator () (int,int,int,int) const;
  inline Complex & operator () (int,int,int,int);
  // 2-2 scattering
  inline Complex   operator () (int,int,int,int,int) const;
  inline Complex & operator () (int,int,int,int,int);
  // 2-3 scattering
  inline Complex   operator () (int,int,int,int,int,int) const;
  inline Complex & operator () (int,int,int,int,int,int);
  // 2-4 scattering
  inline Complex   operator () (int,int,int,int,int,int,int) const;
  inline Complex & operator () (int,int,int,int,int,int,int);
  // 2-5 scattering
  inline Complex   operator () (int,int,int,int,int,int,int,int) const;
  inline Complex & operator () (int,int,int,int,int,int,int,int);
  // 2-6 scattering
  inline Complex   operator () (vector<int>) const;
  inline Complex & operator () (vector<int>);
  // 2-n scattering

public:

  RhoDMatrix calculateDMatrix(int,RhoDMatrix, vector<RhoDMatrix>);
  // calculate the decay matrix for an incoming particle

  RhoDMatrix calculateRhoMatrix(int,RhoDMatrix,RhoDMatrix,
				vector<RhoDMatrix>);
  // calculate the rho matrix for a given outgoing particle
  
public:

  inline void reset(const ProductionMatrixElement &) const;
  // reset the matrix element

public:
  
  static void Init();
  // Standard Init function used to initialize the interfaces.
  
private:
  
  static NoPIOClassDescription<ProductionMatrixElement> initProductionMatrixElement;
  // Describe a concrete class without persistent data.
  
  ProductionMatrixElement & operator=(const ProductionMatrixElement& );
  // Private and non-existent assignment operator.
  
private:
  
  inline void setMESize();
  // set the size of the vector containing the matrix element
  
private:
  
  mutable int _nout;
  // number of outgoing particles
  mutable vector<int> _inspin;
  // spin of the incoming particles as 2s+1
  mutable vector<int> _outspin;
  // spins of the outgoing particles
  mutable vector<Complex> _matrixelement;
  // storage of the matrix element, a vector is better for memory usage
  mutable vector<int> _constants;
  // constants needed to map the index of the vector to a helicity structure 
};

}
}



// CLASSDOC OFF

namespace ThePEG {

  // The following template specialization informs ThePEG about the
  // base class of ProductionMatrixElement.
  template <>
  struct BaseClassTrait<Herwig::Helicity::ProductionMatrixElement,1> {
    typedef Base NthBase;
  };
  
  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::Helicity::ProductionMatrixElement>
    : public ClassTraitsBase<Herwig::Helicity::ProductionMatrixElement> {
    static string className() { return "/Herwig++/Helicity/ProductionMatrixElement"; }
    // Return the class name.
    static string library() { return "hwCorrelations.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
    // (except the base class).
  };
  
}

#include "ProductionMatrixElement.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ProductionMatrixElement.tcc"
#endif

#endif /* HERWIG_ProductionMatrixElement_H */

// -*- C++ -*-
#ifndef HERWIG_HardVertex_H
#define HERWIG_HardVertex_H
//
// This is the declaration of the <!id>HardVertex<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  The  <!id>HardVertex<!!id> class is designed to implement the vertex for a 
//  hard interaction for the Herwig++ spin correlation algorithm. It inherits from the
//  <!id>HelicityVertex<!!id> class of ThePEG and implements the methods to calculate
//  the rho and D matrices.
//
//  The <!id>ProductionMatrixElement<!!id> class is used to store the matrix element
//  and this class performs the calculations of the matrices.
//
//
// CLASSDOC SUBSECTION See also:
//
// <a href="HelicityVertex.html">HelicityVertix.h</a>,
// <a href="ProductionMatrixElement.html">ProductionMatrixElement.h</a>.
// 

#include "ThePEG/Helicity/HelicityVertex.h"
#include "ProductionMatrixElement.h"
// #include "HardVertex.fh"
// #include "HardVertex.xh"

namespace Herwig {
using ThePEG::Helicity::HelicityVertex;

namespace Helicity {
using namespace ThePEG;
    
class HardVertex: public HelicityVertex {
  
public:
  
  inline HardVertex();
  inline HardVertex(const HardVertex &);
  virtual ~HardVertex();
  // Standard ctors and dtor.
  
public:
  
  inline const ProductionMatrixElement & ME() const;
  inline void ME(const ProductionMatrixElement &) const;
  // set and get the matrix element
  
public:
  
  static void Init();
  // Standard Init function used to initialize the interfaces.
  
public:
  
  // methods to calculate the rho and D matrices
  virtual RhoDMatrix getRhoMatrix(int);
  // method to get the rho matrix for a given outgoing particle
  virtual RhoDMatrix getDMatrix(int);
  // method to get the D matrix for an incoming particle
  
private:
  
  static NoPIOClassDescription<HardVertex> initHardVertex;
  // Describe a concrete class without persistent data.
  
  HardVertex & operator=(const HardVertex &);
  // Private and non-existent assignment operator.
  
private:
  
  ProductionMatrixElement _matrixelement;
  vector<double> _rubbish;
  
};

}
}

// CLASSDOC OFF

namespace ThePEG {
  
  // The following template specialization informs ThePEG about the
  // base class of HardVertex.
  template <>
  struct BaseClassTrait<Herwig::Helicity::HardVertex,1> {
    typedef ThePEG::Helicity::HelicityVertex NthBase;
  };
  
  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::Helicity::HardVertex>
    : public ClassTraitsBase<Herwig::Helicity::HardVertex> {
    static string className() { return "/Herwig++/Helicity/HardVertex"; }
    // Return the class name.
    static string library() { return "HwCorrelations.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
    // (except the base class).
  };
  
}

#include "HardVertex.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "HardVertex.tcc"
#endif

#endif /* HERWIG_HardVertex_H */

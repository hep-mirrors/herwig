// -*- C++ -*-
#ifndef HERWIG_DecayVertex_H
#define HERWIG_DecayVertex_H
//
// This is the declaration of the <!id>DecayVertex<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  The <!id>DecayVertex<!!id> class is designed to implement the vertex
//  for a decay for use with the spin correlation algorithm. It inherits from
//  the <!id>HelicityVertex<!!id> class of ThePEG and implements the methods to
//  calculate the rho and D matrices.
//
//  It uses the <!id>DecayMatrixElement<!!id> class to store the matrix element and
//  it is this class which performs the calculations of the matrices.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="HelicityVertex.html">HelicityVertix.h</a>,
// <a href="DecayMatrixElement.html">DecayMatrixElement.h</a>.
// 
// Author: Peter Richardson
//

#include <ThePEG/Helicity/HelicityVertex.h>
#include "DecayMatrixElement.h"
// #include "DecayVertex.fh"
// #include "DecayVertex.xh"

namespace Herwig {
namespace Helicity {

using ThePEG::Helicity::HelicityVertex;

class DecayVertex: public HelicityVertex {
      
public:
  
  inline DecayVertex();
  inline DecayVertex(const DecayVertex &);
  virtual ~DecayVertex();
  // Standard ctors and dtor.
  
public:
  
  inline const DecayMatrixElement & ME() const;
  inline void ME(const DecayMatrixElement &) const;
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
  
  static NoPIOClassDescription<DecayVertex> initDecayVertex;
  // Describe a concrete class without persistent data.
  
  DecayVertex & operator=(const DecayVertex &);
  // Private and non-existent assignment operator.
  
private:
  
  // storage of the decay matrix element
  DecayMatrixElement _matrixelement;
  
};

}
}

// CLASSDOC OFF

namespace ThePEG {
  
  // The following template specialization informs ThePEG about the
  // base class of DecayVertex.
  template <>
  struct BaseClassTrait<Herwig::Helicity::DecayVertex,1> {
    typedef Herwig::Helicity::HelicityVertex NthBase;
  };
  
  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::Helicity::DecayVertex>
    : public ClassTraitsBase<Herwig::Helicity::DecayVertex> {
    static string className() { return "/Herwig++/Helicity/DecayVertex"; }
    // Return the class name.
    static string library() { return "HwCorrelations.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
    // (except the base class).
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

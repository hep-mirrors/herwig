// -*- C++ -*-
#ifndef HERWIG_SMGGGVertex_H
#define HERWIG_SMGGGVertex_H
//
// This is the declaration of the <!id>SMGGGVertex<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  Implementation of the SM triple gluon vertex
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:.html">.h</a>,
// <a href="http:.html">.h</a>.
// 

#include "Herwig++/Helicity/Vertex/Vector/VVVVertex.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Utilities/Rebinder.h"


namespace Herwig {
namespace Helicity{
using namespace ThePEG;

class SMGGGVertex: public VVVVertex {

public:

  inline SMGGGVertex();
  inline SMGGGVertex(const SMGGGVertex &);
  virtual ~SMGGGVertex();
  // Standard ctors and dtor.
  
public:
  
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.
  
  static void Init();
  // Standard Init function used to initialize the interfaces.
  
  void setCoupling(Energy2,tcPDPtr, tcPDPtr, tcPDPtr);
  // calculate the coupling

protected:
  
  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;
  // Standard clone methods.
  
protected:
      
  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void doinitrun();
  inline virtual void dofinish();
  // Standard Interfaced virtual functions.
  
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  // Change all pointers to Interfaced objects to corresponding clones.
  
  inline virtual IVector getReferences();
  // Return pointers to all Interfaced objects refered to by this.
  
private:
  
  static ClassDescription<SMGGGVertex> initSMGGGVertex;
  // Describe a concrete class with persistent data.
  
  SMGGGVertex & operator=(const SMGGGVertex &);
  // Private and non-existent assignment operator.

  // pointer to the Standard Model
  SMPtr _theSM;
  // storage of the couplings
  Complex _couplast;
  Energy2 _q2last;
  
};

}
}
#include "SMGGGVertex.icc"

// CLASSDOC OFF

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of SMGGGVertex.
template <>
struct BaseClassTrait<Herwig::Helicity::SMGGGVertex,1> {
  typedef Herwig::Helicity::VVVVertex NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::Helicity::SMGGGVertex>
  : public ClassTraitsBase<Herwig::Helicity::SMGGGVertex> {
  static string className() { return "/Herwig++/Helicity/SMGGGVertex"; }
  // Return the class name.
  static string library() { return "libHwSMVertex.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#endif /* HERWIG_SMGGGVertex_H */

// -*- C++ -*-
#ifndef HERWIG_SMFFPVertex_H
#define HERWIG_SMFFPVertex_H
//
// This is the declaration of the <!id>SMFFPVertex<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  This is the implementation of the Standard Model fermion-antifermion photon vertex
//
// CLASSDOC SUBSECTION See also:
//
// <a href="FFVVertex.html">FFVVertex.h</a>,
// <a href="VertexBase.html">VertexBase.h</a>.
// 

#include "Herwig++/Helicity/Vertex/Vector/FFVVertex.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Utilities/Rebinder.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

class SMFFPVertex: public FFVVertex {
  
public:
  
  inline SMFFPVertex();
  inline SMFFPVertex(const SMFFPVertex &);
  virtual ~SMFFPVertex();
  // Standard ctors and dtor.

public:
  
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.
  
  static void Init();
  // Standard Init function used to initialize the interfaces.
  
  void setCoupling(Energy2,tcPDPtr, tcPDPtr, tcPDPtr);
  // calculate the couplings
  
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
  
  static ClassDescription<SMFFPVertex> initSMFFPVertex;
  // Describe a concrete class with persistent data.
  
  SMFFPVertex & operator=(const SMFFPVertex &);
  // Private and non-existent assignment operator.
  
private:
  // pointer to the Standard Model object
  SMPtr _theSM;
  // storage of the couplings
  double _charge[17];
  Complex _couplast;
  Energy2 _q2last;
};

}
}
#include "SMFFPVertex.icc"

// CLASSDOC OFF

namespace ThePEG {

  // The following template specialization informs ThePEG about the
  // base class of SMFFPVertex.
  template <>
  struct BaseClassTrait<Herwig::Helicity::SMFFPVertex,1> {
    typedef Herwig::Helicity::FFVVertex NthBase;
  };
  
  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::Helicity::SMFFPVertex>
    : public ClassTraitsBase<Herwig::Helicity::SMFFPVertex> {
    static string className() { return "/Herwig++/Helicity/SMFFPVertex"; }
    // Return the class name.
    static string library() { return "libHwSMVertex.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
  // (except the base class).
  };
  
}


#endif /* HERWIG_SMFFPVertex_H */

// -*- C++ -*-
#ifndef HERWIG_SMFFWVertex_H
#define HERWIG_SMFFWVertex_H
//
// This is the declaration of the <!id>SMFFWVertex<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This is the implementation of the Standard model coupling of the W to the
// fermions.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="FFVVertex.html">FFVVertex.h</a>,
// <a href="VertexBase.html">VertexBase.h</a>.
// 

#include "Herwig++/Helicity/Vertex/Vector/FFVVertex.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/StandardModel/CKMBase.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "Herwig++/Models/StandardModel/StandardCKM.h"
namespace Herwig {
namespace Helicity {
using namespace ThePEG;

class SMFFWVertex: public FFVVertex {
  
public:
  
  inline SMFFWVertex();
  inline SMFFWVertex(const SMFFWVertex &);
  virtual ~SMFFWVertex();
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
  
  static ClassDescription<SMFFWVertex> initSMFFWVertex;
  // Describe a concrete class with persistent data.
  
  SMFFWVertex & operator=(const SMFFWVertex &);
  // Private and non-existent assignment operator.

private:

  // pointer to the Standard Model object
  SMPtr _theSM;
  Ptr<CKMBase>::pointer _theCKM;
  // storage of the couplings
  Complex _ckm[3][3];
  Complex _couplast;
  Energy2 _q2last;
  
};

}
}
#include "SMFFWVertex.icc"

// CLASSDOC OFF

namespace ThePEG {
  
  // The following template specialization informs ThePEG about the
  // base class of SMFFWVertex.
  template <>
  struct BaseClassTrait<Herwig::Helicity::SMFFWVertex,1> {
    typedef Herwig::Helicity::FFVVertex NthBase;
  };
  
  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::Helicity::SMFFWVertex>
    : public ClassTraitsBase<Herwig::Helicity::SMFFWVertex> {
    static string className() { return "/Herwig++/Helicity/SMFFWVertex"; }
    // Return the class name.
    static string library() { return "libHwSMVertex.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
    // (except the base class).
  };
  
}

#endif /* HERWIG_SMFFWVertex_H */

// -*- C++ -*-
#ifndef HERWIG_SMFFHVertex_H
#define HERWIG_SMFFHVertex_H
//
// This is the declaration of the <!id>SMFFHVertex<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  This is the implementation of the vertex coupling the Standard Model Higgs
//  to the Standard Model fermions for helicity amplitude calculations
//
// CLASSDOC SUBSECTION See also:
//
// <a href="FFSVertex.html">FFSVertex.h</a>,
// <a href="VertexBase.html">VertexBase.h</a>.
// 

#include "Herwig++/Helicity/Vertex/Scalar/FFSVertex.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Utilities/Rebinder.h"


namespace Herwig {
namespace Helicity {
using namespace ThePEG;

class SMFFHVertex: public FFSVertex {
  
public:
  
  inline SMFFHVertex();
  inline SMFFHVertex(const SMFFHVertex &);
  virtual ~SMFFHVertex();
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
  
  static ClassDescription<SMFFHVertex> initSMFFHVertex;
  // Describe a concrete class with persistent data.
  
  SMFFHVertex & operator=(const SMFFHVertex &);
  // Private and non-existent assignment operator.

  Ptr<Herwig::StandardModel>::pointer _theSM;
  // pointer to the SM object
  Complex _couplast;
  double _sw;
  int _idlast;
  Energy2 _q2last;
  Energy _masslast,_mw;
  // storage of the couplings
};  

}
} 
#include "SMFFHVertex.icc"

// CLASSDOC OFF

namespace ThePEG {
  
  // The following template specialization informs ThePEG about the
  // base class of SMFFHVertex.
  template <>
  struct BaseClassTrait<Herwig::Helicity::SMFFHVertex,1> {
    typedef Herwig::Helicity::FFSVertex NthBase;
  };
  
  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::Helicity::SMFFHVertex>
    : public ClassTraitsBase<Herwig::Helicity::SMFFHVertex> {
    static string className() { return "/Herwig++/Helicity/SMFFHVertex"; }
    // Return the class name.
    static string library() { return "libHwSMVertex.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
    // (except the base class).
  };
  
}


#endif /* HERWIG_SMFFHVertex_H */

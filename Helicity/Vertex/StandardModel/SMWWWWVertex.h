// -*- C++ -*-
#ifndef HERWIG_SMWWWWVertex_H
#define HERWIG_SMWWWWVertex_H
//
// This is the declaration of the <!id>SMWWWWVertex<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  The <!id>SMWWWWVertex<!!id> class is the implementation of the
//  coupling of four electroweak gauge bosons in the SM. It inherits
//  from VVVVVertex nad implements the setCoupling member
//
// CLASSDOC SUBSECTION See also:
//
// <a href="VVVVVVertex.html">VVVVVertex.h</a>,
// <a href="VertexBase.html">VertexBase.h</a>.
// 
#include "Herwig++/Helicity/Vertex/Vector/VVVVVertex.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Utilities/Rebinder.h"


namespace Herwig {
namespace Helicity {
using namespace ThePEG;
    
class SMWWWWVertex: public VVVVVertex {
  
public:
  
  inline SMWWWWVertex();
  inline SMWWWWVertex(const SMWWWWVertex &);
  virtual ~SMWWWWVertex();
  // Standard ctors and dtor.
  
public:
  
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.
  
  static void Init();
  // Standard Init function used to initialize the interfaces.

  void setCoupling(Energy2,tcPDPtr, tcPDPtr, tcPDPtr,tcPDPtr);
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
  
  static ClassDescription<SMWWWWVertex> initSMWWWWVertex;
  // Describe a concrete class with persistent data.
  
  SMWWWWVertex & operator=(const SMWWWWVertex &);
  // Private and non-existent assignment operator.

  SMPtr _theSM;
  // pointer to the SM object  
  tcPDPtr _gamma,_Z0,_wplus,_wminus;
  // intermediate particles
  Complex _couplast;
  Energy2 _q2last;
  double _vfact[4],_sw2,_cw2;
  // storage of the couplings
};
 
}
}
#include "SMWWWWVertex.icc"

// CLASSDOC OFF

namespace ThePEG {
  
  // The following template specialization informs ThePEG about the
  // base class of SMWWWWVertex.
  template <>
  struct BaseClassTrait<Herwig::Helicity::SMWWWWVertex,1> {
    typedef Herwig::Helicity::VVVVVertex NthBase;
  };
  
  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::Helicity::SMWWWWVertex>
    : public ClassTraitsBase<Herwig::Helicity::SMWWWWVertex> {
    static string className() { return "/Herwig++/Helicity/SMWWWWVertex"; }
    // Return the class name.
    static string library() { return "libHwSMVertex.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
    // (except the base class).
  };
  
}


#endif /* HERWIG_SMWWWWVertex_H */

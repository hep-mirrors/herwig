// -*- C++ -*-
#ifndef HERWIG_SMWWHVertex_H
#define HERWIG_SMWWHVertex_H
//
// This is the declaration of the <!id>SMWWHVertex<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  The <!id>SMWWHVertex<!!id> is the implementation of the
//  coupling of two electroweak gauge bosons to the Higgs in the Standard
//  Model. It inherits from VVSVertex and implements the setCoupling member.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="VVSVertex.html">VVSVertex.h</a>,
// <a href="VertexBase.html">VertexBase.h</a>.
// 

#include "Herwig++/Helicity/Vertex/Scalar/VVSVertex.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/Utilities/Rebinder.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

class SMWWHVertex: public VVSVertex {
  
public:
  
  inline SMWWHVertex();
  inline SMWWHVertex(const SMWWHVertex &);
  virtual ~SMWWHVertex();
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
  
  static ClassDescription<SMWWHVertex> initSMWWHVertex;
  // Describe a concrete class with persistent data.
  
  SMWWHVertex & operator=(const SMWWHVertex &);
  // Private and non-existent assignment operator.
  
  SMPtr _theSM;
  // pointer to he Standard Model object
  Complex _couplast;
  Energy2 _q2last;
  Energy _mw;
  double _zfact,_sw;
  // storage of the couplings
  
};

}
}
#include "SMWWHVertex.icc"

// CLASSDOC OFF

namespace ThePEG {

  // The following template specialization informs ThePEG about the
  // base class of SMWWHVertex.
  template <>
  struct BaseClassTrait<Herwig::Helicity::SMWWHVertex,1> {
    typedef Herwig::Helicity::VVSVertex NthBase;
  };
  
  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::Helicity::SMWWHVertex>
    : public ClassTraitsBase<Herwig::Helicity::SMWWHVertex> {
    static string className() { return "/Herwig++/Helicity/SMWWHVertex"; }
    // Return the class name.
    static string library() { return "libHwSMVertex.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
    // (except the base class).
  };
  
}


#endif /* HERWIG_SMWWHVertex_H */

// -*- C++ -*-
#ifndef HERWIG_RSModelVVVGRVertex_H
#define HERWIG_RSModelVVVGRVertex_H
//
// This is the declaration of the <!id>RSModelVVVGRVertex<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  The <!id>RSModelVVVGRVertex<!!id> class is the implementation of the 
//  triple vector graviton couling in the RS model. It inherits from VVVTVertex
//  and implements the setCoupling member
//
// CLASSDOC SUBSECTION See also:
//
// <a href="VVVTVertex.html">VVVTVertex.h</a>,
// <a href="VertexBase.html">VertexBase.h</a>.
// 

#include "Herwig++/Helicity/Vertex/Tensor/VVVTVertex.h"
#include "Herwig++/Models/RSModel/RSModel.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;
    
class RSModelVVVGRVertex: public VVVTVertex {
  
public:
  
  inline RSModelVVVGRVertex();
  inline RSModelVVVGRVertex(const RSModelVVVGRVertex &);
  virtual ~RSModelVVVGRVertex();
  // Standard ctors and dtor.
  
public:
  
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.
  
  static void Init();
  // Standard Init function used to initialize the interfaces.

  void setCoupling(Energy2,tcPDPtr,tcPDPtr,tcPDPtr,tcPDPtr);
  // set the coupling

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
  
  static ClassDescription<RSModelVVVGRVertex> initRSModelVVVGRVertex;
  // Describe a concrete class with persistent data.
  
  RSModelVVVGRVertex & operator=(const RSModelVVVGRVertex &);
  // Private and non-existent assignment operator.

private:
  SMPtr _theModel;
  // pointer to the model object
  double _theKappa;
  // the graviton coupling
  Complex _couplast[2];
  Energy2 _q2last[2];
  double _zfact;
  // storage of the couplings
};

}
}
#include "RSModelVVVGRVertex.icc"

// CLASSDOC OFF

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of RSModelVVVGRVertex.
template <>
struct BaseClassTrait<Herwig::Helicity::RSModelVVVGRVertex,1> {
  typedef Herwig::Helicity::VVVTVertex NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::Helicity::RSModelVVVGRVertex>
  : public ClassTraitsBase<Herwig::Helicity::RSModelVVVGRVertex> {
  static string className() { return "/Herwig++/Helicity/RSModelVVVGRVertex"; }
  // Return the class name.
  static string library() { return "libHwRSVertex.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}


#endif /* HERWIG_RSModelVVVGRVertex_H */

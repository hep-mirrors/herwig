// -*- C++ -*-
#ifndef HERWIG_RSModelVVGRVertex_H
#define HERWIG_RSModelVVGRVertex_H
//
// This is the declaration of the <!id>RSModelVVGRVertex<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  This is the implementation of the vector-vector-graviton vertex for
//  the RS model
//
// CLASSDOC SUBSECTION See also:
//
// <a href="VVTVertex.html">VVTVertex.h</a>,
// <a href="VertexBase.html">VertexBase.h</a>.
// 

#include "Herwig++/Helicity/Vertex/Tensor/VVTVertex.h"
#include "Herwig++/Models/RSModel/RSModel.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

class RSModelVVGRVertex: public VVTVertex {
  
public:
  
  inline RSModelVVGRVertex();
  inline RSModelVVGRVertex(const RSModelVVGRVertex &);
  virtual ~RSModelVVGRVertex();
  // Standard ctors and dtor.
  
public:
  
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.
  
  static void Init();
  // Standard Init function used to initialize the interfaces.
  
  void setCoupling(Energy2,tcPDPtr,tcPDPtr,tcPDPtr);
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
  
  static ClassDescription<RSModelVVGRVertex> initRSModelVVGRVertex;
  // Describe a concrete class with persistent data.
  
  RSModelVVGRVertex & operator=(const RSModelVVGRVertex &);
  // Private and non-existent assignment operator.

  SMPtr _theModel;
  // pointer to the model object
  double _theKappa;
  // the coupling
  
};

}
}
#include "RSModelVVGRVertex.icc"

// CLASSDOC OFF

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of RSModelVVGRVertex.
template <>
struct BaseClassTrait<Herwig::Helicity::RSModelVVGRVertex,1> {
  typedef Herwig::Helicity::VVTVertex NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::Helicity::RSModelVVGRVertex>
  : public ClassTraitsBase<Herwig::Helicity::RSModelVVGRVertex> {
  static string className() { return "/Herwig++/Helicity/RSModelVVGRVertex"; }
  // Return the class name.
  static string library() { return "libHwRSVertex.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}


#endif /* HERWIG_RSModelVVGRVertex_H */

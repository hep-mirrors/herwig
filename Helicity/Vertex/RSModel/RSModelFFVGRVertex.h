// -*- C++ -*-
#ifndef HERWIG_RSModelFFVGRVertex_H
#define HERWIG_RSModelFFVGRVertex_H
//
// This is the declaration of the <!id>RSModelFFVGRVertex<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  This is the implementation of the fermion-antifermion-vector-graviton
//  vertex for the Randell-Sundrum model
//
// CLASSDOC SUBSECTION See also:
//
// <a href="FFVTVertex.html">FFVTVertex.h</a>,
// <a href="VertexBase.html">VertexBase.h</a>.
// 

#include "Herwig++/Helicity/Vertex/Tensor/FFVTVertex.h"
#include "Herwig++/Models/RSModel/RSModel.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

class RSModelFFVGRVertex: public FFVTVertex {
  
public:
  
  inline RSModelFFVGRVertex();
  inline RSModelFFVGRVertex(const RSModelFFVGRVertex &);
  virtual ~RSModelFFVGRVertex();
  // Standard ctors and dtor.
  
public:
  
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.
  
  static void Init();
  // Standard Init function used to initialize the interfaces.
  
  void setCoupling(Energy2,tcPDPtr,tcPDPtr,tcPDPtr,tcPDPtr);
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
  
  static ClassDescription<RSModelFFVGRVertex> initRSModelFFVGRVertex;
  // Describe a concrete class with persistent data.
  
  RSModelFFVGRVertex & operator=(const RSModelFFVGRVertex &);
  // Private and non-existent assignment operator.

private:
  // pointer to the Standard Model object
  SMPtr _theModel;
  // storage of the couplings
  double _charge[17];
  Complex _couplast[2];
  Energy2 _q2last[2];
  double _theKappa;
};

}
}
#include "RSModelFFVGRVertex.icc"

// CLASSDOC OFF

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of RSModelFFVGRVertex.
template <>
struct BaseClassTrait<Herwig::Helicity::RSModelFFVGRVertex,1> {
  typedef Herwig::Helicity::FFVTVertex NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::Helicity::RSModelFFVGRVertex>
  : public ClassTraitsBase<Herwig::Helicity::RSModelFFVGRVertex> {
  static string className() { return "/Herwig++/Helicity/RSModelFFVGRVertex"; }
  // Return the class name.
  static string library() { return "libHwRSVertex.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}


#endif /* HERWIG_RSModelFFVGRVertex_H */

// -*- C++ -*-
#ifndef HERWIG_RSModelSSGRVertex_H
#define HERWIG_RSModelSSGRVertex_H
//
// This is the declaration of the <!id>RSModelSSGRVertex<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// The <!id>RSModelSSGRVertex<!!id> class is thew implementation of the graviton
// coupling to the Higgs in the RSModel. It inherits from the SSTVertex and
// implements the setCoupling member
//
// CLASSDOC SUBSECTION See also:
//
// <a href="SSTVertex.html">STTVertex.h</a>,
// <a href="VertexBase.html">.VertexBase.h</a>.
// 

#include "Herwig++/Helicity/Vertex/Tensor/SSTVertex.h"
#include "Herwig++/Models/RSModel/RSModel.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;
    
class RSModelSSGRVertex: public SSTVertex {
  
public:
  
  inline RSModelSSGRVertex();
  inline RSModelSSGRVertex(const RSModelSSGRVertex &);
  virtual ~RSModelSSGRVertex();
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
  
  static ClassDescription<RSModelSSGRVertex> initRSModelSSGRVertex;
  // Describe a concrete class with persistent data.
  
  RSModelSSGRVertex & operator=(const RSModelSSGRVertex &);
  // Private and non-existent assignment operator.
  
  SMPtr _theModel;
  // pointer to the Mode
  double _theKappa;
  // coupling
};

}
}    
#include "RSModelSSGRVertex.icc"

// CLASSDOC OFF

namespace ThePEG {
  
  // The following template specialization informs ThePEG about the
  // base class of RSModelSSGRVertex.
  template <>
  struct BaseClassTrait<Herwig::Helicity::RSModelSSGRVertex,1> {
    typedef Herwig::Helicity::SSTVertex NthBase;
  };
  
  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::Helicity::RSModelSSGRVertex>
    : public ClassTraitsBase<Herwig::Helicity::RSModelSSGRVertex> {
    static string className() { return "/Herwig++/Helicity/RSModelSSGRVertex"; }
    // Return the class name.
    static string library() { return "libHwRSVertex.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
    // (except the base class).
  };
  
}


#endif /* HERWIG_RSModelSSGRVertex_H */

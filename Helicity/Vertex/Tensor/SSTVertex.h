// -*- C++ -*-
#ifndef HERWIG_SSTVertex_H
#define HERWIG_SSTVertex_H
//
// This is the declaration of the <!id>SSTVertex<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  The <!id>VVTVertex<!!id> class is the implementation of the scalar-scalar-tensor
//  vertex. It inherits from the VertexBase class for the storage of the particles
//  interacting at the vertex and implements the helicity amplitude calculations.
//
//  All implementations of this vertex should inherit from it and implement the
//  virtual setCoupling member.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="VertexBase.html">VertexBase.h</a>.
// 

#include "Herwig++/Helicity/Vertex/VertexBase.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/TensorWaveFunction.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

class SSTVertex: public VertexBase {
  
public:
  
  inline SSTVertex();
  inline SSTVertex(const SSTVertex &);
  virtual ~SSTVertex();
  // Standard ctors and dtor.
  
public:
  
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.
  
  static void Init();
  // Standard Init function used to initialize the interfaces.

public:
  
  Complex evaluate(Energy2, const ScalarWaveFunction &,
    		       const ScalarWaveFunction &, const TensorWaveFunction &);
  ScalarWaveFunction evaluate(Energy2,int,tcPDPtr, const ScalarWaveFunction &,
    			  const TensorWaveFunction &);
  TensorWaveFunction evaluate(Energy2,int,tcPDPtr, const ScalarWaveFunction &,
    			  const ScalarWaveFunction &);
  virtual void setCoupling(Energy2,tcPDPtr,tcPDPtr,tcPDPtr);
  
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
  
  static AbstractClassDescription<SSTVertex> initSSTVertex;
  // Describe a concrete class with persistent data.
  
  SSTVertex & operator=(const SSTVertex &);
  // Private and non-existent assignment operator.
  
};

}
}
#include "SSTVertex.icc"

// CLASSDOC OFF

namespace ThePEG {
  
  // The following template specialization informs ThePEG about the
  // base class of SSTVertex.
  template <>
  struct BaseClassTrait<Herwig::Helicity::SSTVertex,1> {
    typedef Herwig::Helicity::VertexBase NthBase;
  };
  
  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::Helicity::SSTVertex>
    : public ClassTraitsBase<Herwig::Helicity::SSTVertex> {
    static string className() { return "/Herwig++/Helicity/SSTVertex"; }
    // Return the class name.
    static string library() { return "libHwTVertex.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
    // (except the base class).
  };
  
}


#endif /* HERWIG_SSTVertex_H */

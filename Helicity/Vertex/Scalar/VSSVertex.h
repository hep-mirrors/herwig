// -*- C++ -*-
#ifndef HERWIG_VSSVertex_H
#define HERWIG_VSSVertex_H
//
// This is the declaration of the <!id>VSSVertex<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  The <!id>VSSVertex<!!id> class is the implementation of the vector-scalar-scalar
//  vertex. It inherits from the VertexBase class for storage of the particles and
//  implements the helicity calculations.
//
//  All such vertices should inherit from this class and implement the virtual
//  setCoupling member
//
// CLASSDOC SUBSECTION See also:
//
// <a href="VertexBase.html">VertexBase.h</a>.
// 

#include "Herwig++/Helicity/Vertex/VertexBase.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

class VSSVertex: public VertexBase {
    
public:
  
  inline VSSVertex();
  inline VSSVertex(const VSSVertex &);
  virtual ~VSSVertex();
  // Standard ctors and dtor.
  
public:
  
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.
  
  static void Init();
  // Standard Init function used to initialize the interfaces.
  
public:
  
  // evalaute the vertex
  Complex evaluate(Energy2,const VectorWaveFunction &,
		   const ScalarWaveFunction &, const ScalarWaveFunction &);
  VectorWaveFunction evaluate(Energy2,int,tcPDPtr, const ScalarWaveFunction &,
			      const ScalarWaveFunction &);
  ScalarWaveFunction evaluate(Energy2,int, tcPDPtr, const VectorWaveFunction &,
			      const ScalarWaveFunction &);
  // calculate the couplings
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
  
  static AbstractClassDescription<VSSVertex> initVSSVertex;
  // Describe a concrete class with persistent data.
  
  VSSVertex & operator=(const VSSVertex &);
  // Private and non-existent assignment operator.
  
};

}
}
#include "VSSVertex.icc"

// CLASSDOC OFF

namespace ThePEG {

  // The following template specialization informs ThePEG about the
  // base class of VSSVertex.
  template <>
  struct BaseClassTrait<Herwig::Helicity::VSSVertex,1> {
    typedef Herwig::Helicity::VertexBase NthBase;
  };
  
  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::Helicity::VSSVertex>
    : public ClassTraitsBase<Herwig::Helicity::VSSVertex> {
    static string className() { return "/Herwig++/Helicity/VSSVertex"; }
    // Return the class name.
    static string library() { return "libHwSVertex.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
    // (except the base class).
  };
  
}

#endif /* HERWIG_VSSVertex_H */

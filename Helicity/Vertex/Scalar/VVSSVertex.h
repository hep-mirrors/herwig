// -*- C++ -*-
#ifndef HERWIG_VVSSVertex_H
#define HERWIG_VVSSVertex_H
//
// This is the declaration of the <!id>VVSSVertex<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  The <!id>VVSSVertex<!!id> class is the implementation of the coupling of two
//  vectors and two scalars. It inherits from the VertexBase class for the 
//  storage of the particles and implements the helicity calculations.
//
//  All classes implementing the vertex should inherit from it and implement the
//  virtual setCoupling member.
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

class VVSSVertex: public VertexBase {
  
public:
  
  inline VVSSVertex();
  inline VVSSVertex(const VVSSVertex &);
  virtual ~VVSSVertex();
  // Standard ctors and dtor.
  
public:
  
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.
  
  static void Init();
  // Standard Init function used to initialize the interfaces.
  
public:
  
  // evalaute the vertex
  Complex evaluate(Energy2, const VectorWaveFunction &,
		   const VectorWaveFunction &, const ScalarWaveFunction &,
		   const ScalarWaveFunction &);
  VectorWaveFunction evaluate(Energy2, int,tcPDPtr, const VectorWaveFunction &,
			      const ScalarWaveFunction &,const ScalarWaveFunction &);
  ScalarWaveFunction evaluate(Energy2, int,tcPDPtr,const VectorWaveFunction &,
			      const VectorWaveFunction &,
			      const ScalarWaveFunction &);
  // calculate the couplings
  virtual void setCoupling(Energy2,tcPDPtr,tcPDPtr,tcPDPtr,tcPDPtr);
  
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
  
  static ClassDescription<VVSSVertex> initVVSSVertex;
  // Describe a concrete class with persistent data.
  
  VVSSVertex & operator=(const VVSSVertex &);
  // Private and non-existent assignment operator.
  
};

}
}
#include "VVSSVertex.icc"

// CLASSDOC OFF

namespace ThePEG {
  
  // The following template specialization informs ThePEG about the
  // base class of VVSSVertex.
  template <>
  struct BaseClassTrait<Herwig::Helicity::VVSSVertex,1> {
    typedef Herwig::Helicity::VertexBase NthBase;
  };
  
  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::Helicity::VVSSVertex>
    : public ClassTraitsBase<Herwig::Helicity::VVSSVertex> {
    static string className() { return "/Herwig++/Helicity/VVSSVertex"; }
    // Return the class name.
    static string library() { return "libHwSVertex.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
    // (except the base class).
  };

}


#endif /* HERWIG_VVSSVertex_H */

// -*- C++ -*-
#ifndef HERWIG_FFVTVertex_H
#define HERWIG_FFVTVertex_H
//
// This is the declaration of the <!id>FFVTVertex<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  The <!id>FFVTVertex<!!id> class is the implementation 
//  of the fermion-fermion--vector-tensor
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
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/TensorWaveFunction.h"

namespace Herwig {
using ThePEG::Helicity::DiracRep;
using ThePEG::Helicity::HELASDRep;
using ThePEG::Helicity::HaberDRep;
namespace Helicity {

using namespace ThePEG;

class FFVTVertex: public VertexBase {
      
public:
  
  inline FFVTVertex();
  inline FFVTVertex(const FFVTVertex &);
  virtual ~FFVTVertex();
  // Standard ctors and dtor.
  
public:
  
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.
  
  static void Init();
  // Standard Init function used to initialize the interfaces.
  
public:
  
  // evaluate the vertex
  Complex evaluate(Energy2,const SpinorWaveFunction &,
		   const SpinorBarWaveFunction &,
		   const VectorWaveFunction &, const TensorWaveFunction &);
  TensorWaveFunction evaluate(Energy2,int, tcPDPtr,const SpinorWaveFunction &,
			      const SpinorBarWaveFunction &,
			      const VectorWaveFunction &);
  VectorWaveFunction evaluate(Energy2,int, tcPDPtr, const SpinorWaveFunction &,
			      const SpinorBarWaveFunction &, 
			      const TensorWaveFunction &);
  SpinorWaveFunction evaluate(Energy2,int, tcPDPtr,const SpinorWaveFunction &,
			      const VectorWaveFunction &,
			      const TensorWaveFunction &);
  SpinorBarWaveFunction evaluate(Energy2,int, tcPDPtr,const SpinorBarWaveFunction &,
				 const VectorWaveFunction &,
				 const TensorWaveFunction &);
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
  
  static ClassDescription<FFVTVertex> initFFVTVertex;
  // Describe a concrete class with persistent data.
  
  FFVTVertex & operator=(const FFVTVertex &);
  // Private and non-existent assignment operator.
  
};
}
}
#include "FFVTVertex.icc"

// CLASSDOC OFF

namespace ThePEG {
  
  // The following template specialization informs ThePEG about the
  // base class of FFVTVertex.
  template <>
  struct BaseClassTrait<Herwig::Helicity::FFVTVertex,1> {
    typedef Herwig::Helicity::VertexBase NthBase;
  };
  
  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::Helicity::FFVTVertex>
    : public ClassTraitsBase<Herwig::Helicity::FFVTVertex> {
    static string className() { return "/Herwig++/Helicity/FFVTVertex"; }
    // Return the class name.
    static string library() { return "libHwTVertex.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
    // (except the base class).
  };
  
}


#endif /* HERWIG_FFVTVertex_H */

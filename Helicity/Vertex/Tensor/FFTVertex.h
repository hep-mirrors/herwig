// -*- C++ -*-
#ifndef HERWIG_FFTVertex_H
#define HERWIG_FFTVertex_H
//
// This is the declaration of the <!id>FFTVertex<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  The <!id>FFTVertex<!!id> class is the implementation of the fermion-fermion-tensor
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
#include "Herwig++/Helicity/WaveFunction/TensorWaveFunction.h"

namespace Herwig {
using ThePEG::Helicity::DiracRep;
using ThePEG::Helicity::HELASDRep;
using ThePEG::Helicity::HaberDRep;
using ThePEG::Helicity::defaultDRep;
namespace Helicity {

using namespace ThePEG;

class FFTVertex: public VertexBase {
      
public:
  
  inline FFTVertex();
  inline FFTVertex(const FFTVertex &);
  virtual ~FFTVertex();
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
		   const TensorWaveFunction &);
  // evaluate an off-shell tensor
  TensorWaveFunction evaluate(Energy2, int,tcPDPtr, const SpinorWaveFunction &,
			      const SpinorBarWaveFunction &);

  // evaluate an off-shell spinor (specify Dirac representation) 
  SpinorWaveFunction evaluate(Energy2, int,tcPDPtr, const SpinorWaveFunction &,
			      const TensorWaveFunction &,DiracRep=defaultDRep);

  // evaluate an off-shell spinorbar (specify Dirac representation)
  SpinorBarWaveFunction evaluate(Energy2,int, tcPDPtr,
				 const SpinorBarWaveFunction &,
				 const TensorWaveFunction&,DiracRep=defaultDRep);
  // set the couplings
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
  
  static AbstractClassDescription<FFTVertex> initFFTVertex;
  // Describe a concrete class with persistent data.
  
  FFTVertex & operator=(const FFTVertex &);
  // Private and non-existent assignment operator.
  
};

}
}

#include "FFTVertex.icc"

// CLASSDOC OFF

namespace ThePEG {

  // The following template specialization informs ThePEG about the
  // base class of FFTVertex.
  template <>
  struct BaseClassTrait<Herwig::Helicity::FFTVertex,1> {
    typedef Herwig::Helicity::VertexBase NthBase;
  };

  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::Helicity::FFTVertex>
    : public ClassTraitsBase<Herwig::Helicity::FFTVertex> {
    static string className() { return "/Herwig++/Helicity/FFTVertex"; }
    // Return the class name.
    static string library() { return "libHwTVertex.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
    // (except the base class).
  };
  
}

#endif /* HERWIG_FFTVertex_H */

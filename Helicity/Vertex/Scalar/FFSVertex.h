// -*- C++ -*-
#ifndef HERWIG_FFSVertex_H
#define HERWIG_FFSVertex_H
//
// This is the declaration of the <!id>FFSVertex<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  The <!id>FFSVertex<!!id> class is the implementation of the interact of a
//  scalar boson and a fermion-antifermion pair. It inherits from the VertexBase
//  class for storage of the particles interacting at the vertex and implements
//  the helicity calculations.
//
//  Implementations of specific interactions should inherit from this and implement
//  the virtual setCoupling member.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="VertexBase.html">VertexBase.h</a>.
// 
#include "Herwig++/Helicity/Vertex/VertexBase.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"

namespace Herwig {
using ThePEG::Helicity::DiracRep;
using ThePEG::Helicity::HELASDRep;
using ThePEG::Helicity::HaberDRep;
using ThePEG::Helicity::defaultDRep;

namespace Helicity{
using namespace ThePEG;

class FFSVertex: public VertexBase {
      
public:
  
  inline FFSVertex();
  inline FFSVertex(const FFSVertex &);
  virtual ~FFSVertex();
  // Standard ctors and dtor.
  
public:
  
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.
  
  static void Init();
  // Standard Init function used to initialize the interfaces.
  
public:      
  // set the left and right couplings
  inline void setLeft(Complex);
  inline void setRight(Complex);
  // get the left/right couplings
  inline Complex getLeft();
  inline Complex getRight();
  // evalulate the vertex
  Complex evaluate(Energy2,const SpinorWaveFunction &,
			   const SpinorBarWaveFunction &,
			   const ScalarWaveFunction &);

  // evaluate an off-shell spinor (specific Dirac representation)
  SpinorWaveFunction evaluate(Energy2,int,tcPDPtr,const SpinorWaveFunction &, 
			      const ScalarWaveFunction &,DiracRep=defaultDRep);

  // evalute an off-shell spinor bar (specific Dirac representation)
  SpinorBarWaveFunction evaluate(Energy2,int,tcPDPtr,const SpinorBarWaveFunction &,
				 const ScalarWaveFunction &,DiracRep=defaultDRep);

  // evaluate an off-shell scalar
  ScalarWaveFunction evaluate(Energy2,int,tcPDPtr,const SpinorWaveFunction &, 
			      const SpinorBarWaveFunction &);
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
  
  static AbstractClassDescription<FFSVertex> initFFSVertex;
  // Describe a concrete class with persistent data.
  
  FFSVertex & operator=(const FFSVertex &);
  // Private and non-existent assignment operator.
  
private:
  // storage of the left and right couplings
  Complex _left,_right;
};
}
}

#include "FFSVertex.icc"

// CLASSDOC OFF

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of FFSVertex.
template <>
struct BaseClassTrait<Herwig::Helicity::FFSVertex,1> {
  typedef Herwig::Helicity::VertexBase NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::Helicity::FFSVertex>
  : public ClassTraitsBase<Herwig::Helicity::FFSVertex> {
  static string className() { return "/Herwig++/Helicity/FFSVertex"; }
  // Return the class name.
  static string library() { return "libHwSVertex.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}


#endif /* HERWIG_FFSVertex_H */

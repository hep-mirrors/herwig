// -*- C++ -*-
#ifndef HERWIG_FFVVertex_H
#define HERWIG_FFVVertex_H
//
// This is the declaration of the <!id>FFVVertex<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  The <!id>FFVVertex<!!id> class is the base class for all helicity amplitude
//  vertices which used he renormalisable form for the fermion-fermion-vector
//  vertex. 
//
//  Any such vertices should inherit from this class and implement the virtual
//  setcoupling member function. The base VertexBase class is used to store the
//  particles allowed to interact at the vertex.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="VertexBase.html">VertexBase.h</a>.
// 

#include <Herwig++/Helicity/Vertex/VertexBase.h>
#include <Herwig++/Helicity/WaveFunction/VectorWaveFunction.h>
#include <Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h>
#include <Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h>

namespace Herwig {
 
using ThePEG::Helicity::DiracRep;
using ThePEG::Helicity::HELASDRep;
using ThePEG::Helicity::HaberDRep;
using ThePEG::Helicity::defaultDRep;

namespace Helicity{

using namespace ThePEG; 

class FFVVertex: public VertexBase {
  
public:
  
  inline FFVVertex();
  inline FFVVertex(const FFVVertex &);
  virtual ~FFVVertex();
  // Standard ctors and dtor.
  
public:
  
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.
  
  static void Init();
  // Standard Init function used to initialize the interfaces.
  
public:
  
  // set the left and right couplings
  inline void setLeft(const Complex &);
  inline void setRight(const Complex &);

  // get the left/right couplings
  inline const Complex & getLeft();
  inline const Complex & getRight();

  // evaluate the vertex
  Complex       evaluate(Energy2,const SpinorWaveFunction &,
				 const SpinorBarWaveFunction &, 
				 const VectorWaveFunction &);
  // evaulate an off-shell spinorbar 
  SpinorBarWaveFunction evaluate(Energy2,int,tcPDPtr,
				 const SpinorBarWaveFunction &,
				 const VectorWaveFunction &,DiracRep=defaultDRep);

  // evaluate an off-shell vector
  VectorWaveFunction    evaluate(Energy2,int,tcPDPtr,const SpinorWaveFunction &,
				 const SpinorBarWaveFunction &);

  // evaluate an off-shell spinor 
  SpinorWaveFunction    evaluate(Energy2,int,tcPDPtr,const SpinorWaveFunction &,
				 const VectorWaveFunction &,DiracRep=defaultDRep);

  // return the coupling
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
  
  static AbstractClassDescription<FFVVertex> initFFVVertex;
  // Describe a concrete class with persistent data.
  
  FFVVertex & operator=(const FFVVertex &);
  // Private and non-existent assignment operator.
  
private:
  // left and right couplings
  Complex _left,_right;
  
};
}
}

#include "FFVVertex.icc"

// CLASSDOC OFF

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of FFVVertex.
template <>
struct BaseClassTrait<Herwig::Helicity::FFVVertex,1> {
  typedef Herwig::Helicity::VertexBase NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::Helicity::FFVVertex>
  : public ClassTraitsBase<Herwig::Helicity::FFVVertex> {
  static string className() { return "/Herwig++/Helicity/FFVVertex"; }
  // Return the class name.
  static string library() { return "libHwVVertex.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}


#endif /* HERWIG_FFVVertex_H */

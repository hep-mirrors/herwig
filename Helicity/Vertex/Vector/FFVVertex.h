// -*- C++ -*-
#ifndef HERWIG_FFVVertex_H
#define HERWIG_FFVVertex_H
//
// This is the declaration of the FFVVertex class.

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

/** \ingroup Helicity
 *  The FFVVertex class is the base class for all helicity amplitude
 *  vertices which used he renormalisable form for the 
 *  fermion-fermion-vector vertex. 
 *
 *  Any such vertices should inherit from this class and implement the virtual
 *  setcoupling member function. The base VertexBase class is used to store the
 *  particles allowed to interact at the vertex.
 *
 *  @see VertexBase
 */
class FFVVertex: public VertexBase {
  
public:
  
  /**
   * Standard ctors and dtor.
   */
  inline FFVVertex();
  inline FFVVertex(const FFVVertex &);
  virtual ~FFVVertex();
  
public:
  
  /**
   * Standard functions for writing and reading from persistent streams.
   */
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  
  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();
  
public:
  
  /**
   * Set the left and right couplings.
   */
  inline void setLeft(const Complex &);
  inline void setRight(const Complex &);

  /**
   * Get the left/right couplings.
   */
  inline const Complex & getLeft();
  inline const Complex & getRight();

  /**
   * Evaluate the vertex.
   */
  Complex       evaluate(Energy2,const SpinorWaveFunction &,
				 const SpinorBarWaveFunction &, 
				 const VectorWaveFunction &);

  /**
   * Evaulate an off-shell spinorbar.
   */
  SpinorBarWaveFunction evaluate(Energy2,int,tcPDPtr,
				 const SpinorBarWaveFunction &,
				 const VectorWaveFunction &,DiracRep=defaultDRep);

  /**
   * Evaluate an off-shell vector.
   */
  VectorWaveFunction    evaluate(Energy2,int,tcPDPtr,const SpinorWaveFunction &,
				 const SpinorBarWaveFunction &);

  /**
   * Evaluate an off-shell spinor. 
   */
  SpinorWaveFunction    evaluate(Energy2,int,tcPDPtr,const SpinorWaveFunction &,
				 const VectorWaveFunction &,DiracRep=defaultDRep);

  /**
   * Return the coupling.
   */
  virtual void setCoupling(Energy2,tcPDPtr,tcPDPtr,tcPDPtr);
  
protected:
  
  /**
   * Standard Interfaced virtual functions.
   */
  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void doinitrun();
  inline virtual void dofinish();
  
  /**
   * Change all pointers to Interfaced objects to corresponding clones.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  
  /**
   * Return pointers to all Interfaced objects refered to by this.
   */
  inline virtual IVector getReferences();
  
private:
  
  /**
   * Describe a concrete class with persistent data.
   */
  static AbstractClassDescription<FFVVertex> initFFVVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  FFVVertex & operator=(const FFVVertex &);
  
private:

  /**
   * Left and right couplings.
   */
  Complex _left,_right;
  
};
}
}

#include "FFVVertex.icc"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of FFVVertex.
 */
template <>
struct BaseClassTrait<Herwig::Helicity::FFVVertex,1> {
  typedef Herwig::Helicity::VertexBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::Helicity::FFVVertex>
  : public ClassTraitsBase<Herwig::Helicity::FFVVertex> {

  /**
   * Return the class name.
   */
  static string className() { return "/Herwig++/Helicity/FFVVertex"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwVVertex.so"; }

};

}


#endif /* HERWIG_FFVVertex_H */

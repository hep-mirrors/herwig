// -*- C++ -*-
#ifndef HERWIG_SSSSVertex_H
#define HERWIG_SSSSVertex_H
//
// This is the declaration of the SSSSVertex class.

#include "Herwig++/Helicity/Vertex/VertexBase.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

/** \ingroup Helicity
 * 
 *  The SSSSVertex class is the implementation of the interaction of
 *  four scalars. It inherits from the VertexBase class for the storage 
 *  of the particles interacting at the vertex and implements the 
 *  helicity calculations.
 *
 *  Any classes implementating the vertex should inherit from it and implement
 *  the virtual set Coupling member.
 *
 *  @see VertexBase
 */
class SSSSVertex: public VertexBase {
  
public:
  
  /**
   * Standard ctors and dtor.
   */
  inline SSSSVertex();
  inline SSSSVertex(const SSSSVertex &);
  virtual ~SSSSVertex();
  
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
   * Evaluate the vertex.
   */
  Complex evaluate(Energy2, const ScalarWaveFunction &,
		   const ScalarWaveFunction &, const ScalarWaveFunction &,
		   const ScalarWaveFunction &);
  ScalarWaveFunction evaluate(Energy2,int, tcPDPtr, const ScalarWaveFunction &,
			      const ScalarWaveFunction &,const ScalarWaveFunction &);

  /**
   * Calculate the couplings.
   */
  virtual void setCoupling(Energy2,tcPDPtr,tcPDPtr,tcPDPtr,tcPDPtr);
  
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
  static AbstractClassDescription<SSSSVertex> initSSSSVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  SSSSVertex & operator=(const SSSSVertex &);
  
};
}
}

#include "SSSSVertex.icc"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of SSSSVertex.
 */
template <>
struct BaseClassTrait<Herwig::Helicity::SSSSVertex,1> {
  typedef Herwig::Helicity::VertexBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::Helicity::SSSSVertex>
  : public ClassTraitsBase<Herwig::Helicity::SSSSVertex> {

  /**
   * Return the class name.
   */
  static string className() { return "/Herwig++/Helicity/SSSSVertex"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwSVertex.so"; }

};

}


#endif /* HERWIG_SSSSVertex_H */

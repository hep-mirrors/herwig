// -*- C++ -*-
#ifndef HERWIG_VVSVertex_H
#define HERWIG_VVSVertex_H
//
// This is the declaration of the VVSVertex class.

#include "Herwig++/Helicity/Vertex/VertexBase.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

/** \ingroup Helicity
 *
 * The VVSVertex class is the implementation of the vector-vector-scalar.
 * It inherits from the VertexBase class for the storage of the particles and
 * implements the helicity calculations.
 *
 * All interactions of this type should inherit from it and implement the virtual
 * setCoupling member
 *
 * @see VertexBase
 */
class VVSVertex: public VertexBase {
  
public:
  
  /**
   * Standard ctors and dtor.
   */
  inline VVSVertex();
  inline VVSVertex(const VVSVertex &);
  virtual ~VVSVertex();
  
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
   * Evalaute the vertex.
   */
  Complex evaluate(Energy2,const VectorWaveFunction &,
		   const VectorWaveFunction &, const ScalarWaveFunction &);
  VectorWaveFunction evaluate(Energy2,int,tcPDPtr, const VectorWaveFunction &,
			      const ScalarWaveFunction &);
  ScalarWaveFunction evaluate(Energy2,int, tcPDPtr, const VectorWaveFunction &,
			      const VectorWaveFunction &);

  /**
   * Calculate the couplings.
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
  static AbstractClassDescription<VVSVertex> initVVSVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  VVSVertex & operator=(const VVSVertex &);
  
};

}
}
#include "VVSVertex.icc"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of VVSVertex.
 */
template <>
struct BaseClassTrait<Herwig::Helicity::VVSVertex,1> {
  typedef Herwig::Helicity::VertexBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::Helicity::VVSVertex>
  : public ClassTraitsBase<Herwig::Helicity::VVSVertex> {

  /**
   * Return the class name.
   */
  static string className() { return "/Herwig++/Helicity/VVSVertex"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwSVertex.so"; }

};

}


#endif /* HERWIG_VVSVertex_H */

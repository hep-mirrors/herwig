// -*- C++ -*-
#ifndef HERWIG_VVSSVertex_H
#define HERWIG_VVSSVertex_H
//
// This is the declaration of the VVSSVertex class.

#include "Herwig++/Helicity/Vertex/VertexBase.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"

namespace Herwig {
namespace Helicity {

/** \ingroup Helicity
 *
 *  The VVSSVertex class is the implementation of the coupling of two
 *  vectors and two scalars. It inherits from the VertexBase class for the 
 *  storage of the particles and implements the helicity calculations.
 *
 *  All classes implementing the vertex should inherit from it and implement the
 *  virtual setCoupling member.
 *
 *  @see VertexBase
 */
class VVSSVertex: public VertexBase {
  
public:
  
  /**
   * Standard ctors and dtor.
   */
  inline VVSSVertex();
  inline VVSSVertex(const VVSSVertex &);
  virtual ~VVSSVertex();
  
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
  Complex evaluate(Energy2, const VectorWaveFunction &,
		   const VectorWaveFunction &, const ScalarWaveFunction &,
		   const ScalarWaveFunction &);
  VectorWaveFunction evaluate(Energy2, int,tcPDPtr, const VectorWaveFunction &,
			      const ScalarWaveFunction &,const ScalarWaveFunction &);
  ScalarWaveFunction evaluate(Energy2, int,tcPDPtr,const VectorWaveFunction &,
			      const VectorWaveFunction &,
			      const ScalarWaveFunction &);

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
  static AbstractClassDescription<VVSSVertex> initVVSSVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  VVSSVertex & operator=(const VVSSVertex &);
  
};

}
}

#include "VVSSVertex.icc"

namespace ThePEG {
  
  /**
   * The following template specialization informs ThePEG about the
   * base class of VVSSVertex.
   */
  template <>
  struct BaseClassTrait<Herwig::Helicity::VVSSVertex,1> {
    typedef Herwig::Helicity::VertexBase NthBase;
  };
  
  /**
   * The following template specialization informs ThePEG about the
   * name of this class and the shared object where it is defined.
   */
  template <>
  struct ClassTraits<Herwig::Helicity::VVSSVertex>
    : public ClassTraitsBase<Herwig::Helicity::VVSSVertex> {

    /**
     * Return the class name.
     */
    static string className() { return "/Herwig++/Helicity/VVSSVertex"; }

    /**
     * Return the name of the shared library to be loaded to get
     * access to this class and every other class it uses
     * (except the base class).
     */
    static string library() { return "libHwSVertex.so"; }

  };

}


#endif /* HERWIG_VVSSVertex_H */

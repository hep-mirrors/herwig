// -*- C++ -*-
#ifndef HERWIG_SSTVertex_H
#define HERWIG_SSTVertex_H
//
// This is the declaration of the SSTVertex class.

#include "Herwig++/Helicity/Vertex/VertexBase.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/TensorWaveFunction.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

/** \ingroup Helicity
 *  The VVTVertexclass is the implementation of the 
 *  scalar-scalar-tensor vertex. 
 *  It inherits from the VertexBase class for the storage of the particles
 *  interacting at the vertex and implements the helicity amplitude calculations.
 *
 *  All implementations of this vertex should inherit from it and implement the
 *  virtual setCoupling member.
 *
 *  @see VertexBase
 */
class SSTVertex: public VertexBase {
  
public:
  
  /**
   * Standard ctors and dtor.
   */
  inline SSTVertex();
  inline SSTVertex(const SSTVertex &);
  virtual ~SSTVertex();
  
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
  
  Complex evaluate(Energy2, const ScalarWaveFunction &,
    		       const ScalarWaveFunction &, const TensorWaveFunction &);
  ScalarWaveFunction evaluate(Energy2,int,tcPDPtr, const ScalarWaveFunction &,
    			  const TensorWaveFunction &);
  TensorWaveFunction evaluate(Energy2,int,tcPDPtr, const ScalarWaveFunction &,
    			  const ScalarWaveFunction &);
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
  static AbstractClassDescription<SSTVertex> initSSTVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  SSTVertex & operator=(const SSTVertex &);
  
};

}
}

#include "SSTVertex.icc"

namespace ThePEG {
  
  /**
   * The following template specialization informs ThePEG about the
   * base class of SSTVertex.
   */
  template <>
  struct BaseClassTrait<Herwig::Helicity::SSTVertex,1> {
    typedef Herwig::Helicity::VertexBase NthBase;
  };
  
  /**
   * The following template specialization informs ThePEG about the
   * name of this class and the shared object where it is defined.
   */
  template <>
  struct ClassTraits<Herwig::Helicity::SSTVertex>
    : public ClassTraitsBase<Herwig::Helicity::SSTVertex> {

    /**
     * Return the class name.
     */
    static string className() { return "/Herwig++/Helicity/SSTVertex"; }

    /**
     * Return the name of the shared library to be loaded to get
     * access to this class and every other class it uses
     * (except the base class).
     */
    static string library() { return "libHwTVertex.so"; }

  };
  
}


#endif /* HERWIG_SSTVertex_H */

// -*- C++ -*-
#ifndef HERWIG_VVTVertex_H
#define HERWIG_VVTVertex_H
//
// This is the declaration of the VVTVertex class.

#include "Herwig++/Helicity/Vertex/VertexBase.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/TensorWaveFunction.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

/** \ingroup Helicity
 *
 *  The VVTVertex class is the implementation of the 
 *  vector-vector-tensor vertex. 
 *  It inherits from the VertexBase class for the storage of the particles
 *  interacting at the vertex and implements the helicity amplitude calculations.
 *
 *  All implementations of this vertex should inherit from it and implement the
 *  virtual setCoupling member.
 *
 *  @see VertexBase
 */
class VVTVertex: public VertexBase {
  
public:
  
  /**
   * Standard ctors and dtor.
   */
  inline VVTVertex();
  inline VVTVertex(const VVTVertex &);
  virtual ~VVTVertex();
  
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
  Complex evaluate(Energy2,const VectorWaveFunction &,
    		       const VectorWaveFunction &, const TensorWaveFunction &);
  TensorWaveFunction evaluate(Energy2,int, tcPDPtr, const VectorWaveFunction &,
    			  const VectorWaveFunction &);
  VectorWaveFunction evaluate(Energy2,int, tcPDPtr, const VectorWaveFunction &,
    			  const TensorWaveFunction &);
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
  static AbstractClassDescription<VVTVertex> initVVTVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  VVTVertex & operator=(const VVTVertex &);
  
};

}
}

#include "VVTVertex.icc"

namespace ThePEG {
  
  /**
   * The following template specialization informs ThePEG about the
   * base class of VVTVertex.
   */
  template <>
  struct BaseClassTrait<Herwig::Helicity::VVTVertex,1> {
    typedef Herwig::Helicity::VertexBase NthBase;
  };
  
  /**
   * The following template specialization informs ThePEG about the
   * name of this class and the shared object where it is defined.
   */
  template <>
  struct ClassTraits<Herwig::Helicity::VVTVertex>
    : public ClassTraitsBase<Herwig::Helicity::VVTVertex> {

    /**
     * Return the class name.
     */
    static string className() { return "/Herwig++/Helicity/VVTVertex"; }

    /**
     * Return the name of the shared library to be loaded to get
     * access to this class and every other class it uses
     * (except the base class).
     */
    static string library() { return "libHwTVertex.so"; }

  };
  
}


#endif /* HERWIG_VVTVertex_H */

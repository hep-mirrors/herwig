// -*- C++ -*-
#ifndef HERWIG_VVVTVertex_H
#define HERWIG_VVVTVertex_H
//
// This is the declaration of the VVVTVertex class.

#include "Herwig++/Helicity/Vertex/VertexBase.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/TensorWaveFunction.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;
  
/** \ingroup Helicity
 *
 *  The VVTVertex class is the implementation of the 
 *  vector-vector-vector-tensor vertex. 
 *  It inherits from the VertexBase class for the storage of the particles
 *  interacting at the vertex and implements the helicity amplitude calculations.
 *
 *  All implementations of this vertex should inherit from it and implement the
 *  virtual setCoupling member.
 *
 *  @see VertexBase
 */
class VVVTVertex: public VertexBase {
    
public:
  
  /**
   * Standard ctors and dtor.
   */
  inline VVVTVertex();
  inline VVVTVertex(const VVVTVertex &);
  virtual ~VVVTVertex();
  
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
		   const VectorWaveFunction &,
		   const VectorWaveFunction &, const TensorWaveFunction &);
  TensorWaveFunction evaluate(Energy2,int, tcPDPtr,const VectorWaveFunction &,
			      const VectorWaveFunction &,const VectorWaveFunction &);
  VectorWaveFunction evaluate(Energy2,int, tcPDPtr, const VectorWaveFunction &,
			      const VectorWaveFunction &,
			      const TensorWaveFunction &);
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
  static AbstractClassDescription<VVVTVertex> initVVVTVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  VVVTVertex & operator=(const VVVTVertex &);
  
};

}
}

#include "VVVTVertex.icc"

namespace ThePEG {

/** 
 * The following template specialization informs ThePEG about the
 * base class of VVVTVertex.
 */
template <>
struct BaseClassTrait<Herwig::Helicity::VVVTVertex,1> {
typedef Herwig::Helicity::VertexBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::Helicity::VVVTVertex>
  : public ClassTraitsBase<Herwig::Helicity::VVVTVertex> {

  /**
   * Return the class name.
   */
  static string className() { return "/Herwig++/Helicity/VVVTVertex"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwTVertex.so"; }

};

}

#endif /* HERWIG_VVVTVertex_H */

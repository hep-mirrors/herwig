// -*- C++ -*-
#ifndef HERWIG_VVVVertex_H
#define HERWIG_VVVVertex_H
//
// This is the declaration of the VVVVertex class.

#include "Herwig++/Helicity/Vertex/VertexBase.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"

namespace Herwig {
namespace Helicity{
using namespace ThePEG;
  
/** \ingroup Helicity
 *
 *  The VVVVertex class is the base class for tripe vectro vertices
 *  using the perturbative form in Herwig++. 
 *  It inherits from the VertexBase class for the storage of the 
 *  particles allowed at the vertex.
 *
 *  Classes which implement a specific vertex should inherit from this and
 *  implement the virtual setCoupling member.
 *
 *  @see VertexBase
 */
class VVVVertex: public VertexBase {
    
public:
  
  /**
   * Standard ctors and dtor.
   */
  inline VVVVertex();
  inline VVVVertex(const VVVVertex &);
  virtual ~VVVVertex();
  
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
  Complex evaluate(Energy2, const VectorWaveFunction &,
		   const VectorWaveFunction &, const VectorWaveFunction &);
  VectorWaveFunction evaluate(Energy2,int, tcPDPtr, const VectorWaveFunction &,
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
  static AbstractClassDescription<VVVVertex> initVVVVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  VVVVertex & operator=(const VVVVertex &);
  
};

}
}

#include "VVVVertex.icc"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of VVVVertex.
 */
template <>
struct BaseClassTrait<Herwig::Helicity::VVVVertex,1> {
  typedef Herwig::Helicity::VertexBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::Helicity::VVVVertex>
  : public ClassTraitsBase<Herwig::Helicity::VVVVertex> {

  /**
   * Return the class name.
   */
  static string className() { return "/Herwig++/Helicity/VVVVertex"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwVVertex.so"; }

};

}


#endif /* HERWIG_VVVVertex_H */

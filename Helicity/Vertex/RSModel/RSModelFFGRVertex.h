// -*- C++ -*-
#ifndef HERWIG_RSModelFFGRVertex_H
#define HERWIG_RSModelFFGRVertex_H
//
// This is the declaration of the RSModelFFGRVertex class.

#include "Herwig++/Helicity/Vertex/Tensor/FFTVertex.h"
#include "Herwig++/Models/RSModel/RSModel.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

/** \ingroup Helicity
 *
 *  This is the implementation of the Randell-Sundrum model fermion-antifermion
 *  tensor vertex for helicity amplitude calculations 
 *
 *  @see FFTVertex
 *  @see VertexBase
 */
class RSModelFFGRVertex: public FFTVertex {
  
public:
  
  /**
   * Standard ctors and dtor.
   */
  inline RSModelFFGRVertex();
  inline RSModelFFGRVertex(const RSModelFFGRVertex &);
  virtual ~RSModelFFGRVertex();
  
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

  /**
   * Set the couplings.
   */
  void setCoupling(Energy2,tcPDPtr,tcPDPtr, tcPDPtr);

protected:
  
  /**
   * Standard clone methods.
   */
  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;
  
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
  static ClassDescription<RSModelFFGRVertex> initRSModelFFGRVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  RSModelFFGRVertex & operator=(const RSModelFFGRVertex &);
  
  /**
   * Pinter to the model object.
   */
  SMPtr _theModel;

  /**
   * The coupling.
   */
  double _theKappa;

};
}
}

#include "RSModelFFGRVertex.icc"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of RSModelFFGRVertex.
 */
template <>
struct BaseClassTrait<Herwig::Helicity::RSModelFFGRVertex,1> {
  typedef Herwig::Helicity::FFTVertex NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::Helicity::RSModelFFGRVertex>
  : public ClassTraitsBase<Herwig::Helicity::RSModelFFGRVertex> {

  /**
   * Return the class name.
   */
  static string className() { return "/Herwig++/Helicity/RSModelFFGRVertex"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwRSVertex.so"; }

};

}
#endif /* HERWIG_RSModelFFGRVertex_H */

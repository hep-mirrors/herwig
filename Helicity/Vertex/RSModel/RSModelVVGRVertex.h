// -*- C++ -*-
#ifndef HERWIG_RSModelVVGRVertex_H
#define HERWIG_RSModelVVGRVertex_H
//
// This is the declaration of the RSModelVVGRVertex class.

#include "Herwig++/Helicity/Vertex/Tensor/VVTVertex.h"
#include "Herwig++/Models/RSModel/RSModel.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

/** \ingroup Helicity
 * 
 *  This is the implementation of the vector-vector-graviton vertex for
 *  the RS model
 * 
 *  @see VVTVertex
 *  @see VertexBase
 */
class RSModelVVGRVertex: public VVTVertex {
  
public:
  
  /**
   * Standard ctors and dtor.
   */
  inline RSModelVVGRVertex();
  inline RSModelVVGRVertex(const RSModelVVGRVertex &);
  virtual ~RSModelVVGRVertex();
  
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
   * Calculate the couplings.
   */
  void setCoupling(Energy2,tcPDPtr,tcPDPtr,tcPDPtr);

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
  static ClassDescription<RSModelVVGRVertex> initRSModelVVGRVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  RSModelVVGRVertex & operator=(const RSModelVVGRVertex &);

  /**
   * Pointer to the model object.
   */
  SMPtr _theModel;

  /**
   * The coupling.
   */
  double _theKappa;
  
};

}
}

#include "RSModelVVGRVertex.icc"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of RSModelVVGRVertex.
 */
template <>
struct BaseClassTrait<Herwig::Helicity::RSModelVVGRVertex,1> {
  typedef Herwig::Helicity::VVTVertex NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::Helicity::RSModelVVGRVertex>
  : public ClassTraitsBase<Herwig::Helicity::RSModelVVGRVertex> {

  /**
   * Return the class name.
   */
  static string className() { return "/Herwig++/Helicity/RSModelVVGRVertex"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwRSVertex.so"; }

};

}


#endif /* HERWIG_RSModelVVGRVertex_H */

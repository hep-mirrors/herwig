// -*- C++ -*-
#ifndef HERWIG_RSModelFFVGRVertex_H
#define HERWIG_RSModelFFVGRVertex_H
//
// This is the declaration of the RSModelFFVGRVertex class.

#include "Herwig++/Helicity/Vertex/Tensor/FFVTVertex.h"
#include "Herwig++/Models/RSModel/RSModel.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

/** \ingroup Helicity
 *
 *  This is the implementation of the fermion-antifermion-vector-graviton
 *  vertex for the Randell-Sundrum model
 *
 *  @see FFVTVertex
 *  @see VertexBase
 */
class RSModelFFVGRVertex: public FFVTVertex {
  
public:
  
  /**
   * Standard ctors and dtor.
   */
  inline RSModelFFVGRVertex();
  inline RSModelFFVGRVertex(const RSModelFFVGRVertex &);
  virtual ~RSModelFFVGRVertex();
  
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
  void setCoupling(Energy2,tcPDPtr,tcPDPtr,tcPDPtr,tcPDPtr);

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
  static ClassDescription<RSModelFFVGRVertex> initRSModelFFVGRVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  RSModelFFVGRVertex & operator=(const RSModelFFVGRVertex &);

private:

  /**
   * Pointer to the Standard Model object.
   */
  SMPtr _theModel;

  /**
   * Storage of the couplings.
   */
  double _charge[17];
  Complex _couplast[2];
  Energy2 _q2last[2];
  double _theKappa;

};

}
}

#include "RSModelFFVGRVertex.icc"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of RSModelFFVGRVertex.
 */
template <>
struct BaseClassTrait<Herwig::Helicity::RSModelFFVGRVertex,1> {
  typedef Herwig::Helicity::FFVTVertex NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::Helicity::RSModelFFVGRVertex>
  : public ClassTraitsBase<Herwig::Helicity::RSModelFFVGRVertex> {

  /**
   * Return the class name.
   */
  static string className() { return "/Herwig++/Helicity/RSModelFFVGRVertex"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwRSVertex.so"; }

};

}


#endif /* HERWIG_RSModelFFVGRVertex_H */

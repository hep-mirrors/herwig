// -*- C++ -*-
#ifndef HERWIG_RSModelVVVGRVertex_H
#define HERWIG_RSModelVVVGRVertex_H
//
// This is the declaration of the RSModelVVVGRVertex class.

#include "Herwig++/Helicity/Vertex/Tensor/VVVTVertex.h"
#include "Herwig++/Models/RSModel/RSModel.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;
    
/** \ingroup Helicity
 *
 *  The RSModelVVVGRVertex class is the implementation of the 
 *  triple vector graviton couling in the RS model. 
 *  It inherits from VVVTVertex and implements the setCoupling member.
 *
 *  @see VVVTVertex
 *  @see VertexBase
 */
class RSModelVVVGRVertex: public VVVTVertex {
  
public:
  
  /**
   * Standard ctors and dtor.
   */
  inline RSModelVVVGRVertex();
  inline RSModelVVVGRVertex(const RSModelVVVGRVertex &);
  virtual ~RSModelVVVGRVertex();
  
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
   * Set the coupling.
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
  static ClassDescription<RSModelVVVGRVertex> initRSModelVVVGRVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  RSModelVVVGRVertex & operator=(const RSModelVVVGRVertex &);

private:

  /**
   * Pointer to the model object.
   */
  SMPtr _theModel;

  /**
   * The graviton coupling.
   */
  double _theKappa;

  /**
   * Storage of the couplings.
   */
  Complex _couplast[2];
  Energy2 _q2last[2];
  double _zfact;

};

}
}

#include "RSModelVVVGRVertex.icc"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of RSModelVVVGRVertex. 
 */
template <>
struct BaseClassTrait<Herwig::Helicity::RSModelVVVGRVertex,1> {
  typedef Herwig::Helicity::VVVTVertex NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::Helicity::RSModelVVVGRVertex>
  : public ClassTraitsBase<Herwig::Helicity::RSModelVVVGRVertex> {

  /**
   * Return the class name.
   */
  static string className() { return "/Herwig++/Helicity/RSModelVVVGRVertex"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwRSVertex.so"; }

};

}


#endif /* HERWIG_RSModelVVVGRVertex_H */

// -*- C++ -*-
#ifndef HERWIG_SMGGGGVertex_H
#define HERWIG_SMGGGGVertex_H
//
// This is the declaration of the SMGGGGVertex class.

#include "Herwig++/Helicity/Vertex/Vector/VVVVVertex.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Utilities/Rebinder.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;
    
/** \ingroup Helicity
 *
 *  The SMGGGGVertex class is the implementation of the 
 *  Standard Model quartic gluon vertex. It inherits from 
 *  VVVVVertex and implements the setCoupling member.
 *
 *  @see VVVVVertex
 *  @see VertexBase
 */
class SMGGGGVertex: public VVVVVertex {
  
public:
  
  /**
   * Standard ctors and dtor.
   */
  inline SMGGGGVertex();
  inline SMGGGGVertex(const SMGGGGVertex &);
  virtual ~SMGGGGVertex();
  
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
   * Calculate the coupling.
   */
  void setCoupling(Energy2,tcPDPtr, tcPDPtr, tcPDPtr,tcPDPtr);

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
  static ClassDescription<SMGGGGVertex> initSMGGGGVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  SMGGGGVertex & operator=(const SMGGGGVertex &);
  
  /**
   * Pointer to the Standard Model.
   */
  SMPtr _theSM;

  /**
   * Storage of the couplings.
   */
  Complex _couplast;
  Energy2 _q2last;

};

}
}

#include "SMGGGGVertex.icc"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of SMGGGGVertex.
 */
template <>
struct BaseClassTrait<Herwig::Helicity::SMGGGGVertex,1> {
  typedef Herwig::Helicity::VVVVVertex NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::Helicity::SMGGGGVertex>
  : public ClassTraitsBase<Herwig::Helicity::SMGGGGVertex> {

  /**
   * Return the class name.
   */
  static string className() { return "/Herwig++/Helicity/SMGGGGVertex"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwSMVertex.so"; }

};

}


#endif /* HERWIG_SMGGGGVertex_H */

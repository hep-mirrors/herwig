// -*- C++ -*-
#ifndef HERWIG_SMGGGVertex_H
#define HERWIG_SMGGGVertex_H
//
// This is the declaration of the SMGGGVertex class.

#include "Herwig++/Helicity/Vertex/Vector/VVVVertex.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Utilities/Rebinder.h"


namespace Herwig {
namespace Helicity{
using namespace ThePEG;

/** \ingroup Helicity
 *
 *  Implementation of the SM triple gluon vertex.
 *
 */
class SMGGGVertex: public VVVVertex {

public:

  /**
   * Standard ctors and dtor.
   */
  inline SMGGGVertex();
  inline SMGGGVertex(const SMGGGVertex &);
  virtual ~SMGGGVertex();
  
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
  void setCoupling(Energy2,tcPDPtr, tcPDPtr, tcPDPtr);

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
  static ClassDescription<SMGGGVertex> initSMGGGVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  SMGGGVertex & operator=(const SMGGGVertex &);

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

#include "SMGGGVertex.icc"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of SMGGGVertex.
 */
template <>
struct BaseClassTrait<Herwig::Helicity::SMGGGVertex,1> {
  typedef Herwig::Helicity::VVVVertex NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::Helicity::SMGGGVertex>
  : public ClassTraitsBase<Herwig::Helicity::SMGGGVertex> {

  /**
   * Return the class name.
   */
  static string className() { return "/Herwig++/Helicity/SMGGGVertex"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwSMVertex.so"; }

};

}

#endif /* HERWIG_SMGGGVertex_H */

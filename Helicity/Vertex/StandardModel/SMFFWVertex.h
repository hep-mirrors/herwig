// -*- C++ -*-
#ifndef HERWIG_SMFFWVertex_H
#define HERWIG_SMFFWVertex_H
//
// This is the declaration of the SMFFWVertex class.

#include "Herwig++/Helicity/Vertex/Vector/FFVVertex.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/StandardModel/CKMBase.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "Herwig++/Models/StandardModel/StandardCKM.h"
namespace Herwig {
namespace Helicity {
using namespace ThePEG;

/** \ingroup Helicity
 *
 *  This is the implementation of the Standard model coupling 
 *  of the W to the fermions.
 *
 *  @see FFVVertex
 *  @see VertexBase
 */
class SMFFWVertex: public FFVVertex {
  
public:
  
  /**
   * Standard ctors and dtor.
   */
  inline SMFFWVertex();
  inline SMFFWVertex(const SMFFWVertex &);
  virtual ~SMFFWVertex();
  
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
  static ClassDescription<SMFFWVertex> initSMFFWVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  SMFFWVertex & operator=(const SMFFWVertex &);

private:

  /**
   * Pointer to the Standard Model object.
   */
  SMPtr _theSM;
  Ptr<CKMBase>::pointer _theCKM;

  /**
   * Storage of the couplings.
   */
  Complex _ckm[3][3];
  Complex _couplast;
  Energy2 _q2last;
  
};

}
}
#include "SMFFWVertex.icc"

namespace ThePEG {
  
  /**
   * The following template specialization informs ThePEG about the
   * base class of SMFFWVertex.
   */
  template <>
  struct BaseClassTrait<Herwig::Helicity::SMFFWVertex,1> {
    typedef Herwig::Helicity::FFVVertex NthBase;
  };
  
  /**
   * The following template specialization informs ThePEG about the
   * name of this class and the shared object where it is defined.
   */
  template <>
  struct ClassTraits<Herwig::Helicity::SMFFWVertex>
    : public ClassTraitsBase<Herwig::Helicity::SMFFWVertex> {

    /**
     * Return the class name.
     */
    static string className() { return "/Herwig++/Helicity/SMFFWVertex"; }

    /**
     * Return the name of the shared library to be loaded to get
     * access to this class and every other class it uses
     * (except the base class).
     */
    static string library() { return "libHwSMVertex.so"; }

  };
  
}

#endif /* HERWIG_SMFFWVertex_H */

// -*- C++ -*-
#ifndef HERWIG_SMWWWWVertex_H
#define HERWIG_SMWWWWVertex_H
//
// This is the declaration of the SMWWWWVertex class.

#include "Herwig++/Helicity/Vertex/Vector/VVVVVertex.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Utilities/Rebinder.h"


namespace Herwig {
namespace Helicity {
using namespace ThePEG;
    
/** \ingroup Helicity
 *
 *  The SMWWWWVertex class is the implementation of the
 *  coupling of four electroweak gauge bosons in the SM. 
 *  It inherits from VVVVVertex nad implements the setCoupling member.
 *
 *  @see VVVVVVertex
 *  @see VertexBase
 */
class SMWWWWVertex: public VVVVVertex {
  
public:
  
  /**
   * Standard ctors and dtor. 
   */
  inline SMWWWWVertex();
  inline SMWWWWVertex(const SMWWWWVertex &);
  virtual ~SMWWWWVertex();
  
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
   * Calculate the coupling
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
  static ClassDescription<SMWWWWVertex> initSMWWWWVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  SMWWWWVertex & operator=(const SMWWWWVertex &);

  /**
   * Pointer to the SM object.
   */
  SMPtr _theSM;

  /**
   * Intermediate particles.
   */
  tcPDPtr _gamma,_Z0,_wplus,_wminus;

  /**
   * Storage of the couplings.
   */
  Complex _couplast;
  Energy2 _q2last;
  double _vfact[4],_sw2,_cw2;

};
 
}
}

#include "SMWWWWVertex.icc"

namespace ThePEG {
  
  /**
   * The following template specialization informs ThePEG about the
   * base class of SMWWWWVertex.
   */
  template <>
  struct BaseClassTrait<Herwig::Helicity::SMWWWWVertex,1> {
    typedef Herwig::Helicity::VVVVVertex NthBase;
  };
  
  /**
   * The following template specialization informs ThePEG about the
   * name of this class and the shared object where it is defined.
   */
  template <>
  struct ClassTraits<Herwig::Helicity::SMWWWWVertex>
    : public ClassTraitsBase<Herwig::Helicity::SMWWWWVertex> {

    /**
     * Return the class name.
     */
    static string className() { return "/Herwig++/Helicity/SMWWWWVertex"; }

    /**
     * Return the name of the shared library to be loaded to get
     * access to this class and every other class it uses
     * (except the base class).
     */ 
    static string library() { return "libHwSMVertex.so"; }

  };
  
}


#endif /* HERWIG_SMWWWWVertex_H */

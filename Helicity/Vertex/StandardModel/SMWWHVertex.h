// -*- C++ -*-
#ifndef HERWIG_SMWWHVertex_H
#define HERWIG_SMWWHVertex_H
//
// This is the declaration of the SMWWHVertex class.

#include "Herwig++/Helicity/Vertex/Scalar/VVSVertex.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/Utilities/Rebinder.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

/** \ingroup Helicity
 *
 *  The SMWWHVertex is the implementation of the
 *  coupling of two electroweak gauge bosons to the Higgs in the Standard
 *  Model. It inherits from VVSVertex and implements the setCoupling member.
 *
 *  @see VVSVertex
 *  @see VertexBase
 */
class SMWWHVertex: public VVSVertex {
  
public:
  
  /**
   * Standard ctors and dtor.
   */
  inline SMWWHVertex();
  inline SMWWHVertex(const SMWWHVertex &);
  virtual ~SMWWHVertex();
  
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
  static ClassDescription<SMWWHVertex> initSMWWHVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  SMWWHVertex & operator=(const SMWWHVertex &);
  
  /**
   * Pointer to he Standard Model object.
   */
  SMPtr _theSM;

  /**
   * Storage of the couplings.
   */
  Complex _couplast;
  Energy2 _q2last;
  Energy _mw;
  double _zfact,_sw;
  
};

}
}

#include "SMWWHVertex.icc"

namespace ThePEG {

  /**
   * The following template specialization informs ThePEG about the
   * base class of SMWWHVertex.
   */
  template <>
  struct BaseClassTrait<Herwig::Helicity::SMWWHVertex,1> {
    typedef Herwig::Helicity::VVSVertex NthBase;
  };
  
  /**
   * The following template specialization informs ThePEG about the
   * name of this class and the shared object where it is defined.
   */
  template <>
  struct ClassTraits<Herwig::Helicity::SMWWHVertex>
    : public ClassTraitsBase<Herwig::Helicity::SMWWHVertex> {

    /**
     * Return the class name.
     */
    static string className() { return "/Herwig++/Helicity/SMWWHVertex"; }

    /** 
     * Return the name of the shared library to be loaded to get
     * access to this class and every other class it uses
     * (except the base class).
     */
    static string library() { return "libHwSMVertex.so"; }

  };
  
}


#endif /* HERWIG_SMWWHVertex_H */

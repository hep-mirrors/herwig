// -*- C++ -*-
#ifndef HERWIG_SMFFHVertex_H
#define HERWIG_SMFFHVertex_H
//
// This is the declaration of the SMFFHVertex class.

#include "Herwig++/Helicity/Vertex/Scalar/FFSVertex.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Utilities/Rebinder.h"


namespace Herwig {
namespace Helicity {
using namespace ThePEG;

/** \ingroup Helicity
 *
 *  This is the implementation of the vertex coupling the Standard Model Higgs
 *  to the Standard Model fermions for helicity amplitude calculations
 *
 *  @see FFSVertex
 *  @see VertexBase
 */
class SMFFHVertex: public FFSVertex {
  
public:
  
  /**
   * Standard ctors and dtor.
   */
  inline SMFFHVertex();
  inline SMFFHVertex(const SMFFHVertex &);
  virtual ~SMFFHVertex();
  
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
  static ClassDescription<SMFFHVertex> initSMFFHVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  SMFFHVertex & operator=(const SMFFHVertex &);

  /**
   * Pointer to the SM object.
   */
  Ptr<Herwig::StandardModel>::pointer _theSM;

  /**
   * Storage of the couplings.
   */
  Complex _couplast;
  double _sw;
  int _idlast;
  Energy2 _q2last;
  Energy _masslast,_mw;

};  

}
} 
#include "SMFFHVertex.icc"

namespace ThePEG {
  
  /**
   * The following template specialization informs ThePEG about the
   * base class of SMFFHVertex.
   */
  template <>
  struct BaseClassTrait<Herwig::Helicity::SMFFHVertex,1> {
    typedef Herwig::Helicity::FFSVertex NthBase;
  };
  
  /**
   * The following template specialization informs ThePEG about the
   * name of this class and the shared object where it is defined.
   */
  template <>
  struct ClassTraits<Herwig::Helicity::SMFFHVertex>
    : public ClassTraitsBase<Herwig::Helicity::SMFFHVertex> {

    /**
     * Return the class name.
     */
    static string className() { return "/Herwig++/Helicity/SMFFHVertex"; }

    /**
     * Return the name of the shared library to be loaded to get
     * access to this class and every other class it uses
     * (except the base class).
     */
    static string library() { return "libHwSMVertex.so"; }

  };
  
}


#endif /* HERWIG_SMFFHVertex_H */

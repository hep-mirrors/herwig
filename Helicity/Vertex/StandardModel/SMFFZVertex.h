// -*- C++ -*-
#ifndef HERWIG_SMFFZVertex_H
#define HERWIG_SMFFZVertex_H
//
// This is the declaration of the SMFFZVertex class.

#include "Herwig++/Helicity/Vertex/Vector/FFVVertex.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Utilities/Rebinder.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

/** \ingroup Helicity
 *
 *  This is the implementation of the Standard Model 
 *  fermion-antifermion Z vertex.
 *
 *  @see FFVVertex
 *  @see VertexBase
 */
class SMFFZVertex: public FFVVertex {

public:
  
  /**
   * Standard ctors and dtor.
   */
  inline SMFFZVertex();
  inline SMFFZVertex(const SMFFZVertex &);
  virtual ~SMFFZVertex();
  
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
  static ClassDescription<SMFFZVertex> initSMFFZVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  SMFFZVertex & operator=(const SMFFZVertex &);

private:

  /**
   * Pointer to the Standard Model object. 
   */
  SMPtr _theSM;

  /**
   * Storage of the couplings.
   */
  double _gl[17],_gr[17];
  Complex _couplast;
  Energy2 _q2last;

};

}
}

#include "SMFFZVertex.icc"
    
namespace ThePEG {
  
  /**
   * The following template specialization informs ThePEG about the
   * base class of SMFFZVertex.
   */
  template <>
  struct BaseClassTrait<Herwig::Helicity::SMFFZVertex,1> {
    typedef Herwig::Helicity::FFVVertex NthBase;
  };
  
  /**
   * The following template specialization informs ThePEG about the
   * name of this class and the shared object where it is defined.
   */
  template <>
  struct ClassTraits<Herwig::Helicity::SMFFZVertex>
    : public ClassTraitsBase<Herwig::Helicity::SMFFZVertex> {

    /**
     * Return the class name.
     */
    static string className() { return "/Herwig++/Helicity/SMFFZVertex"; }

    /**
     * Return the name of the shared library to be loaded to get
     * access to this class and every other class it uses
     *  (except the base class).
     */
    static string library() { return "libHwSMVertex.so"; }

  };
  
}

#endif /* HERWIG_SMFFZVertex_H */

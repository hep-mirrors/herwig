// -*- C++ -*-
#ifndef HERWIG_SMHPPVertex_H
#define HERWIG_SMHPPVertex_H
//
// This is the declaration of the SMHPPVertex class.
//

#include "Herwig++/Models/General/SVVLoopVertex.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "SMHPPVertex.fh"

namespace Herwig {
    
/**
 * The <code>SMHPPVertex</code> class implements the
 * setCoupling member for the Standard Model Higgs to  
 * gamma,gamma decay mode.
 */
class SMHPPVertex: public SVVLoopVertex {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline SMHPPVertex();
  //@}

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;
  
  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}
  
  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();
  
  /**
   * Calculate couplings
   *@param q2 Scale at which to evaluate coupling
   *@param part1 ParticleData pointer to first particle
   *@param part2 ParticleData pointer to second particle
   *@param part3 ParticleData pointer to third particle
   */
  virtual void setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
			   tcPDPtr part3);
  
protected:
  
  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const;
  
  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}
  
protected:

  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);
  
private:
  
  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<SMHPPVertex> initSMHPPVertex;
  
  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SMHPPVertex & operator=(const SMHPPVertex &);
  
  /**
   *Storage of couplings
   */
  //@{
  /**
   * Last value of the coupling calculated
   */
  Complex _couplast;
  
  /**
   * The scale \f$q^2\f$ at which coupling was last evaluated
   */
  Energy2 _q2last;
  //@}
  
  /**
   * Pointer to Standard Model object
   */
  tcHwSMPtr _theSM;
  
  /**
   * The mass of the \f$W\f$ boson.
   */
  Energy _mw;

  /**
   * Storage of \f$\sin\theta_W\f$
   */
  double _sw;
  
  /**
   * A pointer to the top quark ParticleData object 
   */
  tcPDPtr _top;

  /**
   * Whether we have calculated the tensor coefficients already
   */
  bool _haveCoeff;  
};
}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SMHPPVertex. */
template <>
struct BaseClassTrait<Herwig::SMHPPVertex,1> {
  /** Typedef of the first base class of SMHPPVertex. */
  typedef Herwig::SVVLoopVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SMHPPVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SMHPPVertex>
  : public ClassTraitsBase<Herwig::SMHPPVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SMHPPVertex"; }
};

/** @endcond */

}

#include "SMHPPVertex.icc"

#endif /* HERWIG_SMHPPVertex_H */

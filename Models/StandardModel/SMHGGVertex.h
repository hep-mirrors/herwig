// -*- C++ -*-
#ifndef HERWIG_SMHGGVertex_H
#define HERWIG_SMHGGVertex_H
//
// This is the declaration of the SMHGGVertex class.
//

#include "Herwig++/Helicity/Vertex/Scalar/SVVLoopVertex.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "SMHGGVertex.fh"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;
    
  /**
   * The <code>SMHGGVertex</code> class implements the
   * setCoupling member for the Standard Model Higgs to  
   * gluon, gluon decay mode.
   */
class SMHGGVertex: public SVVLoopVertex {
  
public:
  
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline SMHGGVertex();
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
   *@param part2 ParticleData pointer to first particle
   *@param part3 ParticleData pointer to first particle
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
  static ClassDescription<SMHGGVertex> initSMHGGVertex;
  
  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SMHGGVertex & operator=(const SMHGGVertex &);
  
  /**
   *Storage of couplings
   */
  //@{
  /**
   *Last value of the coupling calculated
   */
  Complex _couplast;
  
  /**
   *The scale \f$q^2\f$ at which coupling was last evaluated
   */
  Energy2 _q2last;
  //@}
  
  /**
   *Pointer to Standard Model object
   */
  tcHwSMPtr _theSM;
  
  /**
   *Mass of W boson for higgs coupling
   */
  Energy _mw;
  
  /**
   *Storage of \f$\sin\theta_W\f$
   */
  double _sw;
  
  /**
   * Option to turn on b in quark loop
   */
  int _qopt;
};

}
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SMHGGVertex. */
template <>
struct BaseClassTrait<Herwig::Helicity::SMHGGVertex,1> {
  /** Typedef of the first base class of SMHGGVertex. */
  typedef Herwig::Helicity::SVVLoopVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SMHGGVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::Helicity::SMHGGVertex>
  : public ClassTraitsBase<Herwig::Helicity::SMHGGVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::SMHGGVertex"; }
};

/** @endcond */

}

#include "SMHGGVertex.icc"

#endif /* HERWIG_SMHGGVertex_H */

// -*- C++ -*-
#ifndef HERWIG_QuarkoniumDecayer_H
#define HERWIG_QuarkoniumDecayer_H
//
// This is the declaration of the QuarkoniumDecayer class.
//

#include <ThePEG/Config/ThePEG.h>
#include <ThePEG/PDT/Decayer.h>
#include <ThePEG/Handlers/HandlerBase.h>
#include <ThePEG/Interface/Interfaced.h>
#include <ThePEG/PDT/DecayMode.h>
#include <ThePEG/Repository/Strategy.fh>
#include "DecayConfig.h"
#include "Herwig++/Utilities/GlobalParameters.h"
#include <fstream>

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 * <code>QuarkoniumDecayer</code> is a class that defines all the general routines 
 * used in HERWIG++ to imitate the HERWIG 6.4 decays. The goal is to have an exact
 * copy of HERWIG 6.4 decay routines. This will allow for easy 'callibration'
 * of the new C++ code with the old Fortran code.
 *
 *  This particular class is designed for the partonic decay of bottom and charmonium
 *  resonances. In general it is used for decays of the type:
 *    - \f$q,\bar{q}\f$ decay to a quark-antiquark pair generally using phase space,
 *      \e i.e. MECode=0.
 *
 *    - \f$g,g\f$ decay to two gluons normally using phase space,
 *      \e i.e. MECode=0.
 *
 *    - \f$g,g,g\f$ decay to three gluons, this will normally use the Ore-Powell 
 *      matrix element, \e i.e. MECode=130.
 *
 *    - \f$g,g,\gamma\f$ decay to two gluons and a photon, this will normally use
 *      the Ore-Powell matrix element, \e i.e. MECode=130.
 * 
 *
 *  This class supports two values of the MECode variable which can be set using
 *  the interface
 *
 *  - MECode=0   flat-phase space
 *  - MECode=130 The Ore-Powell onium matrix element.
 *
 *  The Ore-Powell matrix element is only for a three-body decay and this
 *  class only supports two and three body phase-space decays.
 *
 * @see HeavyDecayer
 * @see Hw64Decayer
 * @see Decayer
 *
 */
class QuarkoniumDecayer: public Decayer {

public:

  /**
   * Standard ctors and dtor
   */
  inline QuarkoniumDecayer();
  /**
   * Standard ctors and dtor
   */
  inline QuarkoniumDecayer(const QuarkoniumDecayer &);
  /**
   * Standard ctors and dtor
   */
  virtual ~QuarkoniumDecayer();

public:

  /**
   * return true if this decayer can perfom the decay specified by the
   * given decay mode.
   */
  virtual bool accept(const DecayMode &) const;

  /**
   * for a given decay mode and a given particle instance, perform the
   * decay and return the decay products.
   */
  virtual ParticleVector decay(const DecayMode &, const Particle &) const;

  /**
   * Standard Init function used to initialize the interface.
   */
  static void Init();

  /**
   * Standard Interfaced virtual functions
   */
  inline virtual void doupdate() throw(UpdateException);
  /**
   * Standard Interfaced virtual functions
   */
  inline virtual void doinit() throw(InitException);
  /**
   * Standard Interfaced virtual functions
   */
  inline virtual void dofinish();

  /**
   * Standard Persisten stream methods
   */
  void persistentOutput(PersistentOStream &) const;
  /**
   * Standard Persisten stream methods
   */
  void persistentInput(PersistentIStream &, int);

   /**
    * Standard clone methods
    */
protected:
   inline virtual IBPtr clone() const;
   /**
    * Standard clone methods
    */
   inline virtual IBPtr fullclone() const;

   /**
    * Change all pointers of Interfaced objects to corresponding clones
    */
   inline virtual void rebind(const TranslationMap &trans) 
                              throw(RebindException);

   /**
    * Return pointers to all interfaced objects referred to by this class.
    */
   inline virtual IVector getReferences();

private:

  /**
   *  The code for the type of matrix element being used.
   */
  int MECode;

  /**
   * A pointer to the global parameters object/
   */
  Ptr<GlobalParameters>::pointer global;

  /**
   * Describe a concrete class with persistant decay
   */
  static ClassDescription<QuarkoniumDecayer> initQuarkoniumDecayer;

  /**
   * Variable which control the adding of handlers, no longer used
   */
  static long lastAddedNumber;

  /**
   *  Private and non-existent assignment operator.
   */
  const QuarkoniumDecayer & operator=(const QuarkoniumDecayer &);
};

}

namespace ThePEG {

/**
 * This template specialization informs ThePEG about the base class of
 * QuarkoniumDecayer.
 */
template <>
struct BaseClassTrait<Herwig::QuarkoniumDecayer,1> {
  /** Typedef of the base class of QuarkoniumDecayer. */
  typedef HandlerBase NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * QuarkoniumDecayer class.
 */
template <>
struct ClassTraits<Herwig::QuarkoniumDecayer>: public ClassTraitsBase<Herwig::QuarkoniumDecayer> {
  /** Return the class name. */
  static string className() { return "Herwig++::QuarkoniumDecayer"; }
  /** Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwDecay.so"; }
};

}

#include "QuarkoniumDecayer.icc"

#endif /* HERWIG_QuarkoniumDecayer_H */

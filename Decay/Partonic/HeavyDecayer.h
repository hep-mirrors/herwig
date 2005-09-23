// -*- C++ -*-
#ifndef HERWIG_HeavyDecayer_H
#define HERWIG_HeavyDecayer_H
// This is the declaration of the HeavyDecayer class.

#include <ThePEG/Config/ThePEG.h>
#include <ThePEG/PDT/Decayer.h>
#include <ThePEG/Handlers/HandlerBase.h>
#include <ThePEG/Interface/Interfaced.h>
#include <ThePEG/PDT/DecayMode.h>
#include <ThePEG/Repository/Strategy.fh>
#include <fstream>

namespace Herwig {

using namespace ThePEG;

/** \ingroup Decay
 *
 *  <code>HeavyDecayer</code> is a class that defines all the general routines 
 *  used in HERWIG++ to imitate the HERWIG 6.4 decays. The goal is to have an exact
 *  copy of HERWIG 6.4 decay routines. This will allow for easy 'callibration'
 *  of the new C++ code with the old Fortran code.
 *
 *  This class is designed for the partonic decay of a bottom or charm mesons
 *  and baryons. In is
 *  only supports four body partonic decays.
 *
 *  Two types of matrix element are supported for this decay
 *
 *  - MECode=0   flat-phase space
 *  - MECode=100 V-A matrix element for the heavy quark decay in the spectator model.
 *
 * @see QuarkoniumDecayer
 * @see Hw64Decayer
 * @see Decayer
 * 
 */

class HeavyDecayer: public Decayer {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor
   */
  inline HeavyDecayer();

  /**
   * Copy constructor
   */
  inline HeavyDecayer(const HeavyDecayer &);

  /**
   * Destructor
   */
  virtual ~HeavyDecayer();
  //@}

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



  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Check sanity of the object during the setup phase.
   */
  inline virtual void doupdate() throw(UpdateException);

  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in
   * this object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
  //@}


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

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

private:

  /**
   *  The code for the matrix element being used.
   */
  int MECode;

  /**
   * Weighting of phase space for V-A matrix elements
   */
  static double VAWt(double*);

  /**
   *  Describe a concrete class with persistent data.
   */
  static ClassDescription<HeavyDecayer> initHeavyDecayer;

  /**
   *  This control adding handlers but doesn't seem to be used any more
   */
  static long lastAddedNumber;

  /**
   *  Private and non-existent assignment operator.
   */
  const HeavyDecayer & operator=(const HeavyDecayer &);
};

}

namespace ThePEG {

/**
 * This template specialization informs ThePEG about the base class of
 * HeavyDecayer.
 */
template <>
struct BaseClassTrait<Herwig::HeavyDecayer,1> {
  /** Typedef of the base class of HeavyDecayer. */
  typedef HandlerBase NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * HeavyDecayer class.
 */
template <>
struct ClassTraits<Herwig::HeavyDecayer>: public ClassTraitsBase<Herwig::HeavyDecayer> {
  /** Return the class name. */
  static string className() { return "Herwig++::HeavyDecayer"; }
  /** Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwPartonicDecay.so"; }
};

}

#include "HeavyDecayer.icc"

#endif /* HERWIG_HeavyDecayer_H */

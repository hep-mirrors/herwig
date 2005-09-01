// -*- C++ -*-
#ifndef HERWIG_WeakPartonicDecayer_H
#define HERWIG_WeakPartonicDecayer_H
//
// This is the declaration of the WeakPartonicDecayer class.
//

#include "ThePEG/PDT/Decayer.h"
#include "Herwig++/Utilities/GlobalParameters.h"
#include "WeakPartonicDecayer.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The WeakPartonicDecayer class is designed to replace the HeavyDecayer
 * class which implements the partonic decays of hadrons containing a 
 * heavy quark in the same way as in FORTRAN Herwig. 
 *
 * There are a number of major changes
 *
 * - The helicity formalism is used for the decays so that the \f$\tau\f$ lepton
 *   gets the correct correlations (in reality the main effect is to ensure the
 *   \f$\tau\f$ is left-handed.
 *
 * - The particles produced directly by the hadronisation, i.e. the primary hadrons
 *   produced in cluster decay are checked to ensure that none of the exclusive
 *   modes are reproduced.
 *
 * - Two body modes are allowed to try and force baryon production etc. In this case
 *   the colours of the partons are connected.
 * 
 * - Three body modes of the form \f$q g \bar{q}\f$ are supported for penguin mediated
 *   weak decays.
 *
 *  Two types of matrix element are supported for this decay
 *
 *  - MECode=0   flat-phase space.
 *  - MECode=100 V-A matrix element for the heavy quark decay in the spectator model.
 *
 * @see HeavyDecayer
 */
class WeakPartonicDecayer: public Decayer {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline WeakPartonicDecayer();

  /**
   * The copy constructor.
   */
  inline WeakPartonicDecayer(const WeakPartonicDecayer &);

  /**
   * The destructor.
   */
  virtual ~WeakPartonicDecayer();
  //@}

public:

  /** @name Virtual functions required by the Decayer class. */
  //@{
  /**
   * Check if this decayer can perfom the decay specified by the
   * given decay mode.
   * @param dm the DecayMode describing the decay.
   * @return true if this decayer can handle the given mode, otherwise false.
   */
  virtual bool accept(const DecayMode & dm) const;

  /**
   * Perform a decay for a given DecayMode and a given Particle instance.
   * @param dm the DecayMode describing the decay.
   * @param p the Particle instance to be decayed.
   * @return a ParticleVector containing the decay products.
   */
  virtual ParticleVector decay(const DecayMode & dm, const Particle & p) const;
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

public:

  /**
   * Weighting of phase space for V-A matrix elements
   */
  static double VAWt(double*);

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

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Check sanity of the object during the setup phase.
   */
  inline virtual void doupdate() throw(UpdateException);

  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  inline virtual void doinitrun();

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
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<WeakPartonicDecayer> initWeakPartonicDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  WeakPartonicDecayer & operator=(const WeakPartonicDecayer &);

private:

  /**
   *  The code for the matrix element being used.
   */
  int _MECode;

  /**
   * This is a pointer to a Herwig::GlobalParameters object for gluon mass
   */
  Ptr<GlobalParameters>::pointer _globalParameters;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** This template specialization informs ThePEG about the
 *  base classes of WeakPartonicDecayer. */
template <>
struct BaseClassTrait<Herwig::WeakPartonicDecayer,1> {
  /** Typedef of the first base class of WeakPartonicDecayer. */
  typedef Decayer NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the WeakPartonicDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::WeakPartonicDecayer>
  : public ClassTraitsBase<Herwig::WeakPartonicDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::WeakPartonicDecayer"; }
  /** Return the name of the shared library be loaded to get
   *  access to the WeakPartonicDecayer class and every other class it uses
   *  (except the base class). */
  static string library() { return "libHwPartonicDecay.so"; }
};

}

#include "WeakPartonicDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "WeakPartonicDecayer.tcc"
#endif

#endif /* HERWIG_WeakPartonicDecayer_H */
